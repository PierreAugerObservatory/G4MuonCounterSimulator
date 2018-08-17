/**
* @file PMTSimulator.cc
* @class PMTSimulator
* @brief Simulate the PMT signals (spe pulses, ...)
*
* Simulation of the PMT signals on the basis of the information produced in the G4 simulation
* @author S. Riggi
* @date 20/08/2010
*/

#include <PMTSimulator.hh>
#include <AnalysisConsts.hh>

#include <vector>
#include <map>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <cmath>

using namespace std;
using namespace G4MuonCounterSimulatorUSC;


TH1D* PMTSimulator::fTimeTraceHisto;
TH1D* PMTSimulator::fSignalTimeTraceHisto;
TH1D* PMTSimulator::fNoiseTimeTraceHisto;
TH1D* PMTSimulator::fPETimeDistributionHisto;
TH1D* PMTSimulator::fSPETimeTraceHisto;
std::vector<TH1D*> PMTSimulator::fSPETimeTrace;	


PMTSimulator::PMTSimulator(){

	fQE= 0.22;

	fFWHM= 3.0;//ns
  fRiseTime= 2.2;//ns
	fTransTime= 12.0;//25.0 ns
	fTransTimeSpread= 1.85;//3.0 ns
	fGain= 3.3E+06;
	fGainSpread= 10.;//in percentage %

	fPMTTimeTraceBinning= 0.1;//1.;//ns
	fPMTTimeTraceMin= 0.;//ns
  fPMTTimeTraceMax= 100.;//ns

	fPulseModel=kLandauTimePulse;
	fGainModel=kPolya;

	fDarkCurrent= 0.1;//nA

	fNoiseBaseline= 50.0E-05;//2.0E-05 A
	fNoiseRMS= 10.;//in %

  fQDCCounts=0.;
	fCharge2QDCCounts= 0.25E+12;// convert C to QDC channel

	fDiscriminatorThreshold= 5.0E-05;//A
	fArrivalTime=-999;
	
	fPECounts=0;
	fPhotonCounts=0;
	fParametrizationSmearingWidth= 0;

	fOutputFileName="TimeTraces.root";

}//close constructor


PMTSimulator::~PMTSimulator(){
/*
	if(fTimeTraceHisto!=NULL) fTimeTraceHisto->Delete();
	if(fSignalTimeTraceHisto!=NULL) fSignalTimeTraceHisto->Delete();
	if(fNoiseTimeTraceHisto!=NULL) fNoiseTimeTraceHisto->Delete();
	if(fPETimeDistributionHisto!=NULL) fPETimeDistributionHisto->Delete();	
	for(vector<TH1D*>::const_iterator it = fSPETimeTrace.begin(); it != fSPETimeTrace.end(); it++)
	{
		delete *it;
	}
	fSPETimeTrace.clear();
*/

	if(fTimeTraceHisto!=NULL) delete fTimeTraceHisto;
	if(fSignalTimeTraceHisto!=NULL) delete fSignalTimeTraceHisto;
	if(fNoiseTimeTraceHisto!=NULL) delete fNoiseTimeTraceHisto;
	if(fPETimeDistributionHisto!=NULL) delete fPETimeDistributionHisto;	
	for(vector<TH1D*>::const_iterator it = fSPETimeTrace.begin(); it != fSPETimeTrace.end(); it++)
	{
		delete *it;
	}
	fSPETimeTrace.clear();
	
}//close destructor


void PMTSimulator::Init(){

	//cout<<"## PMTSim::Init() ##"<<endl;

	fNTimeBins= (int)((fPMTTimeTraceMax-fPMTTimeTraceMin)/fPMTTimeTraceBinning);
	
	//initialize vectors
	for(int j=0;j<fNTimeBins;j++) {
		fSignalTrace.push_back(0.);//init with null amplitude
		fNoiseTrace.push_back(0.);//init with null amplitude
		fSPETraceContainer.push_back(0.);//init with null amplitude
    fNoiseTraceContainer.push_back(0.);//init with null amplitude
		fTimeBinCenter.push_back(fPMTTimeTraceMin+j*fPMTTimeTraceBinning/2.);
	}//close for j

	//filling QE table map
	for(int i=0;i<NUMENTRIES_QETable;i++){
		fQETable.insert( pair<double,double>(PhotonEnergy[i],QE[i]) );
	}
/*
	//## DEBUG 
	std::map<double,double>::iterator it;
	for(it=fQETable.begin();it!= fQETable.end();it++){
		double x= (*it).first;
		double y= (*it).second;
		cout<<"x="<<x<<"  y="<<y<<endl;
	}
*/

	fPETimeDistributionHisto= new TH1D("PETimeDistributionHisto","PETimeDistributionHisto",fNTimeBins,fPMTTimeTraceMin,fPMTTimeTraceMax);
	fPETimeDistributionHisto->GetXaxis()->SetTitle("t [ns]");
  fPETimeDistributionHisto->GetXaxis()->SetTitleSize(0.05);
	fPETimeDistributionHisto->GetXaxis()->SetTitleOffset(0.8);
	fPETimeDistributionHisto->GetYaxis()->SetTitle("entries");
	fPETimeDistributionHisto->GetYaxis()->SetTitleSize(0.05);

	fSignalTimeTraceHisto= new TH1D("SignalTimeTraceHisto","SignalTimeTraceHisto",fNTimeBins,fPMTTimeTraceMin,fPMTTimeTraceMax);
	fSignalTimeTraceHisto->GetXaxis()->SetTitle("t [ns]");
	fSignalTimeTraceHisto->GetXaxis()->SetTitleSize(0.05);
	fSignalTimeTraceHisto->GetXaxis()->SetTitleOffset(0.8);
	fSignalTimeTraceHisto->GetYaxis()->SetTitle("i(t) [A]");
	fSignalTimeTraceHisto->GetYaxis()->SetTitleSize(0.05);

	fNoiseTimeTraceHisto= new TH1D("NoiseTimeTraceHisto","NoiseTimeTraceHisto",fNTimeBins,fPMTTimeTraceMin,fPMTTimeTraceMax);		
	fNoiseTimeTraceHisto->GetXaxis()->SetTitle("t [ns]");
	fNoiseTimeTraceHisto->GetXaxis()->SetTitleSize(0.05);
	fNoiseTimeTraceHisto->GetXaxis()->SetTitleOffset(0.8);
	fNoiseTimeTraceHisto->GetYaxis()->SetTitle("i(t) [A]");
	fNoiseTimeTraceHisto->GetYaxis()->SetTitleSize(0.05);

	fTimeTraceHisto= new TH1D("TimeTraceHisto","TimeTraceHisto",fNTimeBins,fPMTTimeTraceMin,fPMTTimeTraceMax);
	fTimeTraceHisto->GetXaxis()->SetTitle("t [ns]");
	fTimeTraceHisto->GetXaxis()->SetTitleSize(0.05);
	fTimeTraceHisto->GetXaxis()->SetTitleOffset(0.8);
	fTimeTraceHisto->GetYaxis()->SetTitle("i(t) [A]");
	fTimeTraceHisto->GetYaxis()->SetTitleSize(0.05);

}//close function


bool PMTSimulator::GeneratePE(){

	//cout<<"## PMTSim::GeneratePE() ##"<<endl;
	//## Generate the photoelectrons produced at cathode, according to the PMT QE table
	
	if(fPhotonTime.empty()||fPhotonEnergy.empty()){
		//no photons at cathode...nothing to generate
		cout<<"## No photons at cathode...no photoelectrons to generate ##"<<endl;
		return true;
	}

	//clear PE vectors
	fPETime.clear();
	fPETime.resize(0);
	fPEEnergy.clear();
	fPEEnergy.resize(0);

	int Npe=0;
	int Nph= (int)(fPhotonTime.size());
	gRandom= new TRandom3();
	for(unsigned int s=0;s<fPhotonEnergy.size();s++)
	{
		double rand= gRandom->Uniform(0,1);
		double thisEnergy= fPhotonEnergy[s];
		double thisTime= fPhotonTime[s];
		double thisQE= GetQEFromTable(thisEnergy);
						
		if(rand<thisQE) {
			//create a pe, fill vectors
			Npe++;
			fPEEnergy.push_back(thisEnergy);
			fPETime.push_back(thisTime);
			fPETimeDistributionHisto->Fill(thisTime);
		}
		else continue;
	}//close for photon energy list
	
	cout<<"## "<< Npe <<" pe generated from "<< Nph << "  photons at cathode ##"<<endl;
	fPECounts= Npe;

	return true;

}//close PMTSim::GeneratePE()


bool PMTSimulator::GeneratePhotonsFast(double x,double Edep){

	//##to be filled!
	fPhotonCounts= 0;

	return true;

}//close PMTSim::GeneratePhotons()


bool PMTSimulator::GeneratePEFast(double x,double Edep){


	//### Fast PE Generation: parametrization data
	const int    fNVars  = 2;
	const int    fNCoeff = 9;
	const double fDMean  = 33.9867;

	double fVariables[fNVars];
	fVariables[0]= x;
	fVariables[1]= Edep;

	// Assignment to mean vector.
	const double fXMean[fNVars]= {-1.37668e-14, -0.789684};

	// Assignment to minimum vector.
	const double fXMin[fNVars]= {13.6773, 0.183674};

	// Assignment to maximum vector.
	const double fXMax[fNVars]= {106.423, 22.2794};

	// Assignment to coefficients vector.
	const double fCoefficient[fNCoeff] = {107.585, 136.655, -25.5017, -24.3021, 4.28797,
  																	-7.43121, -4.28537, 2.60583, 4.69381};

	// Assignment to error coefficients vector.
	const double fCoefficientRMS[fNCoeff] = {1e-10, 1e-10, 2.04783e-09,	2.04783e-09, 4.19361e-08, 
																			 4.19361e-08, 0.000360135, 1.75863e-05, 8.58779e-07};

	// Assignment to powers vector.
	// The powers are stored row-wise, that is
	//  p_ij = gPower[i * NVariables + j];
	const int    fPower[] = {1,  1,
  												 1,  2,
  												 2,  1,
  												 2,  2,
  												 3,  1,
  												 3,  3,
  												 6,  1,
  												 5,  1,
  												 4,  3};

	double parameterizedValue = fDMean;

  for(int i=0;i<fNCoeff;i++) {
    
		// Evaluate the ith term in the expansion
    double term = fCoefficient[i];
    for(int j=0;j<fNVars;j++){
      // Evaluate the polynomial in the jth variable.
      int power = fPower[fNVars * i + j]; 
      double p1 = 1, p2 = 0, p3 = 0, r = 0;
      double v =  1 + 2. / (fXMax[j] - fXMin[j]) * (fVariables[j] - fXMax[j]);
      
			// what is the power to use!
      switch(power) {
      case 1: r = 1; break; 
      case 2: r = v; break; 
      default: 
        p2 = v; 
        for(int k=3;k<=power;k++) { 
          p3 = p2 * v;
          p1 = p2; p2 = p3; 
        }
        r = p3;
      }
      // multiply this term by the poly in the jth var
      term *= r; 
    }
    // Add this term to the final result
    parameterizedValue += term;
  }

	bool acceptedPE= false;
	while(!acceptedPE){	
		int Npe= (int)(gRandom->Gaus(0,fParametrizationSmearingWidth)) + (int)(parameterizedValue);
		if(Npe>=0) acceptedPE= true;
	}

	fPECounts= (int)(gRandom->Gaus(0,fParametrizationSmearingWidth)) + (int)(parameterizedValue);

  return true;


}//close PMTSim::GeneratePhotons()


void PMTSimulator::PulseFinder(){

	//cout<<"## PMTSim::PulseFinder() ##"<<endl;

	//###########################################
	//######        GENERATE SIGNAL          ####
	//###########################################
	//## Generate superposition of spe pulses
	
	//loop over all photoelectrons
	for(unsigned int j=0;j<fPETime.size();j++){
		fSPETimeTraceHisto= new TH1D(Form("SPETimeTraceHisto_Npe%d",j+1),Form("SPETimeTraceHisto_Npe%d",j+1),fNTimeBins,fPMTTimeTraceMin,fPMTTimeTraceMax);
		fSPETimeTraceHisto->GetXaxis()->SetTitle("t [ns]");
		fSPETimeTraceHisto->GetXaxis()->SetTitleSize(0.05);
		fSPETimeTraceHisto->GetXaxis()->SetTitleOffset(0.8);
		fSPETimeTraceHisto->GetYaxis()->SetTitle("i(t) [A]");
		fSPETimeTraceHisto->GetYaxis()->SetTitleSize(0.05);	

		double thisTime= fPETime[j];

		//generate SPE pulse
  	SPEPulseSpectrum(thisTime,fPulseModel);
	
		//loop over trace bins
		for(unsigned int k=0;k<fSignalTrace.size();k++){
			fSignalTrace[k]+= fSPETraceContainer[k];// add SPE to total signal
			double bincontent= fSPETraceContainer[k]*ElectronCharge/(1.E-09);
			fSPETimeTraceHisto->SetBinContent(k+1,bincontent);
		}//close loop trace bins

		fSPETimeTrace.push_back(fSPETimeTraceHisto);
			
	}//close loop photoelectrons

	
	//###########################################
	//######        GENERATE NOISE           ####
	//###########################################
	//## generate noise
	NoisePulseSpectrum();
	for(unsigned int k=0;k<fNoiseTrace.size();k++){
		fNoiseTrace[k]= fNoiseTraceContainer[k];
	}//close loop trace bins

	
	//###########################################
	//######    SIGNAL+NOISE                 ####
	//###########################################
	for(int k=0;k<fSignalTimeTraceHisto->GetNbinsX();k++){
		double bincontent= fSignalTrace[k]*ElectronCharge/(1.E-09);
		fSignalTimeTraceHisto->SetBinContent(k+1,bincontent);

		bincontent= fNoiseTrace[k];
		fNoiseTimeTraceHisto->SetBinContent(k+1,bincontent);
	}//close for k
		
	fTimeTraceHisto->Add(fSignalTimeTraceHisto);
	fTimeTraceHisto->Add(fNoiseTimeTraceHisto);

/*
		//Calculate QDC Counts
		// *** integrate TimeTrace ==> Qint
		// *** 1) add nominal pedestal before conversion? ==> Qtot=Qint+Qped
		//				convert into QDC counts ==> Qtot*Charge2QDCCounts
		// *** 2) convert into QDC counts ==> QDCCounts=Qint*Charge2QDCCounts
		//        add measured pedestal   ==> QDCCounts=Qint*Charge2QDCCounts+Pedestal
		double Qint= TimeTraceHisto->Integral()*1.E-09;//in C
		fQDCPedestal= (1+0.03*fPMTTimeTraceBinning)*1.E-12;//in C
		double Qtot= Qint+ fQDCPedestal;
		double QDCCounts= Qtot* fCharge2QDCCounts;
		//cout<<"Qint="<<Qint<<" fQDCPedestal="<<fQDCPedestal<<" Qtot="<<Qtot<<" QDCCounts="<<QDCCounts<<endl;  

		if(i==0) fQDCCounts1= QDCCounts;
		if(i==1) fQDCCounts2= QDCCounts;

		//Calculate Arrival Time
		// *** scan bins and compare to threshold
		double ArrivalTime=-999;
		for(int k=0;k<TimeTraceHisto->GetNbinsX();k++){
			double bincontent=TimeTraceHisto->GetBinContent(k+1);	
			//cout<<"bin "<<k+1<<"bincenter"<<TimeTraceHisto->GetBinCenter(k+1)<<" bincontent="<<bincontent<<" threshold="<<fDiscriminatorThreshold<<endl;
			if(bincontent>=fDiscriminatorThreshold){
				ArrivalTime= TimeTraceHisto->GetBinLowEdge(k+1);
				break;
			}//close if
		}//close for k
		
		if(i==0) fArrivalTime1= ArrivalTime;
		if(i==1) fArrivalTime2= ArrivalTime;
*/


}//close PMTSim::PulseFinder()


void PMTSimulator::SPEPulseSpectrum(double t_cathode,int PulseModel){
	
	//cout<<"Calling PMTSim::SPESpectrum"<<endl;
	delete gRandom;
	gRandom= new TRandom3(0);


	double gain= GainGenerator(fGainModel);
	//double t_transit= gRandom->Gaus(fTransTime,fTransTimeSpread);
	double t_transit= fTransTime;
	double t_anode= t_cathode+t_transit;
		
	
	//#############################
	//##  Define pulse function
	//#############################
	TF1* SpeResponseFcn;

	if(PulseModel==kGaussTimePulse) {	
		SpeResponseFcn=new TF1("SpeResponseFcn","TMath::Gaus(x,[0],[1])*[2]",fPMTTimeTraceMin,fPMTTimeTraceMax);
		double mean= t_anode;
		double sigma = fFWHM/(2.*sqrt(2.*TMath::Log(2.)) ); 
		double norm= 1.;
		SpeResponseFcn->SetParameters(mean,sigma,norm);
	}
	else if(PulseModel==kTruncGaussTimePulse){
		SpeResponseFcn=new TF1("SpeResponseFcn","x/[0]*TMath::Exp(-x/[0])*[1]",fPMTTimeTraceMin,fPMTTimeTraceMax);
		double tau = fFWHM/(2.*sqrt(2.*TMath::Log(2.)) ); ; 
		double norm= 1.;
		SpeResponseFcn->SetParameters(tau,norm);
	}
	else if(PulseModel==kLandauTimePulse){
		SpeResponseFcn=new TF1("SpeResponseFcn","TMath::Landau(x,[0],[1])*[2]",fPMTTimeTraceMin,fPMTTimeTraceMax);	
		double mpshift= -0.22278298;// Landau maximum location
		double width_scale= 0.586;
		double width= width_scale*fFWHM/(2.*sqrt(2.*TMath::Log(2.)) ); 
		
		double mpc = t_anode - mpshift * width;
		double norm= 1.;
		SpeResponseFcn->SetParameters(mpc,width,norm);
	}
	else{
		cerr<<"Error in PulseModel assignement...exit"<<endl;
		exit(1);
	}

	//##############################
	//##  Normalize pulse function
	//##############################
  double t_bin;
	double pulse;
	double t_binLow;
	double t_binUp;
	double t_binCenter;
	double t_sum=0.;

	for(int i=0;i<fNTimeBins;i++){
		t_binLow= fPMTTimeTraceMin+i*fPMTTimeTraceBinning;
		t_binUp= fPMTTimeTraceMin+(i+1)*fPMTTimeTraceBinning;
		t_binCenter= t_binLow+0.5*fPMTTimeTraceBinning;
		pulse= SpeResponseFcn->Eval(t_binCenter);
		t_sum+= pulse*fPMTTimeTraceBinning;
		fSPETraceContainer[i]= pulse;
	}//close for i

	//normalization factor NORM= gain 
	for(int i=0;i<fNTimeBins;i++) fSPETraceContainer[i]*= gain/t_sum;
	
	if(SpeResponseFcn!=NULL) delete SpeResponseFcn;

}//close PMTSim::SPEPulseSpectrum()


void PMTSimulator::NoisePulseSpectrum(){

/*
	double TimeWindow= (fPMTTimeTraceMax-fPMTTimeTraceMin)*1.E-09;//in seconds
	double Bandwidth= 1./(2.*TimeWindow);//in Hz
	double DarkCurrentAtCathode= fDarkCurrent*1.E-09/fGain;//in Amp
	double AverageNpeAtCathode= fNpeAtCathode/((fPMTTimeTraceMax-fPMTTimeTraceMin)*1.E-09);//in Hz
	double SignalCurrentAtCathode= AverageNpeAtCathode*ElectronCharge;//in Amp
	
	double Noise_RMS= fGain*sqrt(2.*ElectronCharge*(DarkCurrentAtCathode+SignalCurrentAtCathode)*Bandwidth);//in Amp

	cout<<"Noise RMS="<<Noise_RMS<<endl;
	cout<<"Bandwidth="<<Bandwidth<<endl;
	cout<<"TimeWindow="<<TimeWindow<<endl;
  cout<<"AverageNpeCath="<<AverageNpeAtCathode<<endl;
	
	cout<<"DarkCurrentCath="<<DarkCurrentAtCathode<<endl;
	cout<<"SignalCurrentCath="<<SignalCurrentAtCathode<<endl;
*/

	double Noise_RMS= fNoiseRMS/100.*fNoiseBaseline;

	for(unsigned int k=0;k<fNoiseTraceContainer.size();k++){	
		fNoiseTraceContainer[k]= fabs( gRandom->Gaus(fNoiseBaseline,Noise_RMS) );					
	}//close for k

}//close PMTSim::NoisePulseSpectrum()


double PMTSimulator::GainGenerator(int GainModel){

	double gain=0.;

	unsigned long int seed1= time(0);
	unsigned long int seed2= clock();
  
	timespec timeStruct;
	int temp;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timeStruct);
	
	unsigned long int seed3= timeStruct.tv_nsec;


  delete gRandom;
	gRandom= new TRandom3(seed1+seed2+seed3);

  double gain_spread_par=0.;	

	if(GainModel==kGauss) {	
		gain_spread_par= fGainSpread*fGain/100.;	
		gain= gRandom->Gaus(fGain,gain_spread_par);
	}
	else if(GainModel==kPolya){
		gain_spread_par= 1./sqrt(fGainSpread/100.);

		TF1* GainDistrFcn= new TF1("GainDistrFcn",PMTSimulator::PolyaFcn,1.E+05,1.E+07,2);
  	GainDistrFcn->SetParameters(gain_spread_par,fGain);

		gain= GainDistrFcn->GetRandom();
		GainDistrFcn->Delete();
	}
	else if(GainModel==kFixed){
		gain= fGain;
	}
	else{
		cerr<<"Error in GainModel assignement...exit"<<endl;
		exit(1);
	}

	return gain;

}//close PMTSim::GainGenerator()


double PMTSimulator::TimeProp(double x){

	double Phi= gRandom->Uniform(-fFiberCriticalAngle,fFiberCriticalAngle);
	double cosPhi= cos(Phi);
	double propTime= x * fFiberCoreRefractiveIndex/(cosPhi*TMath::C());

	return propTime;

}//close PMTSimulator::TimeProp()

double PMTSimulator::TimeDecay(){

	TF1* decayFcn= new TF1("decayFcn"," (TMath::Exp(-x/[0])-TMath::Exp(-x/[1]))/([0]-[1]) ",0,100); 
	decayFcn->SetParameters(fScintillatorDecayTime,fFiberDecayTime);
	
	double decayTime= decayFcn->GetRandom();

	delete decayFcn;

	return decayTime;

}//close PMTSimulator::TimeDecay()

double PMTSimulator::TimeTransit(){

	double transitTime= gRandom->Gaus(fTransTime,fTransTimeSpread);

	return transitTime;

}//close PMTSimulator::TimeTransit()


double PMTSimulator::GeneratePETime(double x){

	double time= TimeDecay() + TimeProp(x) + TimeTransit();

	return time;

}//close PMTSimulator::TimeTransit()

double PMTSimulator::PolyaFcn(double* x,double* par){

	double b=par[0];// (rel err)^2= 1/b 
  double c=par[1];//mean

  double fcn= pow(b*x[0]/c,b-1)*TMath::Exp(-b*x[0]/c)*b/(c*TMath::Gamma(b));

  return fcn;

}//close PMTSim::PolyaFcn()

double PMTSimulator::GetQEFromTable(double E){

	double QEProbability;
	std::map<double,double>::iterator it,itlow,itup;

	//check if given argument exists as key in the map
	it= fQETable.find(E);
	if(it!= fQETable.end()) {
		//key found
		QEProbability= fQETable.find(E)->second;
		return QEProbability;
	}
	
	//if not, do a linear interpolation
  itlow= fQETable.lower_bound(E);  
  itup= fQETable.upper_bound(E);  
	
	//check ranges: if low-bound=end() or up_bound=begin() return 0
	if(itlow== fQETable.end()) {
		//Key outside upper-range...return 0
		return 0.;
	}
	else if(itup== fQETable.begin()) {
		//Key outside lower-range...return 0
		return 0.;
	}
	else{
		itlow--; 
		double QEProb_low= (*itlow).second;
		double QEProb_up= (*itup).second;
		QEProbability= 0.5*(QEProb_low+QEProb_up);
	}

	return QEProbability;

}//close PMTSim::GetQEFromTable()


void PMTSimulator::DoPMTSim(){

	Init();
	GeneratePE();
	PulseFinder();

}//close DoPMTSim()



