/**
* @file G4MuonCounterMCReader.cc
* @class G4MuonCounterMCReader
* @brief Read the output of the GEANT4 simulation
*
* Useful to read the output of the ROOT file produced in the simulation, get and draw sim info
* @author Dr. Simone Riggi, Dr.ssa Enrica Trovato
* @date 14/04/2010
*/

#include "G4MuonCounterMCReader.hh"
#include "TScintHit.hh"
#include "TPMTHit.hh"
#include "TrackPoint.hh"
#include "TrackReco.hh"
#include "AnalysisConsts.hh"
#include "G4UnitsTable.hh"
#include "PMTSimulator.hh"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;
using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterMCReader::G4MuonCounterMCReader(){

	fMCFileName= "SimDataOut.root";
	fOutputFileName= "Output.root";
  fEnergyThreshold= 0.;
	fPEThreshold= 0.;
	fPhotonThreshold= 0.;
	fGeneratePMTPulses= false;
	fGeneratePMTPhotonsFast= false;
	fGeneratePMTPEFast= false;
	fIsFastSimulation= false;
	fFastSimulationSmearing= 0;
	fNevents= 0;
	fVerbosity= 1;
}

G4MuonCounterMCReader::~G4MuonCounterMCReader(){

}


void G4MuonCounterMCReader::Init(){

 
	gROOT->ProcessLine(" gSystem->Load(\"libTHits.so\")");
  cout<<"Reading file "<<fMCFileName<<endl;
	fInputFile=new TFile(fMCFileName);
	if ( fInputFile->IsZombie() ) {
    cerr << " MCReader::ReadInfo() - error opening "
	       << fMCFileName << endl;
    throw std::runtime_error("error opening data file");
  }

	SimTree_pmt= (TTree*)fInputFile->Get("SimTree_pmt");
  SimTree_scint=(TTree*)fInputFile->Get("SimTree_scint");
  SimTree_gen=(TTree*)fInputFile->Get("SimTree_gen");
	SimTree_det=(TTree*)fInputFile->Get("SimTree_det");

	//## det tree
	SimTree_det->SetBranchAddress("Nstrips",&Nstrips);
	SimTree_det->SetBranchAddress("Nplanes",&Nplanes);
	SimTree_det->SetBranchAddress("XYPlaneDistance",&XYPlaneDistance);
	SimTree_det->SetBranchAddress("PlaneDistance",PlaneDistance);
	SimTree_det->SetBranchAddress("PlaneDepth",PlaneDepth);
	SimTree_det->SetBranchAddress("Xstrip",&Xstrip);
	SimTree_det->SetBranchAddress("Ystrip",&Ystrip);
	SimTree_det->SetBranchAddress("Zstrip",&Zstrip);
	SimTree_det->SetBranchAddress("Xcoat",&Xcoat);
	SimTree_det->SetBranchAddress("Ycoat",&Ycoat);
	SimTree_det->SetBranchAddress("Zcoat",&Zcoat);
	

	//## gen tree
  SimTree_gen->SetBranchAddress("X0",&X0);
  SimTree_gen->SetBranchAddress("Y0",&Y0);
  SimTree_gen->SetBranchAddress("Z0",&Z0);
	SimTree_gen->SetBranchAddress("Px",&Px);
  SimTree_gen->SetBranchAddress("Py",&Py);
  SimTree_gen->SetBranchAddress("Pz",&Pz);
	SimTree_gen->SetBranchAddress("KinEnergy",&KinEnergy);
	SimTree_gen->SetBranchAddress("Mass",&Mass);
	SimTree_gen->SetBranchAddress("Time",&Time);
  SimTree_gen->SetBranchAddress("PDGCode",&PDGCode);
  SimTree_gen->SetBranchAddress("Theta",&Theta);
	SimTree_gen->SetBranchAddress("Phi",&Phi);	

	//initialize vectors
	fScintHitCollection = 0;
	fPMTHitCollection = 0;  

	//## pmt tree
	SimTree_pmt->SetBranchAddress("PMTHitCollection",&fPMTHitCollection); 
	
	//## scintillator tree
	SimTree_scint->SetBranchAddress("ScintHitCollection",&fScintHitCollection);

	fOutputFile= new TFile(fOutputFileName,"RECREATE");
	genDir= fOutputFile->mkdir("generatedInfo","generatedInfo");
	recDir= fOutputFile->mkdir("reconstructedInfo","reconstructedInfo");
	recDir->cd();
	stripEdepDir= fOutputFile->mkdir("stripEdepInfo","stripEdepInfo");
	planeEdepDir= fOutputFile->mkdir("planeEdepInfo","planeEdepInfo");

	fNevents= SimTree_scint->GetEntries();

	//## init histo
	fGenEnergyHisto= new TH1D("Energy","Energy",100,0,10);
	fGenEnergyHisto->GetXaxis()->SetTitle("E [GeV]");
	fGenEnergyHisto->GetXaxis()->SetTitleSize(0.06);
	fGenEnergyHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenEnergyHisto->GetYaxis()->SetTitle("entries");
	fGenEnergyHisto->GetYaxis()->SetTitleSize(0.06);
	fGenEnergyHisto->GetYaxis()->SetTitleOffset(0.7);

	fGenMomentumXHisto= new TH1D("Px","Px",100,-10,10);
	fGenMomentumXHisto->GetXaxis()->SetTitle("Px [GeV]");
	fGenMomentumXHisto->GetXaxis()->SetTitleSize(0.06);
	fGenMomentumXHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenMomentumXHisto->GetYaxis()->SetTitle("entries");
	fGenMomentumXHisto->GetYaxis()->SetTitleSize(0.06);
	fGenMomentumXHisto->GetYaxis()->SetTitleOffset(0.7);
	
	fGenMomentumYHisto= new TH1D("Py","Py",100,-10,10);
	fGenMomentumYHisto->GetXaxis()->SetTitle("Py [GeV]");
	fGenMomentumYHisto->GetXaxis()->SetTitleSize(0.06);
	fGenMomentumYHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenMomentumYHisto->GetYaxis()->SetTitle("entries");
	fGenMomentumYHisto->GetYaxis()->SetTitleSize(0.06);
	fGenMomentumYHisto->GetYaxis()->SetTitleOffset(0.7);

	fGenMomentumZHisto= new TH1D("Pz","Pz",100,-10,10);
	fGenMomentumZHisto->GetXaxis()->SetTitle("Pz [GeV]");
	fGenMomentumZHisto->GetXaxis()->SetTitleSize(0.06);
	fGenMomentumZHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenMomentumZHisto->GetYaxis()->SetTitle("entries");
	fGenMomentumZHisto->GetYaxis()->SetTitleSize(0.06);
	fGenMomentumZHisto->GetYaxis()->SetTitleOffset(0.7);
 
  fGenVertexXHisto= new TH1D("X","X",100,-500,500); 
	fGenVertexXHisto->GetXaxis()->SetTitle("X [cm]");
	fGenVertexXHisto->GetXaxis()->SetTitleSize(0.06);
	fGenVertexXHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenVertexXHisto->GetYaxis()->SetTitle("entries");
	fGenVertexXHisto->GetYaxis()->SetTitleSize(0.06);
	fGenVertexXHisto->GetYaxis()->SetTitleOffset(0.7);

  fGenVertexYHisto= new TH1D("Y","Y",100,-500,500); 
	fGenVertexYHisto->GetXaxis()->SetTitle("Y [cm]");
	fGenVertexYHisto->GetXaxis()->SetTitleSize(0.06);
	fGenVertexYHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenVertexYHisto->GetYaxis()->SetTitle("entries");
	fGenVertexYHisto->GetYaxis()->SetTitleSize(0.06);
	fGenVertexYHisto->GetYaxis()->SetTitleOffset(0.7);

  fGenVertexZHisto= new TH1D("Z","Z",100,-500,500);
	fGenVertexZHisto->GetXaxis()->SetTitle("Z [cm]");
	fGenVertexZHisto->GetXaxis()->SetTitleSize(0.06);
	fGenVertexZHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenVertexZHisto->GetYaxis()->SetTitle("entries");
	fGenVertexZHisto->GetYaxis()->SetTitleSize(0.06);
	fGenVertexZHisto->GetYaxis()->SetTitleOffset(0.7);
 
	fGenThetaHisto= new TH1D("Theta","Theta",100,0,360); 
	fGenThetaHisto->GetXaxis()->SetTitle("Theta [deg]");
	fGenThetaHisto->GetXaxis()->SetTitleSize(0.06);
	fGenThetaHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenThetaHisto->GetYaxis()->SetTitle("entries");
	fGenThetaHisto->GetYaxis()->SetTitleSize(0.06);
	fGenThetaHisto->GetYaxis()->SetTitleOffset(0.7);

  fGenPhiHisto= new TH1D("Phi","Phi",100,0,360); 
	fGenPhiHisto->GetXaxis()->SetTitle("Phi [deg]");
	fGenPhiHisto->GetXaxis()->SetTitleSize(0.06);
	fGenPhiHisto->GetXaxis()->SetTitleOffset(0.7);
	fGenPhiHisto->GetYaxis()->SetTitle("entries");
	fGenPhiHisto->GetYaxis()->SetTitleSize(0.06);
	fGenPhiHisto->GetYaxis()->SetTitleOffset(0.7);

	fPhotonCountsMap= new TH2D("PhotonCountsMap","PhotonCountsMap",100,0,120,100,0,10); 	
  fPhotonCountsMap->Sumw2();
	fPhotonCountsMap->GetXaxis()->SetTitle("Hit-PMT Distance [cm]");
	fPhotonCountsMap->GetXaxis()->SetTitleSize(0.06);
	fPhotonCountsMap->GetXaxis()->SetTitleOffset(0.7);
	fPhotonCountsMap->GetYaxis()->SetTitle("Energy deposit [MeV]");
	fPhotonCountsMap->GetYaxis()->SetTitleSize(0.06);
	fPhotonCountsMap->GetYaxis()->SetTitleOffset(0.7);
	fPhotonCountsMap->GetZaxis()->SetTitle("photon Counts");
	fPhotonCountsMap->GetZaxis()->SetTitleSize(0.06);
	fPhotonCountsMap->GetZaxis()->SetTitleOffset(0.7);

	fPECountsMap= new TH2D("PECountsMap","PECountsMap",100,0,120,100,0,10); 	
  fPECountsMap->Sumw2();
	fPECountsMap->GetXaxis()->SetTitle("Hit-PMT Distance [cm]");
	fPECountsMap->GetXaxis()->SetTitleSize(0.06);
	fPECountsMap->GetXaxis()->SetTitleOffset(0.7);
	fPECountsMap->GetYaxis()->SetTitle("Energy deposit [MeV]");
	fPECountsMap->GetYaxis()->SetTitleSize(0.06);
	fPECountsMap->GetYaxis()->SetTitleOffset(0.7);
	fPECountsMap->GetZaxis()->SetTitle("photon Counts");
	fPECountsMap->GetZaxis()->SetTitleSize(0.06);
	fPECountsMap->GetZaxis()->SetTitleOffset(0.7);

	fCountsMapTree= new TTree("CountsMapTree","CountsMapTree");
	fCountsMapTree->Branch("x",&x,"x/D");
	fCountsMapTree->Branch("Edep",&Edep,"Edep/D");
	fCountsMapTree->Branch("phCounts",&phCounts,"phCounts/I");
	fCountsMapTree->Branch("peCounts",&peCounts,"peCounts/I");

	fTimeMapTree= new TTree("TimeMapTree","TimeMapTree");
	fTimeMapTree->Branch("d",&d,"d/D");
	fTimeMapTree->Branch("t_hit",&t_hit,"t_hit/D");
	fTimeMapTree->Branch("t",&t,"t/D");
	

}//close MCReader::Init()


void G4MuonCounterMCReader::Reset(){
    
}//close MCReader::Reset()

void G4MuonCounterMCReader::Fill(){

	for(int i=0;i<fNevents;i++){
		cout<<"*** EVENT no "<<i+1<<"***"<<endl;
		ReadDetectorInfo(i);
		ReadGeneratedInfo(i);
		ReadScintillatorInfo(i);
		if(!fIsFastSimulation) ReadPMTInfo(i);
		cout<<endl;
		//cout<<endl;
	}//close loop events

}//close MCReader::Fill()

void G4MuonCounterMCReader::ReadDetectorInfo(int thisEvent){

	//## print DET INFO	
	SimTree_det->GetEntry(thisEvent);		
	
}//close MCReader::ReadDetectorInfo()


void G4MuonCounterMCReader::ReadGeneratedInfo(int thisEvent){


		SimTree_gen->GetEntry(thisEvent);

		if(fVerbosity>1) 
		{
			cout<<"  X0[cm]="<<X0/cm<<"  Y0[cm]="<<Y0/cm<<"  Z0[cm]="<<Z0/cm<<"  Theta[deg]="<<Theta/deg<<"  Px[GeV]="<<Px/GeV<<"  Py[GeV]="<<Py/GeV<<"  Pz[GeV]="<<Pz/GeV<<"  KinEnergy[GeV]="<<KinEnergy/GeV<<"  Mass[MeV]="<<Mass/MeV<<"  Time[ns]="<<Time<<"  PDGCode="<<PDGCode<<endl;
		}
	
		//Fill histogram with current event	
		fGenEnergyHisto->Fill(KinEnergy/GeV);
		fGenMomentumXHisto->Fill(Px/GeV);
		fGenMomentumYHisto->Fill(Py/GeV);
		fGenMomentumZHisto->Fill(Pz/GeV);
		fGenVertexXHisto->Fill(X0/cm);
		fGenVertexYHisto->Fill(Y0/cm);
		fGenVertexZHisto->Fill(Z0/cm);
		fGenThetaHisto->Fill(Theta/deg);
		fGenPhiHisto->Fill(Phi/deg);

		//calculate geometrical hits for the current generated particle
		//std::vector<int> GeomHitStripId= CalculateGeomHits();
		//fGeomHitStripId.push_back(GeomHitStripId);


}//close MCReader::ReadGeneratedInfo()


void G4MuonCounterMCReader::ReadScintillatorInfo(int thisEvent){
    

		SimTree_scint->GetEntry(thisEvent);
		
		//container vector for plane and strip ID, used to order hits by planeID 
		fHitStripList.clear();
		fHitStripList.resize(0);

		std::vector<int> fPlaneList;
		fPlaneList.clear();
		fPlaneList.resize(0);


		if(fVerbosity>1) cout<< fScintHitCollection->size() << " hit strips in this event" << endl;

		//## loop over all strips hit in this event
		for(unsigned int j=0;j<fScintHitCollection->size();j++){
			int StripId= (fScintHitCollection->at(j)).StripId;
			int PlaneId= (fScintHitCollection->at(j)).PlaneId;
			int SuperPlaneId= (fScintHitCollection->at(j)).SuperPlaneId;
			TVector3 StripPos= (fScintHitCollection->at(j)).StripPosition;
			int AbsPlaneId= 2*SuperPlaneId+PlaneId;
			int AbsStripId= AbsPlaneId*Nstrips+StripId;
			double Etot= (fScintHitCollection->at(j)).Etot;
			std::vector<TVector3> PhotocathodeSurfPos= (fScintHitCollection->at(j)).PhotocathodeSurfacePosition;

			unsigned int Nhits= (fScintHitCollection->at(j)).Edep.size();
			char PlaneType[1];
      if(PlaneId==0) PlaneType[0]='x';
      else PlaneType[0]='y';

			if(fVerbosity>1) cout<< Nhits << " hits for the current strip" << endl;

			TVector3 WeightPos(0,0,0);
			TVector3 currentPos(0,0,0);
			double currentEdep=0.;
			double currentTime=0.;
			double WeightTime=0.;
			for(unsigned int k=0;k<Nhits;k++){	
				currentPos= (fScintHitCollection->at(j)).Position[k];
				currentEdep= (fScintHitCollection->at(j)).Edep[k];
				currentTime= (fScintHitCollection->at(j)).Time[k];
				WeightPos+= currentPos*currentEdep;	
				WeightTime+= currentTime*currentEdep;	
				if(fVerbosity>2) 
				{
					cout<<"hit "<<k+1<<"  Edep="<< currentEdep <<"  Pos=("<< currentPos.X() <<","<< currentPos.Y() <<","<< currentPos.Z()<<")  time="<<currentTime<<endl;
				}
				
			}//close loop vectors
			WeightPos.SetXYZ(WeightPos.X()/Etot,WeightPos.Y()/Etot,WeightPos.Z()/Etot);
			WeightTime/= Etot;  

			//calculate HIT-PMT distance
			std::vector<TVector3> HitPMTDistance;
			HitPMTDistance.clear();
			HitPMTDistance.resize(0);
			std::vector<double> HitPMTHorizontalDistance;
			HitPMTHorizontalDistance.clear();
			HitPMTHorizontalDistance.resize(0);
			std::vector<int> PMTPhotonCounts;
			PMTPhotonCounts.clear();
			PMTPhotonCounts.resize(0);
			std::vector<int> PMTPECounts;
			PMTPECounts.clear();
			PMTPECounts.resize(0);
			

			for(unsigned int l=0;l<PhotocathodeSurfPos.size();l++){
				HitPMTDistance.push_back(PhotocathodeSurfPos[l]-WeightPos);
		
				if(AbsPlaneId%2==0){//plane x
					//take distance in x
					HitPMTHorizontalDistance.push_back( fabs(HitPMTDistance[l].X()) );
				}
				else if(AbsPlaneId%2==1){//plane y
					//take distance in y
					HitPMTHorizontalDistance.push_back( fabs(HitPMTDistance[l].Y()) );
				}
				
				if(fVerbosity>2) 
				{
					cout<<"PhotocathodeSurfPos=("<<PhotocathodeSurfPos[l].X()<<","<<PhotocathodeSurfPos[l].Y()<<","<<PhotocathodeSurfPos[l].Z()<<"  HitPMTDistance=("<<HitPMTDistance[l].X()<<","<<HitPMTDistance[l].Y()<<","<<HitPMTDistance[l].Z()<<"  HitPMTHorizontalDistance="<<HitPMTHorizontalDistance[l]<<endl;
				}

			}//close loop hit-pmt distance


			if(fVerbosity>2) 
			{
				cout<<"SuperPlaneId="<<SuperPlaneId+1<<"  AbsPlaneId="<<AbsPlaneId+1<<"  StripId="<<StripId+1<<"  Etot="<<Etot<<"  StripPos=("<<StripPos.X()<<","<<StripPos.Y()<<","<<StripPos.Z()<<")"<<"  StripHitPos=("<<WeightPos.X()<<","<<WeightPos.Y()<<","<<WeightPos.Z()<<")  HitTime="<<WeightTime<<endl;
			}

			
			//## select hit strips according to total energy deposit and chosen energy threshold
			if(Etot < fEnergyThreshold) continue; //skip this hit strip

			//## do a fast PMT simulation (Nphotons, Npe, timing)
			if(fIsFastSimulation){
				for(unsigned int l=0;l<HitPMTHorizontalDistance.size();l++){
					PMTSimulator* fPMTSim= new PMTSimulator();
					fPMTSim->SetParametrizationSmearingWidth(fFastSimulationSmearing);

					fPMTSim->Init();
						
					//## generate photons from parametrization, given x & Edep
					fPMTSim->GeneratePhotonsFast(HitPMTHorizontalDistance[l],Etot);
				
					//## generate photoelectrons from parametrization, given x & Edep
					fPMTSim->GeneratePEFast(HitPMTHorizontalDistance[l],Etot);

					int PECounts= fPMTSim->GetPECounts();
					int PhotonCounts= fPMTSim->GetPhotonCounts();
					if(fVerbosity>2) 
					{
						cout<<"  SuperPlaneId="<<SuperPlaneId+1<<"  AbsPlaneId="<<AbsPlaneId+1<<"  StripId="<<StripId+1<<"  PMTId "<<l<<"  Etot="<<Etot<<"  x="<<HitPMTHorizontalDistance[l]<<"  Npe="<<PECounts<<endl;
					}

					PMTPECounts.push_back(PECounts);
					PMTPhotonCounts.push_back(PhotonCounts);
			
					delete fPMTSim;
				}//close loop PMT
			}//close if Fast Simulation


			//## Create struct for vector ordering	
			HitStripStruct hitStripStruct;
			hitStripStruct.stripID= StripId;
			hitStripStruct.planeID= PlaneId;
			hitStripStruct.superplaneID= SuperPlaneId;
			hitStripStruct.absplaneID= AbsPlaneId;
			hitStripStruct.absstripID= AbsStripId;
			hitStripStruct.etot= Etot;
			hitStripStruct.time= WeightTime;
			hitStripStruct.stripPosition= StripPos;
			hitStripStruct.stripHitPosition= WeightPos;
			hitStripStruct.stripPMTPosition= PhotocathodeSurfPos;
			hitStripStruct.stripHitPMTDistance= HitPMTHorizontalDistance;
			if(fIsFastSimulation){
				hitStripStruct.stripPMTPECounts= PMTPECounts;
				hitStripStruct.stripPMTPhotonCounts= PMTPhotonCounts;
			}
			fHitStripList.push_back(hitStripStruct);


			//## Create and fill histograms
			bool isHistoAllocated= false;
			int HistoEntryId;
			
			for(unsigned int l=0;l<fSelectedStripList.size();l++){//loop over all hits in current event
				if(AbsStripId==fSelectedStripList[l]){
					isHistoAllocated=true; //this strip module has been already allocated as histograms
					HistoEntryId= l;
					break;
				}
			}//close loop

			if(!isHistoAllocated){
				//allocate and fill histograms
				TString histoName= Form("Edep_superplane%d_plane%d_strip%d",SuperPlaneId,PlaneId,StripId);
				fEnergyDepositHisto= new TH1D(histoName,histoName,100,0,10);
				fEnergyDepositHisto->GetXaxis()->SetTitle("Energy Deposit [MeV]");
				fEnergyDepositHisto->GetXaxis()->SetTitleSize(0.06);
				fEnergyDepositHisto->GetXaxis()->SetTitleOffset(0.7);
				fEnergyDepositHisto->GetYaxis()->SetTitle("entries");
				fEnergyDepositHisto->GetYaxis()->SetTitleSize(0.06);
				fEnergyDepositHisto->GetYaxis()->SetTitleOffset(0.7);
				fEnergyDepositHisto->Fill(Etot/MeV);

				fSelectedStripList.push_back(AbsStripId);
				fEnergyDepositHistoList.push_back(fEnergyDepositHisto);
			}
			else{
				//retrieve and fill existing histograms
				fEnergyDepositHistoList[HistoEntryId]->Fill(Etot/MeV);
			}

			isHistoAllocated= false;
			HistoEntryId= -1;
			
			for(unsigned int l=0;l<fSelectedPlaneList.size();l++){//loop over all hits in current event
				if(AbsPlaneId==fSelectedPlaneList[l]){
					isHistoAllocated=true; //this strip module has been already allocated as histograms
					HistoEntryId= l;
					break;
				}
			}//close loop

			if(!isHistoAllocated){
				//allocate and fill histograms
				TString histoName= Form("Edep_superplane%d_plane%d",SuperPlaneId,PlaneId);
				fEnergyDepositInPlaneHisto= new TH1D(histoName,histoName,100,0,10);
				fEnergyDepositInPlaneHisto->GetXaxis()->SetTitle("Energy Deposit [MeV]");
				fEnergyDepositInPlaneHisto->GetXaxis()->SetTitleSize(0.06);
				fEnergyDepositInPlaneHisto->GetXaxis()->SetTitleOffset(0.7);
				fEnergyDepositInPlaneHisto->GetYaxis()->SetTitle("entries");
				fEnergyDepositInPlaneHisto->GetYaxis()->SetTitleSize(0.06);
				fEnergyDepositInPlaneHisto->GetYaxis()->SetTitleOffset(0.7);
				if(PlaneId%2==0) fEnergyDepositInPlaneHisto->SetFillColor(kGreen);
				else fEnergyDepositInPlaneHisto->SetFillColor(kRed);

				fEnergyDepositInPlaneHisto->Fill(Etot/MeV);
				

				fSelectedPlaneList.push_back(AbsPlaneId);
				fEnergyDepositInPlaneHistoList.push_back(fEnergyDepositInPlaneHisto);

				//allocate and fill histograms
				if(PlaneId%2==0) histoName= Form("Multiplicity_Plane%dX",SuperPlaneId+1);	
				else histoName= Form("Multiplicity_Plane%dY",SuperPlaneId+1);	
				fStripMultiplicityInPlaneHisto= new TH1D(histoName,histoName,21,-0.5,20.5);
				fStripMultiplicityInPlaneHisto->GetXaxis()->SetTitle("Multiplicity");
				fStripMultiplicityInPlaneHisto->GetXaxis()->SetTitleSize(0.06);
				fStripMultiplicityInPlaneHisto->GetXaxis()->SetTitleOffset(0.7);
				fStripMultiplicityInPlaneHisto->GetYaxis()->SetTitle("entries");
				fStripMultiplicityInPlaneHisto->GetYaxis()->SetTitleSize(0.06);
				fStripMultiplicityInPlaneHisto->GetYaxis()->SetTitleOffset(0.7);
				if(PlaneId%2==0) fStripMultiplicityInPlaneHisto->SetFillColor(kGreen);
				else fStripMultiplicityInPlaneHisto->SetFillColor(kRed);

				fStripMultiplicityInPlaneHistoList.push_back(fStripMultiplicityInPlaneHisto);

			}
			else{
				//retrieve and fill existing histograms
				fEnergyDepositInPlaneHistoList[HistoEntryId]->Fill(Etot/MeV);
			}

			
			//## store PMT histo
			if(fIsFastSimulation){
				for(unsigned int l=0;l<HitPMTHorizontalDistance.size();l++){
					isHistoAllocated= false;
					HistoEntryId= -1;
					int PMTId= l;

					for(unsigned int k=0;k<fSelectedPMTStripList.size();k++){//loop over all hits in current event
						if(AbsStripId==fSelectedPMTStripList[k] && PMTId==fSelectedPMTList[k]){
							isHistoAllocated=true; //this strip module has been already allocated as histograms
							HistoEntryId= k;
							break;
						}
					}//close loop

					if(!isHistoAllocated){
						//allocate and fill histograms
						TString histoName= Form("PMTCounts_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
						fPhotonCountsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
						fPhotonCountsHisto->Sumw2();
						fPhotonCountsHisto->Fill(HitPMTHorizontalDistance[l], (double)(PMTPhotonCounts[l]));

						histoName= Form("PMTAverageCounts_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
						fAveragePhotonCountsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
						fAveragePhotonCountsHisto->Sumw2();
						fAveragePhotonCountsHisto->Fill(HitPMTHorizontalDistance[l], (double)(PMTPhotonCounts[l]));

						histoName= Form("PECounts_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
						fPECountsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
						fPECountsHisto->Sumw2();
						fPECountsHisto->Fill(HitPMTHorizontalDistance[l], (double)(PMTPECounts[l]));

						histoName= Form("PEAverageCounts_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
						fAveragePECountsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
						fAveragePECountsHisto->Sumw2();
						fAveragePECountsHisto->Fill(HitPMTHorizontalDistance[l], (double)(PMTPECounts[l]));

						histoName= Form("NumberOfEvents_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
						fNumberOfEventsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
						fNumberOfEventsHisto->Sumw2();
						fNumberOfEventsHisto->Fill(HitPMTHorizontalDistance[l], 1.);

						histoName= Form("PhotonEfficiency_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
						fPhotonEfficiencyHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
						fPhotonEfficiencyHisto->Sumw2();
						if(PMTPhotonCounts[l] >= fPhotonThreshold) fPhotonEfficiencyHisto->Fill(HitPMTHorizontalDistance[l], 1.);
	
						histoName= Form("PEEfficiency_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
						fPEEfficiencyHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
						fPEEfficiencyHisto->Sumw2();
						if(PMTPECounts[l] >= fPEThreshold) fPEEfficiencyHisto->Fill(HitPMTHorizontalDistance[l], 1.);

						fPhotonCountsHistoList.push_back(fPhotonCountsHisto);
						fAveragePhotonCountsHistoList.push_back(fAveragePhotonCountsHisto);
						fPECountsHistoList.push_back(fPECountsHisto);
						fAveragePECountsHistoList.push_back(fAveragePECountsHisto);
						fNumberOfEventsHistoList.push_back(fNumberOfEventsHisto);
						fPhotonEfficiencyHistoList.push_back(fPhotonEfficiencyHisto);
						fPEEfficiencyHistoList.push_back(fPEEfficiencyHisto);

						fSelectedPMTStripList.push_back(AbsStripId);
						fSelectedPMTList.push_back(PMTId);				
					}
					else{
						//retrieve and fill existing histograms
						fPhotonCountsHistoList[HistoEntryId]->Fill(HitPMTHorizontalDistance[l], (double)(PMTPhotonCounts[l]));
						fAveragePhotonCountsHistoList[HistoEntryId]->Fill(HitPMTHorizontalDistance[l], (double)(PMTPhotonCounts[l]));
						fPECountsHistoList[HistoEntryId]->Fill(HitPMTHorizontalDistance[l], (double)(PMTPECounts[l]));
						fAveragePECountsHistoList[HistoEntryId]->Fill(HitPMTHorizontalDistance[l], (double)(PMTPECounts[l]));
						fNumberOfEventsHistoList[HistoEntryId]->Fill(HitPMTHorizontalDistance[l], 1.);
						if(PMTPhotonCounts[l] >= fPhotonThreshold) fPhotonEfficiencyHistoList[HistoEntryId]->Fill(HitPMTHorizontalDistance[l], 1.);
						if(PMTPECounts[l] >= fPEThreshold) fPEEfficiencyHistoList[HistoEntryId]->Fill(HitPMTHorizontalDistance[l], 1.);
					}

					x= HitPMTHorizontalDistance[l];
					Edep= Etot;
					phCounts= PMTPhotonCounts[l];
					peCounts= PMTPECounts[l];
					fCountsMapTree->Fill();

				}//close loop PMTs
			}//close if Fast Simulation

		}//close loop j hit strips in this event


		//## Order the hit strip by super plane Id, plane Id, stripId using the absolute stripID
		std::sort(fHitStripList.begin(), fHitStripList.end(), order_by_absstripID());
		
		//fill list of hit planes: MUST BE ORDERED!
		for(unsigned int k=0; k<fHitStripList.size();k++){
			if(fVerbosity>1){
				cout<<"SuperPlaneId="<<fHitStripList[k].superplaneID+1<<"  AbsPlaneId="<<fHitStripList[k].absplaneID+1<<"  StripId="<<fHitStripList[k].stripID+1<<"  StripPos=("<<(fHitStripList[k].stripPosition).X()<<","<<(fHitStripList[k].stripPosition).Y()<<","<<(fHitStripList[k].stripPosition).Z()<<")"<<endl;
			} 

			int AbsPlaneId= fHitStripList[k].absplaneID;
			bool isPlaneInserted= false;
			for(unsigned int l=0; l<fPlaneList.size();l++){
				if(AbsPlaneId==fPlaneList[l]) {
					isPlaneInserted= true;
					break;
				}
			}
			if(!isPlaneInserted) fPlaneList.push_back(AbsPlaneId); 
		}//close for k
			

		//## Calculate strip multiplicity per plane
		for(unsigned int k=0;k<fPlaneList.size();k++){//loop over all plane hits in current event
			int AbsPlaneId= fPlaneList[k];
			int Multiplicity= GetStripMultiplicity(AbsPlaneId);

			int HistoEntryId= -1;
	
			for(unsigned int l=0;l<fSelectedPlaneList.size();l++){//loop over all plane hits in events
				if(AbsPlaneId==fSelectedPlaneList[l]){
					HistoEntryId= l;
					break;
				}
			}//close loop

			//retrieve and fill existing histograms
			fStripMultiplicityInPlaneHistoList[HistoEntryId]->Fill(Multiplicity);
			
		}//close loop k over all plane hits in current event

		//## create track points and add to track point collection
		CreateTrackPoints();
		cout<<"## TRACK POINTS ##"<<endl;
		for(unsigned int k=0;k<fTrackPointCollection.size();k++){
			cout<<"DET PLANE="<<fTrackPointCollection[k].fDetectorPlaneId<<"  Pos=("<<(fTrackPointCollection[k].fPosition).X()<<","<<(fTrackPointCollection[k].fPosition).Y()<<","<<(fTrackPointCollection[k].fPosition).Z()<<")"<<endl;
		}
	
		//## try track finder
		TrackReco aTrackReco;
		aTrackReco.SetTrackPoints(fTrackPointCollection);
		aTrackReco.TrackReconstruction();


		if(fVerbosity>2)
		{
			for(unsigned int k=0;k<fHitStripList.size();k++){
				cout<<"AbsPlaneId="<<fHitStripList[k].absplaneID+1<<"  SuperPlaneId=  "<<fHitStripList[k].superplaneID+1<<"   PlaneId=  "<<fHitStripList[k].planeID<<"  StripId="<<fHitStripList[k].stripID+1<<"  AbsStripId="<<fHitStripList[k].absstripID+1<<" StripPos="<<fHitStripList[k].stripPosition.X()<<","<<fHitStripList[k].stripPosition.Y()<<","<<fHitStripList[k].stripPosition.Z()<<","<<endl;		
			}//close loop k
		}

		//## Tag this event
		eventTag thisEventTag= EventTagger(fPlaneList);
     
		fEventType.push_back(thisEventTag);
		cout<<"## Event Type:"<<thisEventTag<<endl;


		
		cout<<endl;
		cout<<endl;

}//close ReadScintillatorInfo()


void G4MuonCounterMCReader::ReadPMTInfo(int thisEvent){
	

	SimTree_pmt->GetEntry(thisEvent);
	
	if(fVerbosity>1) cout<< fPMTHitCollection->size() << " hit PMTs in this event" << endl;
		

	//## loop over all PMTs hit in this event
	for(unsigned int j=0;j<fPMTHitCollection->size();j++){
		int StripId= (fPMTHitCollection->at(j)).StripId;
		int PlaneId= (fPMTHitCollection->at(j)).PlaneId;
		int SuperPlaneId= (fPMTHitCollection->at(j)).SuperPlaneId;
		int PMTId= (fPMTHitCollection->at(j)).PMTId;
		int AbsPlaneId= 2*SuperPlaneId+PlaneId;
		int AbsStripId= AbsPlaneId*Nstrips+StripId;

		int PhotonCounts= (fPMTHitCollection->at(j)).PhotonCounts;
		int PhotonCounts_scint= (fPMTHitCollection->at(j)).PhotonCounts_scint;
		int PhotonCounts_cerenk= (fPMTHitCollection->at(j)).PhotonCounts_cerenk;
		int PhotonCounts_wls= (fPMTHitCollection->at(j)).PhotonCounts_wls;

		TVector3 PMTPos= (fPMTHitCollection->at(j)).Position;

		//## calculate horizontal distance hit-PMT
		//find position of this hit strip
		TVector3 StripHitPosition;
		double EnergyDeposit;
		double HitTime;
		for(unsigned int l=0;l<fHitStripList.size();l++){
			int thisAbsStripId= fHitStripList[l].absstripID;
			if(AbsStripId==thisAbsStripId){
				StripHitPosition= fHitStripList[l].stripHitPosition;
				EnergyDeposit= fHitStripList[l].etot;
				HitTime= fHitStripList[l].time;
				break;
			}
		}


		TVector3 HitPMTDistancePos= StripHitPosition-PMTPos;
		double HitPMTDistance=-1.;
		if(AbsPlaneId%2==0){//plane x
			//take distance in x
			HitPMTDistance= fabs(HitPMTDistancePos.X());
		}
		else if(AbsPlaneId%2==1){//plane y
			//take distance in y
			HitPMTDistance= fabs(HitPMTDistancePos.Y());
		}

		if(fVerbosity>1)
		{
			cout<<"HitPos="<<StripHitPosition.X()<<"  PMTPos="<<PMTPos.X()<<endl;
			cout<<"Distance hit-PMT="<<HitPMTDistance<<endl;
			cout<<"Energy Deposit="<<EnergyDeposit<<endl;		
		}
		//fPMTHitDistance[i].push_back(HitPMTDistance);
		

		//Get info for all photons hitting the PMT
		std::vector<double> photonTimeAtCathode= (fPMTHitCollection->at(j)).PhotonTime;
		std::vector<double> photonEnergyAtCathode= (fPMTHitCollection->at(j)).PhotonEnergy;
		std::vector<TVector3> photonPositionAtCathode= (fPMTHitCollection->at(j)).PhotonPosition;

		if(fVerbosity>1)
		{
			cout<<"## PMT Hit no "<< j+1 << endl;
			cout<<"StripId="<<StripId<<endl;
			cout<<"PlaneId="<<PlaneId<<endl;
			cout<<"SuperPlaneId="<<SuperPlaneId<<endl;
			cout<<"PMTId="<<PMTId<<endl;
			cout<<"PhotonCounts="<<PhotonCounts<<endl;
			cout<<"PhotonCounts_scint="<<PhotonCounts_scint<<endl;
			cout<<"PhotonCounts_cerenk="<<PhotonCounts_cerenk<<endl;
			cout<<"PhotonCounts_wls="<<PhotonCounts_wls<<endl;
			cout<<"Position=("<<PMTPos.X()<<","<<PMTPos.Y()<<","<<PMTPos.Z()<<")"<<endl;
		}		

		//## skip PMT below threshold in number of photon counts, otherwise do the PMT signal simulation
		//if(PhotonCounts < fPhotonThreshold) continue;	
		
		//## do PMT signal simulation
		PMTSimulator* fPMTSim= new PMTSimulator();
		fPMTSim->Init();
		
		if(fGeneratePMTPhotonsFast){
			//## generate photons from parametrization, given x & Edep
			fPMTSim->GeneratePhotonsFast(HitPMTDistance,EnergyDeposit);
		}
		else{
			//## use photons from G4 simulation
			fPMTSim->SetPhotonTimeAtCathode(photonTimeAtCathode);
			fPMTSim->SetPhotonEnergyAtCathode(photonEnergyAtCathode);	
			fPMTSim->SetPhotonPositionAtCathode(photonPositionAtCathode);	
		}

		if(fGeneratePMTPhotonsFast){
			//## generate photoelectrons from parametrization, given x & Edep
			fPMTSim->GeneratePEFast(HitPMTDistance,EnergyDeposit);
		}
		else{
			//## generate photoelectrons
			fPMTSim->GeneratePE();
		}

	
		//## generate pulses
		if(fGeneratePMTPulses) fPMTSim->PulseFinder();
		
		fPMTSimCollection.push_back(fPMTSim);

		int PECounts= fPMTSim->GetPECounts();

		//## skip PMT below threshold in number of photoelectrons
		//if(PECounts < fPEThreshold) continue;	

		//## Create struct for vector ordering	
		HitPMTStruct hitPMTStruct;
		hitPMTStruct.pmtID= PMTId;
		hitPMTStruct.stripID= StripId;
		hitPMTStruct.planeID= PlaneId;
		hitPMTStruct.superplaneID= SuperPlaneId;
		hitPMTStruct.absplaneID= AbsPlaneId;
		hitPMTStruct.absstripID= AbsStripId;
		hitPMTStruct.photonCounts= PhotonCounts;
		hitPMTStruct.peCounts= PECounts;
		hitPMTStruct.pmtPosition= PMTPos;
		hitPMTStruct.pmtHitPosition= HitPMTDistancePos;
		hitPMTStruct.pmtSim= fPMTSim;
		fHitPMTList.push_back(hitPMTStruct);
		
		//## Create and fill histograms
		bool isHistoAllocated= false;
		int HistoEntryId;
			
		for(unsigned int l=0;l<fSelectedPMTStripList.size();l++){//loop over all hits in current event
			if(AbsStripId==fSelectedPMTStripList[l] && PMTId==fSelectedPMTList[l]){
				isHistoAllocated=true; //this strip module has been already allocated as histograms
				HistoEntryId= l;
				break;
			}
		}//close loop

		if(!isHistoAllocated){
			//allocate and fill histograms
			TString histoName= Form("PMTCounts_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
			fPhotonCountsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
			fPhotonCountsHisto->Sumw2();
			fPhotonCountsHisto->Fill(HitPMTDistance, (double)(PhotonCounts));

			histoName= Form("PMTAverageCounts_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
			fAveragePhotonCountsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
			fAveragePhotonCountsHisto->Sumw2();
			fAveragePhotonCountsHisto->Fill(HitPMTDistance, (double)(PhotonCounts));

			histoName= Form("PECounts_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
			fPECountsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
			fPECountsHisto->Sumw2();
			fPECountsHisto->Fill(HitPMTDistance, (double)(PECounts));

			histoName= Form("PEAverageCounts_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
			fAveragePECountsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
			fAveragePECountsHisto->Sumw2();
			fAveragePECountsHisto->Fill(HitPMTDistance, (double)(PECounts));

			histoName= Form("NumberOfEvents_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
			fNumberOfEventsHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
			fNumberOfEventsHisto->Sumw2();
			fNumberOfEventsHisto->Fill(HitPMTDistance, 1.);

			histoName= Form("PhotonEfficiency_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
			fPhotonEfficiencyHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
			fPhotonEfficiencyHisto->Sumw2();
			if(PhotonCounts >= fPhotonThreshold) fPhotonEfficiencyHisto->Fill(HitPMTDistance, 1.);

			histoName= Form("PEEfficiency_superplane%d_plane%d_strip%d_pmt%d",SuperPlaneId,PlaneId,StripId,PMTId);
			fPEEfficiencyHisto= new TH1D(histoName,histoName,101,-0.5,100.5);
			fPEEfficiencyHisto->Sumw2();
			if(PECounts >= fPEThreshold) fPEEfficiencyHisto->Fill(HitPMTDistance, 1.);

			fSelectedPMTStripList.push_back(AbsStripId);
			fSelectedPMTList.push_back(PMTId);

			fPhotonCountsHistoList.push_back(fPhotonCountsHisto);
			fAveragePhotonCountsHistoList.push_back(fAveragePhotonCountsHisto);
			fPECountsHistoList.push_back(fPECountsHisto);
			fAveragePECountsHistoList.push_back(fAveragePECountsHisto);
			fNumberOfEventsHistoList.push_back(fNumberOfEventsHisto);
			fPhotonEfficiencyHistoList.push_back(fPhotonEfficiencyHisto);
			fPEEfficiencyHistoList.push_back(fPEEfficiencyHisto);
		}
		else{
			//retrieve and fill existing histograms
			fPhotonCountsHistoList[HistoEntryId]->Fill(HitPMTDistance, (double)(PhotonCounts));
			fAveragePhotonCountsHistoList[HistoEntryId]->Fill(HitPMTDistance, (double)(PhotonCounts));
			fPECountsHistoList[HistoEntryId]->Fill(HitPMTDistance, (double)(PECounts));
			fAveragePECountsHistoList[HistoEntryId]->Fill(HitPMTDistance, (double)(PECounts));
			fNumberOfEventsHistoList[HistoEntryId]->Fill(HitPMTDistance, 1.);
			if(PhotonCounts >= fPhotonThreshold) fPhotonEfficiencyHistoList[HistoEntryId]->Fill(HitPMTDistance, 1.);
			if(PECounts >= fPEThreshold) fPEEfficiencyHistoList[HistoEntryId]->Fill(HitPMTDistance, 1.);
		}

		fPhotonCountsMap->Fill(HitPMTDistance,EnergyDeposit,PhotonCounts);
		fPECountsMap->Fill(HitPMTDistance,EnergyDeposit,PECounts);
		x= HitPMTDistance;
		Edep= EnergyDeposit;
		phCounts= PhotonCounts;
		peCounts= PECounts;
		fCountsMapTree->Fill();

		t_hit= HitTime;
		for(unsigned int k=0;k<photonTimeAtCathode.size();k++){
			t= photonTimeAtCathode[k];
			TVector3 currentPhotonPos= photonPositionAtCathode[k];
			TVector3 currentDistance= currentPhotonPos - StripHitPosition;
			d= currentDistance.Mag();
			//cout<<"phPos=("<<currentPhotonPos.X()<<","<<currentPhotonPos.Y()<<","<<currentPhotonPos.Z()<<")   hitPos="<<StripHitPosition.X()<<","<<StripHitPosition.Y()<<","<<StripHitPosition.Z()<<")  d="<<d<<endl;
			fTimeMapTree->Fill();
		}

		delete fPMTSim;

	}//close loop PMT hits in this event


}//close ReadPMTInfo()



int G4MuonCounterMCReader::GetStripMultiplicity(int AbsPlaneId){  
  
	//multiplicity per detector plane
	int multiplicity= 0;

  for(unsigned int k=0; k<fHitStripList.size(); k++){ 
    int HitAbsPlaneId= fHitStripList[k].absplaneID;
		if(HitAbsPlaneId == AbsPlaneId) multiplicity++;
	}//loop on Selected Strip

 	cout<<"Multiplicity for plane= "<<AbsPlaneId+1<<" ==>"<<multiplicity<<endl;
         
	return multiplicity;  
    
}//close MCReader::GetStripMultiplicity()



std::vector<int> G4MuonCounterMCReader::CalculateGeomHits(){

	// calculate hit strips according to particle geometry
  // HitCoord is y-coordinate for a x-plane and viceversa for a x-plane
  double HitCoord_geom_up=0.;//hit in superior surface
  double HitCoord_geom_down=0.;//hit in bottom surface
  double HitCoord_geom_center=0.;//hit at center of the track

  int HitId_geom_up=-999;
  int HitId_geom_down=-999;
  int HitId_geom_center=-999;   

	std::vector<int> GeomHitStripId;
	GeomHitStripId.clear();
	GeomHitStripId.resize(2*Nplanes);
	for(unsigned int i=0;i<GeomHitStripId.size();i++) GeomHitStripId[i]= -1;

 
  //find geometrical intersections track-plane and store strip id 
  int n=0;
  for(unsigned int s=0;s<GeomHitStripId.size();s++){
  	
  	if(s%2==0){
    	//x plane (calculate y coord.)
      HitCoord_geom_up= (PlaneDepth[s-n]-Zcoat/2.)*tan(Theta*TMath::Pi()/180.)*sin(Phi*TMath::Pi()/180.) + Y0; 
      HitCoord_geom_down= (PlaneDepth[s]+Zcoat/2.)*tan(Theta*TMath::Pi()/180.)*sin(Phi*TMath::Pi()/180.) + Y0; 
                        
      //check if current track geom intersects more than one strip
      HitId_geom_up= (int)(ceil( (HitCoord_geom_up+Xcoat/2.)/Ycoat ) );        
      HitId_geom_down= (int)(ceil( (HitCoord_geom_down+Xcoat/2.)/Ycoat ) );        
      if(HitId_geom_up==HitId_geom_down){
      	GeomHitStripId[s]= HitId_geom_up;
      }
      else{
      	//find hottest strip according to greater track length        
      	HitCoord_geom_center= (PlaneDepth[s-n])*tan(Theta*TMath::Pi()/180.)*sin(Phi*TMath::Pi()/180.) + Y0; 
        HitId_geom_center= (int)(ceil( (HitCoord_geom_center+Xcoat/2.)/Ycoat ) );
        GeomHitStripId[s]= HitId_geom_center;        
      }
    }//close if plane x


    if(s%2==1){
    	//y plane (calculate x coord.)
      HitCoord_geom_up= (PlaneDepth[s-(n+1)]+Zcoat/2.+XYPlaneDistance)*tan(Theta*TMath::Pi()/180.)*cos(Phi*TMath::Pi()/180.) + X0; 
 	    HitCoord_geom_down= (PlaneDepth[s-(n+1)]+1.5*Zcoat+XYPlaneDistance)*tan(Theta*TMath::Pi()/180.)*cos(Phi*TMath::Pi()/180.) + X0; 
                        
     	//check if current track geom intersects more than one strip
      HitId_geom_up= (int)(ceil( (HitCoord_geom_up+Xcoat/2.)/Ycoat ) );        
      HitId_geom_down= (int)(ceil( (HitCoord_geom_down+Xcoat/2.)/Ycoat ) );        
      if(HitId_geom_up==HitId_geom_down){
      	GeomHitStripId[s]= HitId_geom_up;
      }
      else{
      	//find hottest strip according to greater track length        
       	HitCoord_geom_center= (PlaneDepth[s-(n+1)]+Zcoat+XYPlaneDistance)*tan(Theta*TMath::Pi()/180.)*cos(Phi*TMath::Pi()/180.) + X0; 
        HitId_geom_center= (int)(ceil( (HitCoord_geom_center+Xcoat/2.)/Ycoat ) );
        GeomHitStripId[s]= HitId_geom_center;
      }

      n++;
    }//close plane y

    if(fVerbosity>1)
		{  
    	cout<<"HitCoord_geom_up="<<HitCoord_geom_up<<" HitCoord_geom_down="<<HitCoord_geom_down;
    	cout<<"  HitCoord_geom_center="<<HitCoord_geom_center<<" GeomHitStripId="<<GeomHitStripId[s]<<endl;     
		}
	}//close plane s loop

 	return GeomHitStripId;

}//close MCReader::CalculateGeomHits()



void G4MuonCounterMCReader::DrawInfo(){
    
}

void G4MuonCounterMCReader::StoreInfo(){
    

	//##################################
	//#### STORE REC HISTO       #######
	//##################################
	//normalize average histograms, calculate efficiency histograms, write to file
	recDir->cd();
	for(unsigned int k=0;k<fAveragePhotonCountsHistoList.size();k++){
		fAveragePhotonCountsHistoList[k]->Divide(fAveragePhotonCountsHistoList[k],fNumberOfEventsHistoList[k],1,1);
		fAveragePECountsHistoList[k]->Divide(fAveragePECountsHistoList[k],fNumberOfEventsHistoList[k],1,1);
		fPhotonEfficiencyHistoList[k]->Divide(fPhotonEfficiencyHistoList[k],fNumberOfEventsHistoList[k],1,1,"B");
		fPEEfficiencyHistoList[k]->Divide(fPEEfficiencyHistoList[k],fNumberOfEventsHistoList[k],1,1,"B");

		//write to file	
		fPhotonCountsHistoList[k]->Write();
		fAveragePhotonCountsHistoList[k]->Write();
		fPECountsHistoList[k]->Write();
		fAveragePECountsHistoList[k]->Write();
		fNumberOfEventsHistoList[k]->Write();
		fPhotonEfficiencyHistoList[k]->Write();
		fPEEfficiencyHistoList[k]->Write();
	}
	fPhotonCountsMap->Write();
	fPECountsMap->Write();
	fCountsMapTree->Write();
	fTimeMapTree->Write();

	
	stripEdepDir->cd();
	for(unsigned int k=0;k<fEnergyDepositHistoList.size();k++){
		fEnergyDepositHistoList[k]->Write();
	}//close loop histo

	
	planeEdepDir->cd();
	for(unsigned int k=0;k<fEnergyDepositInPlaneHistoList.size();k++){
		fEnergyDepositInPlaneHistoList[k]->Write();
	}//close loop histo

	for(unsigned int k=0;k<fStripMultiplicityInPlaneHistoList.size();k++){
		fStripMultiplicityInPlaneHistoList[k]->Write();
	}//close loop histo


	//##################################
	//#### STORE GEN HISTO       #######
	//##################################
	genDir->cd();
	fGenEnergyHisto->Write();
	fGenMomentumXHisto->Write();
	fGenMomentumYHisto->Write();
	fGenMomentumZHisto->Write();
	fGenVertexXHisto->Write();
	fGenVertexYHisto->Write();
	fGenVertexZHisto->Write();
	fGenThetaHisto->Write();
	fGenPhiHisto->Write();

	fOutputFile->Close();

}//close MCReader::StoreInfo()


double G4MuonCounterMCReader::GetQEFromTable(double E){

	double QEProbability;
	std::map<double,double>::iterator it,itlow,itup;

	//check if given argument exists as key in the map
	it= QETable.find(E);
	if(it!=QETable.end()) {
		//key found
		QEProbability= QETable.find(E)->second;
		return QEProbability;
	}
	
	//if not, do a linear interpolation
  itlow=QETable.lower_bound(E);  
  itup=QETable.upper_bound(E);  
	
	//check ranges: if low-bound=end() or up_bound=begin() return 0
	if(itlow==QETable.end()) {
		//Key outside upper-range...return 0
		return 0.;
	}
	else if(itup==QETable.begin()) {
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

}//close MCReader::GetQEFromTable()


void G4MuonCounterMCReader::CreateTrackPoints(){

	//## Create collection of track points
	fTrackPointCollection.clear();
	fTrackPointCollection.resize(0);
	

	for(unsigned int k=0;k<fHitStripList.size();k++){//loop over all hits in current event		
		int AbsPlaneId= fHitStripList[k].absplaneID;
		int SuperPlaneId= fHitStripList[k].superplaneID;
			
		if(AbsPlaneId%2==0){//plane X
			//search corresponding position in Y plane						
			for(unsigned int j=k;j<fHitStripList.size();j++){//loop over all following hits in current event
				int AbsPlaneId_next= fHitStripList[j].absplaneID;
				if(AbsPlaneId_next%2==1 && (AbsPlaneId_next-AbsPlaneId)==1) {//select consecutive X-Y planes
					TVector3 XPosition= fHitStripList[k].stripPosition;
					TVector3 YPosition= fHitStripList[j].stripPosition;
					TVector3 PointPosition= TVector3(YPosition.X(), XPosition.Y(), 0.5*(XPosition.Z()+YPosition.Z()) );

					TrackPoint aTrackPoint;
					aTrackPoint.fPosition= PointPosition;	
					aTrackPoint.fDetectorPlaneId= SuperPlaneId;
	
					fTrackPointCollection.push_back(aTrackPoint);					
				}
				else continue;//skip plane X entry or not consecutive X-Y planes	
			}	
		}//close if plane X
		else continue;//skip plane Y entry

	}//close loop k


}//CreateTrackPoints()



G4MuonCounterMCReader::eventTag G4MuonCounterMCReader::EventTagger(std::vector<int> hitPlaneList){

	G4MuonCounterMCReader::eventTag thisEventTag;		
	
	int q=0;
  int index_x;
  int index_y;
	int planeMultiplicity=0;
		
	while(q< hitPlaneList.size()){
		//check for X-Y combinations			
		if( (hitPlaneList[q]%2==0) && (hitPlaneList[q+1]%2==1) && (hitPlaneList[q+1]-hitPlaneList[q]==1) ){
			planeMultiplicity++;
      index_x=q;
      index_y=q+1;
      q+=2;
		}
		else q+=1;
	}//close while	

	//Determine event type
	switch(planeMultiplicity){
		case kONEFOLD:
		{
			thisEventTag= kONEFOLD;
			break;
		}
		case kTWOFOLD:
		{
			thisEventTag= kTWOFOLD;
			break;
		}
		case kTHREEFOLD:
		{
			thisEventTag= kTHREEFOLD;
			break;
		}
		default:
			thisEventTag= kEMPTY;
			break;	
	}//close switch

	return thisEventTag;

}//close MCReader::EventTagger()

