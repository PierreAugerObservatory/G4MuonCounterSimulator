/**
* @file EventFileRecorder.cc
* @class EventFileRecorder
* @brief Record the event gen/sim/rec information in a ROOT file
*
* @author S. Riggi
* @date 18/01/2011
*/

static const char CVSId[] = 
"$Id$";

#include <EventFileRecorder.h>
#include <G4MuonCounterSimulator.hh>
#include <G4MuonCounterReconstructor.hh>
#include <MuonDetector.hh>
#include <TEventSimData.hh>

#include <fwk/CentralConfig.h>
#include <det/Detector.h>
#include <io/EventFile.h>
#include <evt/Event.h>
#include <evt/ShowerSimData.h>

#include <sevt/SEvent.h>
#include <sevt/Header.h>
#include <sevt/EventTrigger.h>

#include <fevt/FEvent.h>
#include <fevt/Eye.h>
#include <fevt/EyeHeader.h>

#include <utl/Particle.h>
#include <utl/ErrorLogger.h>
#include <utl/Reader.h>
#include <utl/TimeStamp.h>

#include <EyeEvent.hh>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;


#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TSystem.h>
#include <Rtypes.h>
#include <TObject.h>
#include <TRint.h>
#include <TMath.h>

#include <vector>
#include <cstddef>
#include <iostream>
#include <cstddef>
#include <sstream>

using namespace G4MuonCounterSimulatorUSC;
using namespace G4MuonCounterReconstructorUSC;
using namespace fwk;
using namespace io;
using namespace evt;
using namespace utl;
using namespace std;

EventFileRecorder::EventFileRecorder() : 
  fOutputFile(0),
	fGenInfo(0),
	fDetInfo(0),
	fSimInfo(0),
	fRecInfo(0),
	fSaveGenInfo(0),
	fSaveDetInfo(0),
	fSaveSimInfo(0),
	fSaveRecInfo(0), 
  fEventsWritten(0) {
}

EventFileRecorder::~EventFileRecorder() {
  //delete fOutputFile;
}

VModule::ResultFlag EventFileRecorder::Init() {

  CentralConfig *theConfig = CentralConfig::GetInstance();
  Branch topB = theConfig->GetTopBranch("EventFileRecorder");
  
	Branch saveSimInfoB = topB.GetChild("SaveInfo");
  if (saveSimInfoB) {
    saveSimInfoB.GetChild("OutFileName").GetData(fOutputFileName);
		saveSimInfoB.GetChild("GenInfo").GetData(fSaveGenInfo);	
		saveSimInfoB.GetChild("DetInfo").GetData(fSaveDetInfo);	
		saveSimInfoB.GetChild("SimInfo").GetData(fSaveSimInfo);	
		saveSimInfoB.GetChild("RecInfo").GetData(fSaveRecInfo);		
  }


  fOutputFile = new TFile(fOutputFileName.c_str(), "RECREATE");
	if(fSaveGenInfo){
		fMuonLateralDensity= new TH1D("MuonLateralDensity","MuonLateralDensity",40,0,2000);
		fMuonLateralDensity->Sumw2();
		fMuonLateralDensity->GetXaxis()->SetTitle("r [m]");
		fMuonLateralDensity->GetXaxis()->SetTitleSize(0.06);
		fMuonLateralDensity->GetXaxis()->SetTitleOffset(0.7);
		fMuonLateralDensity->GetYaxis()->SetTitle("#rho_{#mu} [m^{-2}]");
		fMuonLateralDensity->GetYaxis()->SetTitleSize(0.06);
		fMuonLateralDensity->GetYaxis()->SetTitleOffset(0.7);

		fMuonLateralDensityGrd= new TH1D("MuonLateralDensityGrd","MuonLateralDensityGrd",40,0,2000);
		fMuonLateralDensityGrd->Sumw2();
		fMuonLateralDensityGrd->GetXaxis()->SetTitle("r [m]");
		fMuonLateralDensityGrd->GetXaxis()->SetTitleSize(0.06);
		fMuonLateralDensityGrd->GetXaxis()->SetTitleOffset(0.7);
		fMuonLateralDensityGrd->GetYaxis()->SetTitle("#rho_{#mu} [m^{-2}]");
		fMuonLateralDensityGrd->GetYaxis()->SetTitleSize(0.06);
		fMuonLateralDensityGrd->GetYaxis()->SetTitleOffset(0.7);

		fMuonLateralNumber= new TH1D("MuonLateralNumber","MuonLateralNumber",40,0,2000);
		fMuonLateralNumber->Sumw2();
		fMuonLateralNumber->GetXaxis()->SetTitle("r [m]");
		fMuonLateralNumber->GetXaxis()->SetTitleSize(0.06);
		fMuonLateralNumber->GetXaxis()->SetTitleOffset(0.7);
		fMuonLateralNumber->GetYaxis()->SetTitle("N_{#mu}");
		fMuonLateralNumber->GetYaxis()->SetTitleSize(0.06);
		fMuonLateralNumber->GetYaxis()->SetTitleOffset(0.7);

		fMuonLateralNumberGrd= new TH1D("MuonLateralNumberGrd","MuonLateralNumberGrd",40,0,2000);
		fMuonLateralNumberGrd->Sumw2();
		fMuonLateralNumberGrd->GetXaxis()->SetTitle("r [m]");
		fMuonLateralNumberGrd->GetXaxis()->SetTitleSize(0.06);
		fMuonLateralNumberGrd->GetXaxis()->SetTitleOffset(0.7);
		fMuonLateralNumberGrd->GetYaxis()->SetTitle("N_{#mu}");
		fMuonLateralNumberGrd->GetYaxis()->SetTitleSize(0.06);
		fMuonLateralNumberGrd->GetYaxis()->SetTitleOffset(0.7);

		fGenInfo= new TTree("GenInfo","GenInfo");	
		fGenInfo->Branch("lgE",&lgE,"lgE/D");
	  fGenInfo->Branch("Theta",&Theta,"Theta/D");
		fGenInfo->Branch("Phi",&Phi,"Phi/D");
		fGenInfo->Branch("Rp",&Rp,"Rp/D");
		fGenInfo->Branch("Psi",&Psi,"Psi/D");
		fGenInfo->Branch("GPSTimeSec",&GPSTimeSec,"GPSTimeSec/D");	
		fGenInfo->Branch("GPSTimeNanoSec",&GPSTimeNanoSec,"GPSTimeNanoSec/D");
		fGenInfo->Branch("CoreX",&CoreX,"CoreX/D");
		fGenInfo->Branch("CoreY",&CoreY,"CoreY/D");
		fGenInfo->Branch("nMu",&nMu,"nMu/D");
		fGenInfo->Branch("nMuPlus",&nMuPlus,"nMuPlus/D");
		fGenInfo->Branch("nMuMinus",&nMuMinus,"nMuMinus/D");
		fGenInfo->Branch("nE",&nE,"nE/D");
		fGenInfo->Branch("nEPlus",&nEPlus,"nEPlus/D");
		fGenInfo->Branch("nEMinus",&nEMinus,"nEMinus/D");
		fGenInfo->Branch("nGamma",&nGamma,"nGamma/D");
		fGenInfo->Branch("nProton",&nProton,"nProton/D");
		fGenInfo->Branch("nNeutron",&nNeutron,"nNeutron/D");
		fGenInfo->Branch("HadronicModel",&HadronicModel,"HadronicModel/I");
		fGenInfo->Branch("PrimaryParticle",&PrimaryParticle,"PrimaryParticle/I");
		fGenInfo->Branch("ShowerNumber",&ShowerNumber,"ShowerNumber/I");
		fGenInfo->Branch("ShowerRunId",&ShowerRunId,"ShowerRunId/C");

		fGenInfo->Branch("MuonLateralDensity","TH1D",&fMuonLateralDensity);
		fGenInfo->Branch("MuonLateralDensityGrd","TH1D",&fMuonLateralDensityGrd);
		fGenInfo->Branch("MuonLateralNumber","TH1D",&fMuonLateralNumber);
		fGenInfo->Branch("MuonLateralNumberGrd","TH1D",&fMuonLateralNumberGrd);
	}

	if(fSaveDetInfo){
		fDetInfo= new TTree("DetInfo","DetInfo");	
		fDetInfo->Branch("Detector","MuonDetector",&fMuonDetector);
	}

	if(fSaveSimInfo){
		fSimInfo= new TTree("SimInfo","SimInfo");	
		fSimInfo->Branch("EventSimData","TEventSimData",&fEventSimData);
	}

	if(fSaveRecInfo){
		fRecInfo= new TTree("RecInfo","RecInfo");	
		//fRecInfo->Branch("EventRecData","TEventRecData",&fCurrentEventSimData);
	}
	

  return eSuccess;
}

VModule::ResultFlag EventFileRecorder::Run(evt::Event& theEvent) {

  if(fSaveGenInfo) RecordGenInfo(theEvent);
	if(fSaveDetInfo && fEventsWritten==0) RecordDetInfo(); //save detector info just once
	if(fSaveSimInfo) RecordSimInfo();
	if(fSaveRecInfo) RecordRecInfo();

  fEventsWritten++;
  
  return eSuccess;

}

VModule::ResultFlag EventFileRecorder::RecordGenInfo(evt::Event& event) {

	if(!event.HasSimShower()) return eContinueLoop;  
  if(!event.HasSEvent()) return eContinueLoop;
  
	fMuonLateralDensity->Reset();
	fMuonLateralDensityGrd->Reset();
	fMuonLateralNumber->Reset();
	fMuonLateralNumberGrd->Reset();

  const CoordinateSystemPtr siteCS = det::Detector::GetInstance().GetSiteCoordinateSystem();
  const ShowerSimData & simShower = event.GetSimShower();  
  const CoordinateSystemPtr showerCS = simShower.GetShowerCoordinateSystem();
	const CoordinateSystemPtr localCS = simShower.GetLocalCoordinateSystem();
  

	//## get shower vars
  lgE = log10(simShower.GetEnergy()/utl::eV); 
  Theta = simShower.GetZenith()/utl::degree;
  Phi   = simShower.GetAzimuth()/utl::degree; 
  CoreX  = simShower.GetPosition().GetX(siteCS)/utl::m; 
  CoreY  = simShower.GetPosition().GetY(siteCS)/utl::m; 
	Rp  = simShower.GetDirection().GetRho(showerCS)/utl::m; 
	Psi = simShower.GetDirection().GetPhi(showerCS)/utl::degree; 
  rCut = simShower.GetMaxRadiusCut()/utl::m;
	GPSTimeSec     = simShower.GetTimeStamp().GetGPSSecond();
	GPSTimeNanoSec = simShower.GetTimeStamp().GetGPSNanoSecond();
	PrimaryParticle = simShower.GetPrimaryParticle();
	ShowerNumber = simShower.GetShowerNumber();
	ShowerRunId = simShower.GetShowerRunId().c_str();

	//## loop over ground particles
  nMu       = 0;
	nMuPlus   = 0;
	nMuMinus  = 0;
  nE        = 0;
	nEPlus    = 0;
	nEMinus   = 0;
	nGamma    = 0;
	nProton   = 0;
	nNeutron  = 0;
	
	int particleCounter= 0;

	for(ShowerParticleIterator particleIt = simShower.GroundParticlesBegin(); particleIt != simShower.GroundParticlesEnd(); ++particleIt) {
    
		Point particlePosition = particleIt->GetPosition();
    int particleType = particleIt->GetType();
		double weight= particleIt->GetWeight();
		double xGrd = particlePosition.GetX(localCS);
		double yGrd = particlePosition.GetY(localCS);  
		double zGrd = particlePosition.GetZ(localCS);
		double rGrd = particlePosition.GetRho(localCS);
		double psiGrd= particlePosition.GetPhi(localCS);
		
		//ground coordinates in shower plane system  
		double x = particlePosition.GetX(showerCS); 
		double y = particlePosition.GetY(showerCS);  
		double z = particlePosition.GetZ(showerCS);
   	double r = particlePosition.GetRho(showerCS);
    double psi = particlePosition.GetPhi(showerCS);

		//find area
		double rbinGrd= fMuonLateralDensityGrd->GetXaxis()->FindBin(rGrd);
		double rDeltaGrd= fMuonLateralDensityGrd->GetXaxis()->GetBinWidth(rbinGrd);
		double r1Grd= fMuonLateralDensityGrd->GetXaxis()->GetBinLowEdge(rbinGrd);
		double r2Grd= r1Grd + rDeltaGrd;		
		double rAreaGrd= TMath::Pi()* (r2Grd*r2Grd-r1Grd*r1Grd);

		double rbin= fMuonLateralDensity->GetXaxis()->FindBin(r);
		double rDelta= fMuonLateralDensity->GetXaxis()->GetBinWidth(rbin);
		double r1= fMuonLateralDensity->GetXaxis()->GetBinLowEdge(rbin);
		double r2= r1 + rDelta;		
		double rArea= TMath::Pi()* (r2*r2-r1*r1);
	

		switch(particleType){
      	case utl::Particle::eMuon:
				{
					nMuMinus+= weight;
					fMuonLateralDensityGrd->Fill(rGrd,weight/rAreaGrd);
					fMuonLateralDensity->Fill(r,weight/rArea);
					fMuonLateralNumberGrd->Fill(rGrd,weight);
					fMuonLateralNumber->Fill(r,weight);

					break;
				}
				case utl::Particle::eAntiMuon:
				{
					nMuPlus+= weight;
					fMuonLateralDensityGrd->Fill(rGrd,weight/rAreaGrd);
					fMuonLateralDensity->Fill(r,weight/rArea);
					fMuonLateralNumberGrd->Fill(rGrd,weight);
					fMuonLateralNumber->Fill(r,weight);

					break;
				}
				case utl::Particle::eElectron:
				{
					nEMinus+= weight;
					break;
				}
				case utl::Particle::ePositron:
				{
					nEPlus+= weight;
					break;
				}
				case utl::Particle::ePhoton:
				{
					nGamma+= weight;
					break;
				}
				case utl::Particle::eProton:
				{
					nProton+= weight;
					break;
				}
				case utl::Particle::eNeutron:
				{
					nNeutron+= weight;
					break;
				}
				default:
				{
					continue;
				}
				break;
		}//end switch

		particleCounter++;

	}//end loop particles

	nMu= nMuMinus + nMuPlus;
	nE= nEMinus + nEPlus;
	
	

	//## filling tree
	fGenInfo->Fill();

	return eSuccess;
}

VModule::ResultFlag EventFileRecorder::RecordDetInfo() {

	//get pointer to current TEventSimData from G4MuonCounterSimulator module
	fMuonDetector= G4MuonCounterSimulator::GetCurrentMuonDetector();
	
	//if simulator module not present, cannot get pointer, throw error
	if(!fMuonDetector){
		ostringstream errinfo;
  	errinfo << "Cannot get current MuonDetector from G4MuonCounterSimulator module";         
  	ERROR(errinfo.str().c_str());
		return eFailure;
	}

	//fill tree with current TEventSimData
	fDetInfo->Fill();

	return eSuccess;
}

VModule::ResultFlag EventFileRecorder::RecordSimInfo() {

	//get pointer to current TEventSimData from G4MuonCounterSimulator module
	fEventSimData= G4MuonCounterSimulator::GetCurrentEventSimData();
	
	//if simulator module not present, cannot get pointer, throw error
	if(!fEventSimData){
		ostringstream errinfo;
  	errinfo << "Cannot get current TEventSimData from G4MuonCounterSimulator module";         
  	ERROR(errinfo.str().c_str());
		return eFailure;
	}

	//fill tree with current TEventSimData
	fSimInfo->Fill();

	return eSuccess;
}

VModule::ResultFlag EventFileRecorder::RecordRecInfo() {

	return eSuccess;
}


VModule::ResultFlag EventFileRecorder::Finish() {
  
  std::ostringstream info;
  
  info << "Wrote " << fEventsWritten << " events.";  
  INFO(info); 
  
	fOutputFile->cd();
	if(fSaveGenInfo) fGenInfo->Write();
	if(fSaveDetInfo) fDetInfo->Write();
	if(fSaveSimInfo) fSimInfo->Write();
	if(fSaveRecInfo) fRecInfo->Write();

  fOutputFile->Close();
  
  return eSuccess;

}

