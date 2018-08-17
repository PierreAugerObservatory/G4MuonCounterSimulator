
/**
* @file G4MuonCounterReconstructor.cc
* @class module
* @brief module used to handle the reconstruction of simulated data
*
* @author S. Riggi
* @date 05/04/2010
*/

#include <G4MuonCounterReconstructor.hh>

#include <RVersion.h>

#include <TEventSimData.hh>
#include <TParticleSimData.hh>
#include <TStationSimData.hh>
#include <TScintHit.hh>
#include <TPMTHit.hh>
#include <TrackPoint.hh>
#include <TrackReco.hh>
#include <AnalysisConsts.hh>
#include <G4UnitsTable.hh>
#include <PMTSimulator.hh>
#include <TrackFinder.hh>
#include <KFTrackFinder.hh>

#include <Utilities.hh>


#include <fwk/CentralConfig.h>
#include <fwk/RandomEngineRegistry.h>
#include <fwk/RunController.h>

#include <det/Detector.h>

#include <evt/Event.h>
#include <evt/ShowerSimData.h>

#include <sdet/SDetector.h>
#include <sdet/Station.h>

#include <sevt/PMT.h>
#include <sevt/PMTSimData.h>
#include <sevt/SEvent.h>
#include <sevt/Station.h>
#include <sevt/StationSimData.h>

#include <utl/ErrorLogger.h>
#include <utl/Reader.h>
#include <utl/Particle.h>
#include <utl/ShowerParticleIterator.h>
#include <utl/TimeDistribution.h>
#include <utl/TimeDistributionAlgorithm.h>

#include <CLHEP/Random/Random.h>


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
#include <TColor.h>


#include <TRandom.h>
#include <Rtypes.h>
#include <TObject.h>
#include <TRint.h>

#include <vector>
#include <cstddef>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>


using namespace utl;
using namespace fwk;
using namespace std;
using namespace sevt;
using namespace sdet;
using namespace evt;
using namespace det;

using namespace G4MuonCounterSimulatorUSC;
using namespace G4MuonCounterReconstructorUSC;

G4MuonCounterReconstructor::G4MuonCounterReconstructor() :
	fSaveRecInfo(true),
  fReadSimInfoFromFile(false),
  fEnergyThreshold(0),
	fPEThreshold(0),
  fPhotonThreshold(0),
	fNevents(0),
	fCurrentEventNo(0),
	fCurrentStation(0),
	fVerbosity(1),
	fProcessSelectedHits(0),

	fTrackMuons(0),
	fTrackingAlgo(2),
	fTrackWithMSCorrection(1),
	fUseAverageEnergyInMS(1),
	fAverageEnergyInMS(5000),
	fUseTrueEnergyInMS(0),
	fUseEnergyLossCorrectionInMS(0),
	fIncludeMomentumInTracking(0),
	fUseClusterInTracking(1),
	fSplitClustersInTracking(0),	
	fSplitClusterThreshold(2),
	fRemoveFakeHitsInTracking(0),
	fUsePID(0)

{
	
}


G4MuonCounterReconstructor::~G4MuonCounterReconstructor()
{
	//if(currentSimParticle) delete currentSimParticle;
	
}


VModule::ResultFlag
G4MuonCounterReconstructor::Init()
{

	INFO("Geant4 muon counter rec (beta-version)");
	//##################################
	//## get some XML config
	//##################################
  Branch topB = CentralConfig::GetInstance()->GetTopBranch("G4MuonCounterReconstructorUSC");
	if (!topB) {
    ERROR("No XML configuration found for G4MuonCounterReconstructorUSC!");
    return eFailure;
  }
	
	map<string, string> useAtt;
  useAtt["use"] = string("yes");

	Branch readSimInfoB = topB.GetChild("readSimInfoFromFile");
	if(!readSimInfoB){
		ERROR("No readSimInfoFromFile branch found for G4MuonCounterReconstructorUSC!");
    return eFailure;
	} 
  readSimInfoB.GetData(fReadSimInfoFromFile);	
	
	Branch inputFileB = topB.GetChild("InputFileName");
	if(! inputFileB){
		ERROR("No InputFileName branch found for G4MuonCounterReconstructorUSC!");
    return eFailure;
	}
	inputFileB.GetData(fInputFileName);	
  
	Branch processHitsB = topB.GetChild("processHits", useAtt);
	if(processHitsB){
		fProcessHits= true;
		processHitsB.GetChild("processSelectedHits").GetData(fProcessSelectedHits);	
		processHitsB.GetChild("energyDepositThreshold").GetData(fEnergyThreshold);	
  	processHitsB.GetChild("nPhotonThreshold").GetData(fPhotonThreshold);	
  	processHitsB.GetChild("nPEThreshold").GetData(fPEThreshold);	
	}	
	else fProcessHits= false;

		
	Branch mergeParticlesInEventB = topB.GetChild("mergeParticlesInEvent", useAtt);
	if(mergeParticlesInEventB){
		fMergeSimParticles= true;
		mergeParticlesInEventB.GetChild("mergeTimeInterval").GetData(fMergeTimeInterval);	
		mergeParticlesInEventB.GetChild("samplingTime").GetData(fSamplingTime);	
	}	
	else fMergeSimParticles= false;
	

	Branch trackMuonsB = topB.GetChild("trackMuons", useAtt);
	if(trackMuonsB){
		fTrackMuons= true;
		trackMuonsB.GetChild("TrackingAlgo").GetData(fTrackingAlgo);
		trackMuonsB.GetChild("UseClusters").GetData(fUseClusterInTracking);
		trackMuonsB.GetChild("SplitClusters").GetData(fSplitClustersInTracking);
		trackMuonsB.GetChild("SplitClusterThreshold").GetData(fSplitClusterThreshold);
		trackMuonsB.GetChild("RemoveFakeHits").GetData(fRemoveFakeHitsInTracking);
		trackMuonsB.GetChild("UsePID").GetData(fUsePID);	
		trackMuonsB.GetChild("NNCutInPID").GetData(fNNCutInPID);
		
		Branch MSModeB= trackMuonsB.GetChild("UseMSCorrection", useAtt);
		if(MSModeB){
			fTrackWithMSCorrection= true;
			MSModeB.GetChild("UseAverageEnergy").GetData(fUseAverageEnergyInMS);
			MSModeB.GetChild("AverageEnergy").GetData(fAverageEnergyInMS);
			MSModeB.GetChild("UseTrueEnergy").GetData(fUseTrueEnergyInMS);
			MSModeB.GetChild("UseEnergyLossCorrection").GetData(fUseEnergyLossCorrectionInMS);
		}
		else fTrackWithMSCorrection= false;
  	
  	trackMuonsB.GetChild("IncludeMomentum").GetData(fIncludeMomentumInTracking);	
	}		
	else fTrackMuons= false;

	cout<<"*** TRACKING STEERINGS ***"<<endl;
	cout<<"fTrackMuons? "<< fTrackMuons<<endl;
	cout<<"fTrackingAlgo="<<fTrackingAlgo<<endl;
	cout<<"fUseClusterInTracking? "<<fUseClusterInTracking<<endl;
	cout<<"fSplitClustersInTracking? "<<fSplitClustersInTracking<<endl;
	cout<<"fRemoveFakeHitsInTracking? "<<fRemoveFakeHitsInTracking<<endl;
	cout<<"fTrackWithMSCorrection? "<<fTrackWithMSCorrection<<endl;
	cout<<"fUsePID? "<<fUsePID<<endl;
	cout<<"fNNCutInPID? "<<fNNCutInPID<<endl;
	
	if(fTrackWithMSCorrection){
		cout<<"fUseAverageEnergyInMS? "<<fUseAverageEnergyInMS<<endl;
		cout<<"fAverageEnergyInMS [MeV]= "<<fAverageEnergyInMS/utl::MeV<<endl;
		cout<<"fUseTrueEnergyInMS ?"<<fUseTrueEnergyInMS<<endl;
		cout<<"fUseEnergyLossCorrectionInMS? "<<fUseEnergyLossCorrectionInMS<<endl;
	}
	cout<<"IncludeMomentum? "<<fIncludeMomentumInTracking<<endl;
	
	
		
	Branch saveRecInfoB = topB.GetChild("SaveRecInfo", useAtt);
  if (saveRecInfoB) {
    fSaveRecInfo = true;
		saveRecInfoB.GetChild("saveRecTree").GetData(fSaveRecTree);	
		saveRecInfoB.GetChild("saveTrackTree").GetData(fSaveTrackTree);	
		saveRecInfoB.GetChild("saveHisto").GetData(fSaveHisto);	
    saveRecInfoB.GetChild("OutFileName").GetData(fOutputFileName);
  }


	cout<<"Saving rec output in "<< fOutputFileName.c_str()<< endl;
	fOutputFile= new TFile(fOutputFileName.c_str(),"RECREATE");

	RecTree= new TTree("RecTree","RecTree");
	RecTree->Branch("Id",&Id,"Id/I");
	RecTree->Branch("E",&E,"E/D");
	RecTree->Branch("Theta",&Theta,"Theta/D");
	RecTree->Branch("Phi",&Phi,"Phi/D");	
	RecTree->Branch("ThetaRec",&ThetaRec,"ThetaRec/D");
	RecTree->Branch("PhiRec",&PhiRec,"PhiRec/D");		
	RecTree->Branch("Time",&Time,"Time/D");
	RecTree->Branch("X",&X,"X/D");
	RecTree->Branch("Y",&Y,"Y/D");
	RecTree->Branch("StationId",&StationId,"StationId/I");
	RecTree->Branch("RStat",&RStat,"RStat/D");
	RecTree->Branch("RStatGrd",&RStatGrd,"RStatGrd/D");
	RecTree->Branch("R",&R,"R/D");
	RecTree->Branch("RGrd",&RGrd,"RGrd/D");
	RecTree->Branch("Psi",&Psi,"Psi/D");
	RecTree->Branch("PsiGrd",&PsiGrd,"PsiGrd/D");
	RecTree->Branch("EventTag",&EventTag,"EventTag/I");
	RecTree->Branch("GeomEventTag",&GeomEventTag,"GeomEventTag/I");
	RecTree->Branch("HasGeomCrossedFirstPlane",&HasGeomCrossedFirstPlane,"HasGeomCrossedFirstPlane/I");
	RecTree->Branch("nHitPlanes",&nHitPlanes,"nHitPlanes/I");
	RecTree->Branch("AverageXHit",AverageXHit,"AverageXHit[nHitPlanes]/D");
	RecTree->Branch("AverageYHit",AverageYHit,"AverageYHit[nHitPlanes]/D");
	RecTree->Branch("MultiplicityX",MultiplicityX,"MultiplicityX[nHitPlanes]/D");
	RecTree->Branch("MultiplicityY",MultiplicityY,"MultiplicityY[nHitPlanes]/D");
	RecTree->Branch("RX",RX,"RX[nHitPlanes]/D");
	RecTree->Branch("RY",RY,"RY[nHitPlanes]/D");
	RecTree->Branch("nMuonsX",nMuonsX,"nMuonsX[nHitPlanes]/D");
	RecTree->Branch("nMuonsAlbedoX",nMuonsAlbedoX,"nMuonsAlbedoX[nHitPlanes]/D");
	RecTree->Branch("nMuonsY",nMuonsY,"nMuonsY[nHitPlanes]/D");
	RecTree->Branch("nMuonsAlbedoY",nMuonsAlbedoY,"nMuonsAlbedoY[nHitPlanes]/D");
	RecTree->Branch("nEmX",nEmX,"nEmX[nHitPlanes]/D");
	RecTree->Branch("nEmAlbedoX",nEmAlbedoX,"nEmAlbedoX[nHitPlanes]/D");
	RecTree->Branch("nEmY",nEmY,"nEmY[nHitPlanes]/D");
	RecTree->Branch("nEmAlbedoY",nEmAlbedoY,"nEmAlbedoY[nHitPlanes]/D");
	RecTree->Branch("nHadronsX",nHadronsX,"nHadronsX[nHitPlanes]/D");
	RecTree->Branch("nHadronsAlbedoX",nHadronsAlbedoX,"nHadronsAlbedoX[nHitPlanes]/D");
	RecTree->Branch("nHadronsY",nHadronsY,"nHadronsY[nHitPlanes]/D");
	RecTree->Branch("nHadronsAlbedoY",nHadronsAlbedoY,"nHadronsAlbedoY[nHitPlanes]/D");

	TrackTree= new TTree("TrackTree","TrackTree");
	TrackTree->Branch("StationId",&StationId,"StationId/I");
	TrackTree->Branch("StationR",&RStat,"StationR/D");
	TrackTree->Branch("nTrueTrack",&nTrueTrack,"nTrueTrack/I");
	//TrackTree->Branch("nTrueMuonTrack",&nTrueMuonTrack,"nTrueMuonTrack/I");
	TrackTree->Branch("GenE",GenE,"GenE[nTrueTrack]/D");
	TrackTree->Branch("GenTheta",GenTheta,"GenTheta[nTrueTrack]/D");
	TrackTree->Branch("GenPhi",GenPhi,"GenPhi[nTrueTrack]/D");
	TrackTree->Branch("GenX",GenX,"GenX[nTrueTrack]/D");
	TrackTree->Branch("GenY",GenY,"GenY[nTrueTrack]/D");
	TrackTree->Branch("GenTx",GenTx,"GenTx[nTrueTrack]/D");
	TrackTree->Branch("GenTy",GenTy,"GenTy[nTrueTrack]/D");
	TrackTree->Branch("GenR",GenR,"GenR[nTrueTrack]/D");	
	TrackTree->Branch("GenMaxRelAngle",GenMaxRelAngle,"GenMaxRelAngle[nTrueTrack]/D");	
	
	TrackTree->Branch("nTrack",&nTrack,"nTrack/I");
	TrackTree->Branch("RecOpeningAngle",RecOpeningAngle,"RecOpeningAngle[nTrack]/D");
	TrackTree->Branch("RecTheta",RecTheta,"RecTheta[nTrack]/D");
	TrackTree->Branch("RecPhi",RecPhi,"RecPhi[nTrack]/D");
	TrackTree->Branch("RecX",RecX,"RecX[nTrack]/D");
	TrackTree->Branch("RecY",RecY,"RecY[nTrack]/D");
	TrackTree->Branch("RecTx",RecTx,"RecTx[nTrack]/D");
	TrackTree->Branch("RecTy",RecTy,"RecTy[nTrack]/D");
	TrackTree->Branch("RecThetaErr",RecThetaErr,"RecThetaErr[nTrack]/D");
	TrackTree->Branch("RecPhiErr",RecPhiErr,"RecPhiErr[nTrack]/D");
	TrackTree->Branch("RecXErr",RecXErr,"RecXErr[nTrack]/D");
	TrackTree->Branch("RecYErr",RecYErr,"RecYErr[nTrack]/D");
	TrackTree->Branch("RecTxErr",RecTxErr,"RecTxErr[nTrack]/D");
	TrackTree->Branch("RecTyErr",RecTyErr,"RecTyErr[nTrack]/D");
	TrackTree->Branch("TrackRecoStatus",&TrackRecoStatus,"TrackRecoStatus/I");
	TrackTree->Branch("TrackFitStatus",TrackFitStatus,"TrackFitStatus[nTrack]/I");
	TrackTree->Branch("TrackFitChi2",TrackFitChi2,"TrackFitChi2[nTrack]/D");
	TrackTree->Branch("RecThetaStart",RecThetaStart,"RecThetaStart[nTrack]/D");
	TrackTree->Branch("RecPhiStart",RecPhiStart,"RecPhiStart[nTrack]/D");
	TrackTree->Branch("RecXStart",RecXStart,"RecXStart[nTrack]/D");
	TrackTree->Branch("RecYStart",RecYStart,"RecYStart[nTrack]/D");
	TrackTree->Branch("RecTxStart",RecTxStart,"RecTxStart[nTrack]/D");
	TrackTree->Branch("RecTyStart",RecTyStart,"RecTyStart[nTrack]/D");
	TrackTree->Branch("RecThetaStartErr",RecThetaStartErr,"RecThetaStartErr[nTrack]/D");
	TrackTree->Branch("RecPhiStartErr",RecPhiStartErr,"RecPhiStartErr[nTrack]/D");
	TrackTree->Branch("RecXStartErr",RecXStartErr,"RecXStartErr[nTrack]/D");
	TrackTree->Branch("RecYStartErr",RecYStartErr,"RecYStartErr[nTrack]/D");
	TrackTree->Branch("RecTxStartErr",RecTxStartErr,"RecTxStartErr[nTrack]/D");
	TrackTree->Branch("RecTyStartErr",RecTyStartErr,"RecTyStartErr[nTrack]/D");
	TrackTree->Branch("nHit",&nHit,"nHit/I");
	TrackTree->Branch("HitX",HitX,"HitX[nHit]/D");
	TrackTree->Branch("HitY",HitY,"HitY[nHit]/D");
	TrackTree->Branch("HitZ",HitZ,"HitZ[nHit]/D");	
	TrackTree->Branch("HitTimeX",HitTimeX,"HitTimeX[nHit]/D");
	TrackTree->Branch("HitEdepX",HitEdepX,"HitEdepX[nHit]/D");
	TrackTree->Branch("HitTimeY",HitTimeY,"HitTimeY[nHit]/D");
	TrackTree->Branch("HitEdepY",HitEdepY,"HitEdepY[nHit]/D");
	TrackTree->Branch("IsMuonHit",IsMuonHit,"IsMuonHit[nHit]/I");
	TrackTree->Branch("HitPlaneId",HitPlaneId,"HitPlaneId[nHit]/I");
	TrackTree->Branch("nClusterHit",&nClusterHit,"nClusterHit/I");
	TrackTree->Branch("ClusterHitX",ClusterHitX,"ClusterHitX[nHit]/D");
	TrackTree->Branch("ClusterHitY",ClusterHitY,"ClusterHitY[nHit]/D");
	TrackTree->Branch("ClusterHitZ",ClusterHitZ,"ClusterHitZ[nHit]/D");
	TrackTree->Branch("nTrackHit",&nTrackHit,"nTrackHit/I");
	TrackTree->Branch("TrackHitX",TrackHitX,"TrackHitX[nTrackHit]/D");
	TrackTree->Branch("TrackHitY",TrackHitY,"TrackHitY[nTrackHit]/D");
	TrackTree->Branch("TrackHitZ",TrackHitZ,"TrackHitZ[nTrackHit]/D");
	TrackTree->Branch("TrackHitTimeX",TrackHitTimeX,"TrackHitTimeX[nTrackHit]/D");
	TrackTree->Branch("TrackHitEdepX",TrackHitEdepX,"TrackHitEdepX[nTrackHit]/D");
	TrackTree->Branch("TrackHitTimeY",TrackHitTimeY,"TrackHitTimeY[nTrackHit]/D");
	TrackTree->Branch("TrackHitEdepY",TrackHitEdepY,"TrackHitEdepY[nTrackHit]/D");
	TrackTree->Branch("TrackHitChi2",TrackHitChi2,"TrackHitChi2[nTrackHit]/D");
	TrackTree->Branch("TrackHitExpChi2",TrackHitExpChi2,"TrackHitExpChi2[nTrackHit]/D");
	TrackTree->Branch("TrackId",TrackId,"TrackId[nTrackHit]/I");
	TrackTree->Branch("nHitPlanes",&nHitPlanes,"nHitPlanes/I");
	TrackTree->Branch("TrackTxAtPlaneX",MuonTrackTxAtPlaneX,"TrackTxAtPlaneX[nHitPlanes]/D");
	TrackTree->Branch("TrackTyAtPlaneX",MuonTrackTyAtPlaneX,"TrackTyAtPlaneX[nHitPlanes]/D");
	TrackTree->Branch("TrackTxAtPlaneY",MuonTrackTxAtPlaneY,"TrackTxAtPlaneY[nHitPlanes]/D");
	TrackTree->Branch("TrackTyAtPlaneY",MuonTrackTyAtPlaneY,"TrackTyAtPlaneY[nHitPlanes]/D");
	TrackTree->Branch("IsMuon",&IsMuon,"IsMuon/I");

	//init data structures
	if(fReadSimInfoFromFile){//init from input file
		InitInfoFromFile();
	}
	
	//## init histo
	for(int i=0;i<nMaxPlanes;i++){
		int superplaneIndex= i/2;

		TString histo1Name;
		TString histo2Name;
		TString graph1Name;
		TString graph2Name;
		TString graph3Name;
		TString graph4Name;	
		int histoColor;
		if(i%2==0){
			histo1Name= Form("EnergyDeposit_Plane%dX",superplaneIndex+1);	
			histo2Name= Form("StripMultiplicity_Plane%dX",superplaneIndex+1);		
			graph1Name= Form("MuonDensity_Plane%dX",superplaneIndex+1);	
			graph2Name= Form("EmDensity_Plane%dX",superplaneIndex+1);	
			graph3Name= Form("MuonAlbedoDensity_Plane%dX",superplaneIndex+1);	
			graph4Name= Form("EmAlbedoDensity_Plane%dX",superplaneIndex+1);	
			histoColor= kGreen;
		}
		else{
			histo1Name= Form("EnergyDeposit_Plane%dY",superplaneIndex+1);
			histo2Name= Form("StripMultiplicity_Plane%dY",superplaneIndex+1);
			graph1Name= Form("MuonDensity_Plane%dY",superplaneIndex+1);	
			graph2Name= Form("EmDensity_Plane%dY",superplaneIndex+1);	
			graph3Name= Form("MuonAlbedoDensity_Plane%dY",superplaneIndex+1);	
			graph4Name= Form("EmAlbedoDensity_Plane%dY",superplaneIndex+1);		
			histoColor= kRed;
		}
		
		fEnergyDepositInPlaneHisto[i]= new TH1D(histo1Name,histo1Name,100,0,10);
		fEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
		fEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleSize(0.06);
		fEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleOffset(0.7);
		fEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitle("entries");
		fEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleSize(0.06);
		fEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleOffset(0.7);
		fEnergyDepositInPlaneHisto[i]->SetFillColor(histoColor);

		fElectronEnergyDepositInPlaneHisto[i]= new TH1D(TString("Electron")+histo1Name,TString("Electron")+histo1Name,100,0,10);
		fElectronEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
		fElectronEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleSize(0.06);
		fElectronEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleOffset(0.7);
		fElectronEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitle("entries");
		fElectronEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleSize(0.06);
		fElectronEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleOffset(0.7);
		fElectronEnergyDepositInPlaneHisto[i]->SetFillColor(histoColor);

		fMuonEnergyDepositInPlaneHisto[i]= new TH1D(TString("Muon")+histo1Name,TString("Muon")+histo1Name,100,0,10);
		fMuonEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
		fMuonEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleSize(0.06);
		fMuonEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleOffset(0.7);
		fMuonEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitle("entries");
		fMuonEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleSize(0.06);
		fMuonEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleOffset(0.7);
		fMuonEnergyDepositInPlaneHisto[i]->SetFillColor(histoColor);

		fGammaEnergyDepositInPlaneHisto[i]= new TH1D(TString("Gamma")+histo1Name,TString("Gamma")+histo1Name,100,0,10);
		fGammaEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
		fGammaEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleSize(0.06);
		fGammaEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleOffset(0.7);
		fGammaEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitle("entries");
		fGammaEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleSize(0.06);
		fGammaEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleOffset(0.7);
		fGammaEnergyDepositInPlaneHisto[i]->SetFillColor(histoColor);

		fProtonEnergyDepositInPlaneHisto[i]= new TH1D(TString("Proton")+histo1Name,TString("Proton")+histo1Name,100,0,10);
		fProtonEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
		fProtonEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleSize(0.06);
		fProtonEnergyDepositInPlaneHisto[i]->GetXaxis()->SetTitleOffset(0.7);
		fProtonEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitle("entries");
		fProtonEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleSize(0.06);
		fProtonEnergyDepositInPlaneHisto[i]->GetYaxis()->SetTitleOffset(0.7);
		fProtonEnergyDepositInPlaneHisto[i]->SetFillColor(histoColor);
				
		fStripMultiplicityInPlaneHisto[i]= new TH1D(histo2Name,histo2Name,21,-0.5,20.5);
		fStripMultiplicityInPlaneHisto[i]->GetXaxis()->SetTitle("Multiplicity");
		fStripMultiplicityInPlaneHisto[i]->GetXaxis()->SetTitleSize(0.06);
		fStripMultiplicityInPlaneHisto[i]->GetXaxis()->SetTitleOffset(0.7);
		fStripMultiplicityInPlaneHisto[i]->GetYaxis()->SetTitle("entries");
		fStripMultiplicityInPlaneHisto[i]->GetYaxis()->SetTitleSize(0.06);
		fStripMultiplicityInPlaneHisto[i]->GetYaxis()->SetTitleOffset(0.7);
		fStripMultiplicityInPlaneHisto[i]->SetFillColor(histoColor);

		fMuonDensityInPlaneGraph[i]= new TGraphErrors;
		fMuonDensityInPlaneGraph[i]->SetName(graph1Name);
		fMuonDensityInPlaneGraph[i]->GetXaxis()->SetTitle("r [m]");
		fMuonDensityInPlaneGraph[i]->GetXaxis()->SetTitleSize(0.06);
		fMuonDensityInPlaneGraph[i]->GetXaxis()->SetTitleOffset(0.7);
		fMuonDensityInPlaneGraph[i]->GetYaxis()->SetTitle("#rho_{#mu} [m^{-2}]");
		fMuonDensityInPlaneGraph[i]->GetYaxis()->SetTitleSize(0.06);
		fMuonDensityInPlaneGraph[i]->GetYaxis()->SetTitleOffset(0.7);
		fMuonDensityInPlaneGraph[i]->SetMarkerColor(kBlack);
		fMuonDensityInPlaneGraph[i]->SetMarkerStyle(8);
		fMuonDensityInPlaneGraph[i]->SetLineColor(kBlack);

		fMuonDensityAlbedoInPlaneGraph[i]= new TGraphErrors;
		fMuonDensityAlbedoInPlaneGraph[i]->SetName(graph3Name);
		fMuonDensityAlbedoInPlaneGraph[i]->GetXaxis()->SetTitle("r [m]");
		fMuonDensityAlbedoInPlaneGraph[i]->GetXaxis()->SetTitleSize(0.06);
		fMuonDensityAlbedoInPlaneGraph[i]->GetXaxis()->SetTitleOffset(0.7);
		fMuonDensityAlbedoInPlaneGraph[i]->GetYaxis()->SetTitle("#rho_{#mu} [m^{-2}]");
		fMuonDensityAlbedoInPlaneGraph[i]->GetYaxis()->SetTitleSize(0.06);
		fMuonDensityAlbedoInPlaneGraph[i]->GetYaxis()->SetTitleOffset(0.7);
		fMuonDensityAlbedoInPlaneGraph[i]->SetMarkerColor(kGray+1);
		fMuonDensityAlbedoInPlaneGraph[i]->SetMarkerStyle(8);
		fMuonDensityAlbedoInPlaneGraph[i]->SetLineColor(kGray+1);

		nMuonDensityPointsCounter[i]= 0;
		nMuonDensityAlbedoPointsCounter[i]= 0;
		
		fEmDensityInPlaneGraph[i]= new TGraphErrors;
		fEmDensityInPlaneGraph[i]->SetName(graph2Name);
		fEmDensityInPlaneGraph[i]->GetXaxis()->SetTitle("r [m]");
		fEmDensityInPlaneGraph[i]->GetXaxis()->SetTitleSize(0.06);
		fEmDensityInPlaneGraph[i]->GetXaxis()->SetTitleOffset(0.7);
		fEmDensityInPlaneGraph[i]->GetYaxis()->SetTitle("#rho_{em} [m^{-2}]");
		fEmDensityInPlaneGraph[i]->GetYaxis()->SetTitleSize(0.06);
		fEmDensityInPlaneGraph[i]->GetYaxis()->SetTitleOffset(0.7);
		fEmDensityInPlaneGraph[i]->SetMarkerColor(kBlack);
		fEmDensityInPlaneGraph[i]->SetMarkerStyle(8);
		fEmDensityInPlaneGraph[i]->SetLineColor(kBlack);

		fEmDensityAlbedoInPlaneGraph[i]= new TGraphErrors;
		fEmDensityAlbedoInPlaneGraph[i]->SetName(graph4Name);
		fEmDensityAlbedoInPlaneGraph[i]->GetXaxis()->SetTitle("r [m]");
		fEmDensityAlbedoInPlaneGraph[i]->GetXaxis()->SetTitleSize(0.06);
		fEmDensityAlbedoInPlaneGraph[i]->GetXaxis()->SetTitleOffset(0.7);
		fEmDensityAlbedoInPlaneGraph[i]->GetYaxis()->SetTitle("#rho_{em} [m^{-2}]");
		fEmDensityAlbedoInPlaneGraph[i]->GetYaxis()->SetTitleSize(0.06);
		fEmDensityAlbedoInPlaneGraph[i]->GetYaxis()->SetTitleOffset(0.7);
		fEmDensityAlbedoInPlaneGraph[i]->SetMarkerColor(kGray+1);
		fEmDensityAlbedoInPlaneGraph[i]->SetMarkerStyle(8);
		fEmDensityAlbedoInPlaneGraph[i]->SetLineColor(kGray+1);

		nEmDensityPointsCounter[i]= 0;
		nEmDensityAlbedoPointsCounter[i]= 0;
	}

	fTrigger3FoldEnergyHisto= new TH1D("Trigger3FoldEnergy","Trigger3FoldEnergy",100,-2,5);
	fTrigger3FoldEnergyHisto->GetXaxis()->SetTitle("Energy [GeV]");
	fTrigger3FoldEnergyHisto->GetXaxis()->SetTitleSize(0.06);
	fTrigger3FoldEnergyHisto->GetXaxis()->SetTitleOffset(0.7);
	fTrigger3FoldEnergyHisto->GetYaxis()->SetTitle("entries");
	fTrigger3FoldEnergyHisto->GetYaxis()->SetTitleSize(0.06);
	fTrigger3FoldEnergyHisto->GetYaxis()->SetTitleOffset(0.7);

	fTrigger2FoldEnergyHisto= new TH1D("Trigger2FoldEnergy","Trigger2FoldEnergy",100,-2,5);
	fTrigger2FoldEnergyHisto->GetXaxis()->SetTitle("Energy [GeV]");
	fTrigger2FoldEnergyHisto->GetXaxis()->SetTitleSize(0.06);
	fTrigger2FoldEnergyHisto->GetXaxis()->SetTitleOffset(0.7);
	fTrigger2FoldEnergyHisto->GetYaxis()->SetTitle("entries");
	fTrigger2FoldEnergyHisto->GetYaxis()->SetTitleSize(0.06);
	fTrigger2FoldEnergyHisto->GetYaxis()->SetTitleOffset(0.7);
	

	//currentSimParticle= new TParticleSimData;

  return eSuccess;
}


void G4MuonCounterReconstructor::InitInfoFromFile(){
 
	gROOT->ProcessLine(" gSystem->Load(\"libTHits.so\")");
	
	fInputFile= new TFile(fInputFileName.c_str(),"READ");

	if ( fInputFile->IsZombie() || !fInputFile || !fInputFile->IsOpen() ) {
		ostringstream errinfo;
  	errinfo << "Error opening sim data file " << fInputFileName;         
  	ERROR(errinfo.str().c_str());
		//return eFailure;
		exit(1);
  }

	ostringstream info;
  info << "Reading sim data from file " << fInputFileName;         
  INFO(info.str().c_str());
	
	GenTree= (TTree*)fInputFile->Get("GenInfo");
	DetTree= (TTree*)fInputFile->Get("DetInfo");
	SimTree= (TTree*)fInputFile->Get("SimInfo");
  
	if(!GenTree){
		ostringstream warninfo;
  	warninfo << "Cannot get tree with generated info from file "<< fInputFileName;         
  	WARNING(warninfo.str().c_str());	
	}

	if(!DetTree || !SimTree){
		ostringstream errinfo;
  	errinfo << "Cannot get tree with detector or simulation info from file "<< fInputFileName;         
  	ERROR(errinfo.str().c_str());	
		//return false;
		return;
	}

	//initialize vectors
	fDetectorData= 0;
	fEventSimData= 0;

	//## det info tree
	DetTree->SetBranchAddress("Detector",&fDetectorData); 

	//## gen info tree
	if(GenTree){	
		GenTree->SetBranchAddress("lgE",&lgE);
  	GenTree->SetBranchAddress("Theta",&Theta);
  	GenTree->SetBranchAddress("Phi",&Phi);
		GenTree->SetBranchAddress("Rp",&Rp);
  	GenTree->SetBranchAddress("Psi",&Psi);
  	GenTree->SetBranchAddress("GPSTimeSec",&GPSTimeSec);
		GenTree->SetBranchAddress("GPSTimeNanoSec",&GPSTimeNanoSec);
		GenTree->SetBranchAddress("CoreX",&CoreX);
		GenTree->SetBranchAddress("CoreY",&CoreY);
  	GenTree->SetBranchAddress("nMu",&nMu);
  	GenTree->SetBranchAddress("nMuPlus",&nMuPlus);
		GenTree->SetBranchAddress("nMuMinus",&nMuMinus);
		GenTree->SetBranchAddress("nE",&nE);
		GenTree->SetBranchAddress("nEPlus",&nEPlus);
		GenTree->SetBranchAddress("nEMinus",&nEMinus);
  	GenTree->SetBranchAddress("nGamma",&nGamma);
  	GenTree->SetBranchAddress("nProton",&nProton);
		GenTree->SetBranchAddress("nNeutron",&nNeutron);
		GenTree->SetBranchAddress("HadronicModel",&HadronicModel);
		GenTree->SetBranchAddress("PrimaryParticle",&PrimaryParticle);
		GenTree->SetBranchAddress("ShowerNumber",&ShowerNumber);
  	GenTree->SetBranchAddress("ShowerRunId",&ShowerRunId);
  }

	//## sim info tree
	SimTree->SetBranchAddress("EventSimData",&fEventSimData); 

	//return eSuccess;
	//return true;

}//close G4MuonCounterReconstructor::InitInfoFromFile()

VModule::ResultFlag
G4MuonCounterReconstructor::Run(Event& event)
{
	
	INFO("Running muon counter reconstruction");
	
	//## Run from simulation sequence or from ROOT file?
	(fReadSimInfoFromFile) ?
		RunFromFile() :
		RunFromSimulator(event);
			
  return eSuccess;

}//close G4MuonCounterReconstructor::Run()


VModule::ResultFlag
G4MuonCounterReconstructor::RunFromSimulator(Event& event)
{

	INFO("Retrieve sim info from simulator class");
	
	//## Retrieve sim info from simulator class
	fDetectorData= G4MuonCounterSimulator::GetCurrentMuonDetector();
	fEventSimData= G4MuonCounterSimulator::GetCurrentEventSimData();

	//## Get Detector info
	if(!fDetectorData){
		ERROR("Cannot get detector data from G4MuonCounterSimulator module!");
		return eFailure;
	}
	
	//## Get EventSimData
	if(!fEventSimData){
		ERROR("Cannot get event sim data from G4MuonCounterSimulator module!");
		return eFailure;
	}

	FillDetectorInfo();
	
	fStationSimDataCollection= fEventSimData->fStationSimDataCollection;

	unsigned int nSimStation= fStationSimDataCollection.size();
	cout<<"nSimStation="<<nSimStation<<endl;

	for(unsigned int k=0;k<nSimStation;k++){
		currentStationSimData= &(fStationSimDataCollection.at(k));
		fCurrentStation= currentStationSimData->fId;
		
		//## Init station info
		if(fMergeSimParticles && fProcessHits){
			//## Init station info
			InitStationHits();
		
			//## Find station trigger time
			FindStationStartTime();

			//## Analyze station data		
			ProcessStationHits();	
		}
		else{
			//## Loop over all particles
			std::vector<TParticleSimData> currentParticleSimDataCollection= currentStationSimData->fParticleSimDataCollection;
			unsigned int nSimParticle= currentParticleSimDataCollection.size();
			cout<<"nSimParticle="<<nSimParticle<<endl;
			for(unsigned int j=0;j<nSimParticle;j++){
				if(j%1000==0) cout<<"Reconstructing particle no. "<<j<<endl;
				currentSimParticle= &(currentParticleSimDataCollection.at(j));	
	
				//## Process particle hits?
				if(fProcessHits) ProcessParticleHits();
				
				//## Analyze particle data		
				FillSimParticleInfo();

			}//end loop sim particles
			cout<<endl;

		}//close else not merge particle

	}//end loop stations
	cout<<endl;

	return eSuccess;

}//close G4MuonCounterReconstructor::RunFromSimulator()


VModule::ResultFlag
G4MuonCounterReconstructor::RunFromFile(){

	INFO("Retrieve sim info from input ROOT file");
	
	//## Check Detector info
	if(!DetTree){
		ERROR("Cannot get detector info from input file");
		return eFailure;
	}

	//## Check Sim info
	if(!SimTree){
		ERROR("Cannot get tree with detector or simulation info from input file");	
		return eFailure;
	}	

	//## Check entries in sim file
	fNevents= SimTree->GetEntries();
	if(fNevents<=0){
		ostringstream warnMess;
  	warnMess << "No simulated events stored in file " << fInputFileName;         
  	WARNING(warnMess.str().c_str());
		return eContinueLoop;
		//return false;
	}
	
	ostringstream infoMess;
  infoMess << "Reading "<<fNevents<<" simulated events from input file " << fInputFileName;         
  INFO(infoMess.str().c_str());
	
	//## Get detector info
  DetTree->GetEntry(0);
	FillDetectorInfo();

	nMuonTracks= 0;
	
	//## Reconstruct events in input ROOT file
	for(int i=0;i<fNevents;i++){
		cout<<"==> EVENT no. "<<i+1<<endl;
		fCurrentEventNo= i;
		if(GenTree) GenTree->GetEntry(fCurrentEventNo);			
		SimTree->GetEntry(fCurrentEventNo);	

		fStationSimDataCollection= fEventSimData->fStationSimDataCollection;
		unsigned int nSimStation= fStationSimDataCollection.size();
		cout<<"nSimStations "<<nSimStation<<endl;
		//## Looping over all sim stations for this event
		for(unsigned int k=0;k<nSimStation;k++){
			cout<<"Station no. "<<k<<endl;
			currentStationSimData= &(fStationSimDataCollection.at(k));
			fCurrentStation= currentStationSimData->fId;

			//## Init station info
			if(fMergeSimParticles && fProcessHits){
				//## Init station info
				InitStationHits();
		
				//## Find station trigger time
				FindStationStartTime();

				//## Analyze station data		
				ProcessStationHits();	
			}
			else{
				//## Loop over all particles
				std::vector<TParticleSimData> currentParticleSimDataCollection= currentStationSimData->fParticleSimDataCollection;
				unsigned int nSimParticle= currentParticleSimDataCollection.size();
				cout<<"nSimParticle="<<nSimParticle<<endl;
				for(unsigned int j=0;j<nSimParticle;j++){
					if(j%1000==0) cout<<"Reconstructing particle no. "<<j<<endl;
					currentSimParticle= &(currentParticleSimDataCollection.at(j));	

					if(fProcessHits) ProcessParticleHits();
					//if(fMergeSimParticles) MergeParticleHits();//OLD
					//else FillSimParticleInfo();//OLD
		
					//## Analyze particle data		
					FillSimParticleInfo();

				}//end loop sim particles
				cout<<endl;

			}//close else not merge particles
			
		}//end loop stations
		//cout<<endl;
	}//close loop events in file


	return eSuccess;

}//close RunFromFile() 



void G4MuonCounterReconstructor::FillDetectorInfo(){


	//## Fill detector info
	nStrips= fDetectorData->GetNumberOfStrips();
	nPlanes= fDetectorData->GetNumberOfPlanes();
	double StripCoatSizeY= (fDetectorData->GetStripCoatingSize()).Y()/CLHEP::m;
	double StripLength= (fDetectorData->GetStripSize()).X()/CLHEP::m;
	double SensAreaSizeY= nStrips*StripCoatSizeY;
	double SensAreaSizeZ= StripLength;
	fDetectorAreaSize= pow(min(SensAreaSizeY,SensAreaSizeZ),2);//in m^2
	double StripSizeY= (fDetectorData->GetStripSize()).Y()/CLHEP::cm;
	double StripCoatHousingSizeY= fDetectorData->GetStripCoatingThickness()/CLHEP::cm;
	double SuperPlaneModuleSizeZ= fDetectorData->GetSuperPlaneHeight()/CLHEP::cm;
	double StripModuleSizeY= (StripSizeY + 2.* StripCoatHousingSizeY);
	fHitSeparationInCluster= 3./2.*StripModuleSizeY;
	//fHitSeparationInCluster= 2.*StripModuleSizeY;
	StripSeparation= StripModuleSizeY;

	fHitSigmaX= StripModuleSizeY/sqrt(12.);
	fHitSigmaY= StripModuleSizeY/sqrt(12.);
	fHitSigmaZ= SuperPlaneModuleSizeZ/sqrt(12.);
	cout<<"StripModuleSizeY="<<StripModuleSizeY<<"  StripSizeY="<<StripSizeY<<"  StripCoatThickness="<<StripCoatHousingSizeY<<endl;
	cout<<"StripSizeY="<<StripSizeY<<"  StripCoatHousingSizeY="<<StripCoatHousingSizeY<<endl;
	cout<<"StripModuleSizeY="<<StripModuleSizeY<<"  SuperPlaneModuleSizeZ="<<SuperPlaneModuleSizeZ<<endl; 
	cout<<"fHitSigmaX="<<fHitSigmaX<<"  fHitSigmaY="<<fHitSigmaY<<"  fHitSigmaZ="<<fHitSigmaZ<<"  fHitSeparationInCluster="<<fHitSeparationInCluster<<endl;

	fSuperPlaneDepth= fDetectorData->GetSuperPlaneDepth();
	fSuperPlaneSizeZ= fDetectorData->GetSuperPlaneHeight()/CLHEP::cm;
	cout<<"PlaneDepths"<<endl;
	for(unsigned int i=0;i<fSuperPlaneDepth.size();i++){
		cout<< "Plane "<<i<<" = "<<fSuperPlaneDepth[i]<< endl;
	}
	cout<<"fSuperPlaneSizeZ="<<fSuperPlaneSizeZ<<endl;

	//return true;
			
}//close G4MuonCounterReconstructor::FillDetectorInfo()



void G4MuonCounterReconstructor::FillSimParticleInfo(){

	//cout<<"G4MuonCounterReconstructor::FillSimParticleInfo()"<<endl;
	
	StationId= currentStationSimData->fId;
	RStat= currentStationSimData->fRadius;
	RStatGrd= currentStationSimData->fRadiusGrd;

	if(!currentSimParticle->fDetectorData) 
		currentSimParticle->fDetectorData= fDetectorData;

	Id= currentSimParticle->fId;
	std::string ParticleName= currentSimParticle->fName;		
	E= currentSimParticle->fEnergy;
	TVector3 particleTrueDirection= TVector3(-currentSimParticle->fDirection.X(),-currentSimParticle->fDirection.Y(),-currentSimParticle->fDirection.Z());
	Theta= particleTrueDirection.Theta()*180./TMath::Pi();
	Phi= particleTrueDirection.Phi()*180./TMath::Pi();
	Tx= particleTrueDirection.X()/particleTrueDirection.Z();
	Ty= particleTrueDirection.Y()/particleTrueDirection.Z();
	//Theta= currentSimParticle->fTheta;
	//Phi= currentSimParticle->fPhi;
	Time= currentSimParticle->fTime;
	X= (currentSimParticle->fPosition).X();
	Y= (currentSimParticle->fPosition).Y();
	R= currentSimParticle->fRadius;
	RGrd= currentSimParticle->fRadiusGrd;
	Psi= currentSimParticle->fPsi;
	PsiGrd= currentSimParticle->fPsiGrd;
	EventTag= currentSimParticle->fParticleSimDataTag;
	nHitPlanes= (int)((currentSimParticle->fStripMultiplicityX).size());
	
	double MaxRelAngle= CalculateMaxRelAngle(currentSimParticle);
	
	GeomEventTag= GeomEventTagger(currentSimParticle);

	cout<<"EventTag= "<<EventTag<<"  GeomEventTag="<<GeomEventTag<<"  HasGeomCrossedFirstPlane? "<<HasGeomCrossedFirstPlane<<endl;
	
	std::vector<int> currentStripMultiplicity= currentSimParticle->fStripMultiplicity;
	std::vector<int> currentStripMultiplicityX= currentSimParticle->fStripMultiplicityX;
	std::vector<int> currentStripMultiplicityY= currentSimParticle->fStripMultiplicityY;
	std::vector<double> currentStripAverageXHit= currentSimParticle->fAverageXHit;
	std::vector<double> currentStripAverageYHit= currentSimParticle->fAverageYHit;
	std::vector<double> currentStripWeightRX= currentSimParticle->fRX;		
	std::vector<double> currentStripWeightRY= currentSimParticle->fRY;

	
	//cout<<"Find trigger times per plane"<<endl;
	//## Find trigger times per plane
	//## Looping over all hit strips produced by this particle
	double TriggerTime[2*nHitPlanes];
	TVector3 MuonTrackDirection[2*nHitPlanes];

	for(int i=0;i<nHitPlanes;i++){
		TriggerTimeX[i]= 1.e+99; 
		TriggerTimeY[i]= 1.e+99;
		MuonTrackTxAtPlaneX[i]= 1.e+99;
		MuonTrackTxAtPlaneY[i]= 1.e+99;
		MuonTrackTyAtPlaneX[i]= 1.e+99;
		MuonTrackTyAtPlaneY[i]= 1.e+99;
	}
	for(int i=0;i<2*nHitPlanes;i++){
		TriggerTime[i]= 1.e+99;
		MuonTrackDirection[i]= TVector3(1.e+99,1.e+99,1.e+99);
	}
		
	std::vector<TScintHit> currentScintHitCollection= currentSimParticle->fSelScintHitCollection;
	//cout<<"currentScintHitCollection size="<<currentScintHitCollection.size()<<endl;
	for(unsigned int s=0;s<currentScintHitCollection.size();s++){
		currentScintHit= &(currentScintHitCollection.at(s));
		int StripId= currentScintHit->StripId;
		int PlaneId= currentScintHit->PlaneId;
		int SuperPlaneId= currentScintHit->SuperPlaneId;
		int AbsPlaneId= 2*SuperPlaneId+PlaneId;
		double TotEnergyDep= currentScintHit->Etot;

		//## Fill energy deposit histo
		fEnergyDepositInPlaneHisto[AbsPlaneId]->Fill(TotEnergyDep/CLHEP::MeV);
		if(abs(Id) == utl::Particle::eElectron)
			fElectronEnergyDepositInPlaneHisto[AbsPlaneId]->Fill(TotEnergyDep/CLHEP::MeV);
		else if(abs(Id) == utl::Particle::eMuon)
			fMuonEnergyDepositInPlaneHisto[AbsPlaneId]->Fill(TotEnergyDep/CLHEP::MeV);	
		else if(Id == utl::Particle::ePhoton)
			fGammaEnergyDepositInPlaneHisto[AbsPlaneId]->Fill(TotEnergyDep/CLHEP::MeV);	
		else if(Id == utl::Particle::eProton || Id == utl::Particle::eNeutron)
			fProtonEnergyDepositInPlaneHisto[AbsPlaneId]->Fill(TotEnergyDep/CLHEP::MeV);						

		//## Looping over all hits produced in this strip
		unsigned int nHitsForThisStrip= (currentScintHit->Edep).size();
		//cout<<"nHitsForThisStrip="<<nHitsForThisStrip<<endl;
		for(unsigned int ss=0;ss<nHitsForThisStrip;ss++){	
			double Edep= currentScintHit->Edep[ss];
			double Time= currentScintHit->Time[ss];			
			int TrackId= currentScintHit->TrackId[ss];
			int ParticleType= currentScintHit->ParticleType[ss];
			TVector3 TrackDirection= currentScintHit->TrackDirection[ss];
			
			//cout<<"Time "<<Time<<"  TriggerTime[AbsPlaneId]="<<TriggerTime[AbsPlaneId]<<endl;
			//cout<<"TrackId="<<TrackId<<"  TrackDirection=";
			//TrackDirection.Print();
			//cout<<endl;
			if(Time<TriggerTime[AbsPlaneId]) TriggerTime[AbsPlaneId]= Time;	
			if(TrackId==1 && abs(Id) == utl::Particle::eMuon) MuonTrackDirection[AbsPlaneId]= TrackDirection;	
		}//end loop nHits for this strip
	}//end loop strip hits for this particle
	
	for(int i=0;i<2*nHitPlanes;i++){
		int index= i/2;
		TVector3 oppositeTrackDir= TVector3(-MuonTrackDirection[i].X(),-MuonTrackDirection[i].Y(),-MuonTrackDirection[i].Z());
		if(i%2==0) {
			TriggerTimeX[index]= TriggerTime[i];
			MuonTrackTxAtPlaneX[i]= oppositeTrackDir.X()/oppositeTrackDir.Z(); 
			MuonTrackTyAtPlaneX[i]= oppositeTrackDir.Y()/oppositeTrackDir.Z(); 
		}
		else {
			TriggerTimeY[index]= TriggerTime[i];		
			MuonTrackTxAtPlaneY[i]= oppositeTrackDir.X()/oppositeTrackDir.Z(); 
			MuonTrackTyAtPlaneY[i]= oppositeTrackDir.Y()/oppositeTrackDir.Z(); 	
		}	
		//cout<<"Plane "<<i<<"  TriggerTime="<<TriggerTime[i]<<endl;
	}
	

	//## Fill average position, multiplicity, ...
	for(int i=0;i<nHitPlanes;i++){
		AverageXHit[i]= (currentSimParticle->fAverageXHit).at(i);
		AverageYHit[i]= (currentSimParticle->fAverageYHit).at(i);
		MultiplicityX[i]= (currentSimParticle->fStripMultiplicityX).at(i);
		MultiplicityY[i]= (currentSimParticle->fStripMultiplicityY).at(i);	
		RX[i]= (currentSimParticle->fRX).at(i);
		RY[i]= (currentSimParticle->fRY).at(i);
		//cout<<"plane no. "<<i<<"  MultiplicityX,Y="<<MultiplicityX[i]<<"  "<<MultiplicityY[i]<<endl;
	}



	//## Fill multiplicity histo
	for(unsigned t=0;t<currentStripMultiplicity.size();t++)
		fStripMultiplicityInPlaneHisto[t]->Fill(currentStripMultiplicity[t]);
			
	
	//## Fill number of particles per plane
	for(int i=0;i<2*nHitPlanes;i++){
		int index= (int)(i/2);
		if(i%2==0){
			nMuonsX[index]= (currentSimParticle->fMuonNumber).at(i);
			nMuonsAlbedoX[index]= (currentSimParticle->fMuonAlbedoNumber).at(i);
	
			nEmX[index]= (currentSimParticle->fEmNumber).at(i);
			nEmAlbedoX[index]= (currentSimParticle->fEmAlbedoNumber).at(i);

			nHadronsX[index]= (currentSimParticle->fHadronNumber).at(i);
			nHadronsAlbedoX[index]= (currentSimParticle->fHadronAlbedoNumber).at(i);
		}
		else{
			nMuonsY[index]= (currentSimParticle->fMuonNumber).at(i);
			nMuonsAlbedoY[index]= (currentSimParticle->fMuonAlbedoNumber).at(i);

			nEmY[index]= (currentSimParticle->fEmNumber).at(i);
			nEmAlbedoY[index]= (currentSimParticle->fEmAlbedoNumber).at(i);

			nHadronsY[index]= (currentSimParticle->fHadronNumber).at(i);
			nHadronsAlbedoY[index]= (currentSimParticle->fHadronAlbedoNumber).at(i);
		}
	}

	lgE= log10(E*1.e-3);//convert in GeV
	if(EventTag==3) fTrigger3FoldEnergyHisto->Fill(lgE);
	if(EventTag>=2) fTrigger2FoldEnergyHisto->Fill(lgE);


	//## Fill tree
	RecTree->Fill();

	//## Check tracking conditions
	if(!fTrackMuons){
		ostringstream infoMess;
  	infoMess << "Tracking option not selected...skip tracking!";         
  	INFO(infoMess.str().c_str());
		//return true;
		return;
	}

	if(EventTag<TStationSimData::eThreeFold){
		ostringstream infoMess;
  	infoMess << "No 3-Fold events for this station...skip tracking!";         
  	INFO(infoMess.str().c_str());
		//return true;
		return;
	}



	//## Track reconstructor (make sense for muon tracks & at least 3-fold events)
	ThetaRec= -999;
	PhiRec= -999;


	//## Do particle identification?
	bool IsMuonInPID= false;
	if(fUsePID && EventTag>=TParticleSimData::eThreeFold){
		// get parameters for this particle
		double m1X= MultiplicityX[0];
		double m1Y= MultiplicityY[0];
		double m2X= MultiplicityX[1];
		double m2Y= MultiplicityY[1];
		double m3X= MultiplicityX[2];
		double m3Y= MultiplicityY[2];
		double r1X= RX[0];
		double r1Y= RY[0];
		double r2X= RX[1];
		double r2Y= RY[1];
		double r3X= RX[2];
		double r3Y= RY[2];

		cout<<"m1X="<<m1X<<"  m1Y="<<m1Y<<"  m2X="<<m2X<<"  m2Y="<<m2Y<<"  m3X="<<m3X<<"  m3Y="<<m3Y<<endl;
		cout<<"r1X="<<r1X<<"  r1Y="<<r1Y<<"  r2X="<<r2X<<"  r2Y="<<r2Y<<"  r3X="<<r3X<<"  r3Y="<<r3Y<<endl;

		//## normalize inputs
		double m_min= 0.;
		double m1X_max= 33.;
		double m1Y_max= 38.;
		double m2X_max= 31.;
		double m2Y_max= 35.;
		double m3X_max= 20.;
		double m3Y_max= 18.;
		double r_min= 0.;
		double r_max= 1.;

		// normalize in range [a',b'] from range [a,b]
		// formula: x' = ((x - a) / (b - a)) * (b' - a') + a'
		double a_new= -1.;
		double b_new= 1.;
		double m1X_norm= ((m1X-m_min)/(m1X_max-m_min)) * (b_new-a_new) + a_new;
		double m1Y_norm= ((m1Y-m_min)/(m1Y_max-m_min)) * (b_new-a_new) + a_new;
		double m2X_norm= ((m2X-m_min)/(m2X_max-m_min)) * (b_new-a_new) + a_new;
		double m2Y_norm= ((m2Y-m_min)/(m2Y_max-m_min)) * (b_new-a_new) + a_new;
		double m3X_norm= ((m3X-m_min)/(m3X_max-m_min)) * (b_new-a_new) + a_new;
		double m3Y_norm= ((m3Y-m_min)/(m3Y_max-m_min)) * (b_new-a_new) + a_new;

		double r1X_norm= ((r1X-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;
		double r1Y_norm= ((r1Y-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;	
		double r2X_norm= ((r2X-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;
		double r2Y_norm= ((r2Y-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;
		double r3X_norm= ((r3X-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;
		double r3Y_norm= ((r3Y-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;

		//cout<<endl;
		//cout<<"m1X_norm="<<m1X_norm<<"  m1Y_norm="<<m1Y_norm<<"  m2X_norm="<<m2X_norm<<"  m2Y_norm="<<m2Y_norm<<"  m3X_norm="<<m3X_norm<<"  m3Y_norm="<<m3Y_norm<<endl;
		//cout<<"r1X_norm="<<r1X_norm<<"  r1Y_norm="<<r1Y_norm<<"  r2X_norm="<<r2X_norm<<"  r2Y_norm="<<r2Y_norm<<"  r3X_norm="<<r3X_norm<<"  r3Y_norm="<<r3Y_norm<<endl;
		//cout<<endl;

		const int nParsNN= 12;
		std::vector<double> NNInputData(nParsNN,0);
		
		NNInputData[0]= m1X_norm;
		NNInputData[1]= m1Y_norm;
		NNInputData[2]= m2X_norm;
		NNInputData[3]= m2Y_norm;
		NNInputData[4]= m3X_norm;
		NNInputData[5]= m3Y_norm;
		NNInputData[6]= r1X_norm;
		NNInputData[7]= r1Y_norm;
		NNInputData[8]= r2X_norm;
		NNInputData[9]= r2Y_norm;
		NNInputData[10]= r3X_norm;
		NNInputData[11]= r3Y_norm;

		IsMuonInPID= IsMuonEvent(NNInputData);
		//cout<<"IsMuonInPID? "<<IsMuonInPID<<endl;
	}
	else{
		if(abs(Id)==utl::Particle::eMuon) 
			IsMuonInPID= true;	
	}

	IsMuon= 0;
	if(abs(Id)==utl::Particle::eMuon){ 
		IsMuon= 1;	
	}
	
	nTrueTrack= 0;
	nTrueMuonTrack= 0;

	if(!IsMuonInPID){
		ostringstream infoMess;
  	infoMess << "No muons in this event...skip tracking!";         
  	INFO(infoMess.str().c_str());
		//return true;
		return;
	}


	//################################
	//##      TRACKING
	//###############################
	cout<<"*** Muon Track Reco ***"<<endl;

	std::vector<TrackPoint> currentTrackPointCollection= currentSimParticle->fTrackPointCollection;
			
	//## Get true muon track info
		
	GenE[nTrueTrack]= E;
	GenTheta[nTrueTrack]= Theta;
	GenPhi[nTrueTrack]= Phi;
	GenX[nTrueTrack]= X;
	GenY[nTrueTrack]= Y;
	GenTx[nTrueTrack]= Tx;
	GenTy[nTrueTrack]= Ty;
	GenR[nTrueTrack]= R;
	GenMaxRelAngle[nTrueTrack]= MaxRelAngle;

	nTrueTrack++;
	nTrueMuonTrack++;
				
	std::vector<Track*> RecTrackCollection;
	std::vector<TrackPoint> ClusterTrackPointCollection;
	std::vector<TrackPoint> FinalTrackPointCollection;
		
	TrackRecoStatus= -1;

	if(fTrackingAlgo == eLeastSquare){ 
		TrackFinder theTrackFinder;
		theTrackFinder.TrackWithMSCorrection(fTrackWithMSCorrection);
		theTrackFinder.UseEnergyLossCorrectionInMS(fUseEnergyLossCorrectionInMS);
		theTrackFinder.UseAverageEnergyInMS(fUseAverageEnergyInMS);
		theTrackFinder.SetAverageEnergyInMS(fAverageEnergyInMS/utl::MeV);
		theTrackFinder.UseTrueEnergyInMS(fUseTrueEnergyInMS);
		theTrackFinder.IncludeMomentumInTracking(fIncludeMomentumInTracking);
		
		theTrackFinder.UseClusterInTracking(fUseClusterInTracking);
		theTrackFinder.SplitClustersInTracking(fSplitClustersInTracking);
		theTrackFinder.SetSplitClusterThreshold(fSplitClusterThreshold);
		theTrackFinder.RemoveFakeHitsInTracking(fRemoveFakeHitsInTracking);
		
		theTrackFinder.SetTrackPoints(currentTrackPointCollection);
		theTrackFinder.SetHitSeparationInCluster(fHitSeparationInCluster);
		theTrackFinder.SetHitSigmaX(fHitSigmaX);
		theTrackFinder.SetHitSigmaY(fHitSigmaY);
		theTrackFinder.SetHitSigmaZ(fHitSigmaZ);
		theTrackFinder.SetSuperPlaneDepth(fSuperPlaneDepth);
		theTrackFinder.SetSuperPlaneSizeZ(fSuperPlaneSizeZ);
		theTrackFinder.SetVerbosity(1);
		theTrackFinder.SetGenEnergy(E);
		theTrackFinder.SetGenTheta(Theta);
		theTrackFinder.SetGenPhi(Phi);
		theTrackFinder.SetGenTx(Tx);
		theTrackFinder.SetGenTx(Ty);
		theTrackFinder.SetGenVertexX(X*100.);
		theTrackFinder.SetGenVertexY(Y*100.);
		theTrackFinder.SetRunNumber(nMuonTracks);
		TrackRecoStatus= theTrackFinder.TrackReconstructor();

		RecTrackCollection= theTrackFinder.GetRecTrackCollection();
		ClusterTrackPointCollection= theTrackFinder.GetClusterTrackPoints();
		FinalTrackPointCollection= theTrackFinder.GetTrackPoints();
	}
	else if(fTrackingAlgo == eKalmanFilter){
		KFTrackFinder theKFTrackFinder;
		theKFTrackFinder.TrackWithMSCorrection(fTrackWithMSCorrection);
		theKFTrackFinder.UseEnergyLossCorrectionInMS(fUseEnergyLossCorrectionInMS);
		theKFTrackFinder.UseAverageEnergyInMS(fUseAverageEnergyInMS);
		theKFTrackFinder.SetAverageEnergyInMS(fAverageEnergyInMS/utl::MeV);
		theKFTrackFinder.UseTrueEnergyInMS(fUseTrueEnergyInMS);
		theKFTrackFinder.IncludeMomentumInTracking(fIncludeMomentumInTracking);
			
		theKFTrackFinder.UseClusterInTracking(fUseClusterInTracking);
		theKFTrackFinder.SplitClustersInTracking(fSplitClustersInTracking);
		theKFTrackFinder.SetSplitClusterThreshold(fSplitClusterThreshold);
		theKFTrackFinder.RemoveFakeHitsInTracking(fRemoveFakeHitsInTracking);
			
		theKFTrackFinder.SetTrackPoints(currentTrackPointCollection);
		theKFTrackFinder.SetHitSeparationInCluster(fHitSeparationInCluster);
		theKFTrackFinder.SetHitSigmaX(fHitSigmaX);
		theKFTrackFinder.SetHitSigmaY(fHitSigmaY);
		theKFTrackFinder.SetHitSigmaZ(fHitSigmaZ);
		theKFTrackFinder.SetSuperPlaneDepth(fSuperPlaneDepth);
		theKFTrackFinder.SetSuperPlaneSizeZ(fSuperPlaneSizeZ);
		theKFTrackFinder.SetGenEnergy(E);
		theKFTrackFinder.SetGenTheta(Theta);
		theKFTrackFinder.SetGenPhi(Phi);
		theKFTrackFinder.SetGenTx(Tx);
		theKFTrackFinder.SetGenTx(Ty);
		theKFTrackFinder.SetGenVertexX(X*100.);
		theKFTrackFinder.SetGenVertexY(Y*100.);
		theKFTrackFinder.SetVerbosity(1);
		theKFTrackFinder.SetRunNumber(nMuonTracks);
		TrackRecoStatus= theKFTrackFinder.TrackReconstructor();

		RecTrackCollection= theKFTrackFinder.GetRecTrackCollection();
		ClusterTrackPointCollection= theKFTrackFinder.GetClusterTrackPoints();
		FinalTrackPointCollection= theKFTrackFinder.GetTrackPoints();
	}
	else{
		cerr<<"G4MuonCounterReconstructor::FillSimParticleInfo(): Invalid tracking algorithm...exit"<<endl;
		exit(1);
	}

	cout<<"***********   GEN PARTICLE INFO  ************"<<endl;
	cout<<"GenTheta [deg]="<< Theta <<"  GenPhi [deg]="<<Phi<<endl;
	cout<<"GenX [cm]="<< X*100. <<"  GenY [cm]="<<Y*100.<<endl;
	cout<<"GenTx="<< Tx <<"  GenTy="<< Ty<<endl;
	cout<<"*********************************************"<<endl;
		
	/*	
	for(int i=0;i<nHitPlanes;i++){
		cout<<"TxAtPlaneX="<<MuonTrackTxAtPlaneX[i]<<endl;
		cout<<"TyAtPlaneX="<<MuonTrackTyAtPlaneX[i]<<endl;
		cout<<"TxAtPlaneY="<<MuonTrackTxAtPlaneY[i]<<endl;
		cout<<"TyAtPlaneY="<<MuonTrackTyAtPlaneY[i]<<endl;
	}
	*/

	//## Get rec muon track info and cluster points
	nTrack= (int)(RecTrackCollection.size());
	nClusterHit= (int)(ClusterTrackPointCollection.size());
	nTrackHit= 0;
	cout<<"INFO: Tracking(): nTrack true/rec="<<nTrueTrack<<"/"<<nTrack<<endl;
	
	for(int s=0;s<nTrack;s++){
		Track* currentRecTrack= RecTrackCollection[s];
			
		TVector3 RecDirection= currentRecTrack->GetDirection();
		RecOpeningAngle[s]= RecDirection.Angle(particleTrueDirection)*180./TMath::Pi();

		RecTheta[s]= currentRecTrack->GetTheta();
		RecPhi[s]= currentRecTrack->GetPhi();
		RecX[s]= (currentRecTrack->GetVertexPos()).X();
		RecY[s]= (currentRecTrack->GetVertexPos()).Y();
		RecTx[s]= currentRecTrack->GetTx();
		RecTy[s]= currentRecTrack->GetTy();
		RecThetaErr[s]= currentRecTrack->GetThetaErr();
		RecPhiErr[s]= currentRecTrack->GetPhiErr();
		RecXErr[s]= (currentRecTrack->GetVertexPosErr()).X();
		RecYErr[s]= (currentRecTrack->GetVertexPosErr()).Y();
		RecTxErr[s]= currentRecTrack->GetTxErr();
		RecTyErr[s]= currentRecTrack->GetTyErr();
		TrackFitStatus[s]= currentRecTrack->GetFitStatus();	
		TrackFitChi2[s]= currentRecTrack->GetFitChi2();
		RecThetaStart[s]= currentRecTrack->GetThetaStart();
		RecPhiStart[s]= currentRecTrack->GetPhiStart();
		RecXStart[s]= (currentRecTrack->GetVertexPosStart()).X();
		RecYStart[s]= (currentRecTrack->GetVertexPosStart()).Y();
		RecTxStart[s]= currentRecTrack->GetTxStart();
		RecTyStart[s]= currentRecTrack->GetTyStart();
		
		RecThetaStartErr[s]= currentRecTrack->GetThetaStartErr();
		RecPhiStartErr[s]= currentRecTrack->GetPhiStartErr();
		RecXStartErr[s]= (currentRecTrack->GetVertexPosStartErr()).X();
		RecYStartErr[s]= (currentRecTrack->GetVertexPosStartErr()).Y();
		RecTxStartErr[s]= currentRecTrack->GetTxStartErr();
		RecTyStartErr[s]= currentRecTrack->GetTyStartErr();
	
		double x,y,z,t;
		TrackPoint aTrackPoint;
		int id= currentRecTrack->GetId();
			
		for(int p=0;p<currentRecTrack->GetNpoints();p++){
			currentRecTrack->GetPoint(p,x,y,z,t);
			currentRecTrack->GetTrackPoint(p,aTrackPoint);
			TrackHitX[nTrackHit]= x;
			TrackHitY[nTrackHit]= y;
			TrackHitZ[nTrackHit]= z;
			TrackId[nTrackHit]= id;
			TrackHitEdepX[nTrackHit]= aTrackPoint.fEdepX;
			TrackHitEdepY[nTrackHit]= aTrackPoint.fEdepY;	
			TrackHitTimeX[nTrackHit]= aTrackPoint.fTimeX;
			TrackHitTimeY[nTrackHit]= aTrackPoint.fTimeY; 
				
			currentRecTrack->GetTrackHitChi2(p,TrackHitChi2[nTrackHit]);	
			currentRecTrack->GetTrackHitExpChi2(p,TrackHitExpChi2[nTrackHit]);
			nTrackHit++;
		}
	}//end loop rec tracks	

	for(int s=0;s<nClusterHit;s++){
		TrackPoint currentClusterPoint= ClusterTrackPointCollection[s];	
		ClusterHitX[s]= (currentClusterPoint.fPosition).X();
		ClusterHitY[s]= (currentClusterPoint.fPosition).Y();
		ClusterHitZ[s]= (currentClusterPoint.fPosition).Z();
		//cout<<"ClusterHitX[s]="<<ClusterHitX[s]<<"  ClusterHitY[s]="<<ClusterHitY[s]<<"  ClusterHitZ[s]="<<ClusterHitZ[s]<<endl;
	}

	if(nTrack>0){
		ThetaRec= RecTrackCollection[0]->GetTheta();
		PhiRec= RecTrackCollection[0]->GetPhi();		
	}

	//nHit= (int)(currentTrackPointCollection.size());//## ALL SELECTED HITS
	nHit= (int)(FinalTrackPointCollection.size());//## ALL HITS PASSED TO TRACKING STEP

	//cout<<"nHit="<<nHit<<endl;
	for(int i=0;i<nHit;i++){
		//HitX[i]= (currentTrackPointCollection[i].fPosition).X();
		//HitY[i]= (currentTrackPointCollection[i].fPosition).Y();
		//HitZ[i]= (currentTrackPointCollection[i].fPosition).Z();
		//IsMuonHit[i]= currentTrackPointCollection[i].fIsMuon;
		//HitPlaneId[i]= currentTrackPointCollection[i].fDetectorPlaneId;
		//HitEdepX[i]= currentTrackPointCollection[i].fEdepX;
		//HitEdepY[i]= currentTrackPointCollection[i].fEdepY;	
		//HitTimeX[i]= currentTrackPointCollection[i].fTimeX;
		//HitTimeY[i]= currentTrackPointCollection[i].fTimeY;

		HitX[i]= (FinalTrackPointCollection[i].fPosition).X();
		HitY[i]= (FinalTrackPointCollection[i].fPosition).Y();
		HitZ[i]= (FinalTrackPointCollection[i].fPosition).Z();
		IsMuonHit[i]= FinalTrackPointCollection[i].fIsMuon;
		HitPlaneId[i]= FinalTrackPointCollection[i].fDetectorPlaneId;
		HitEdepX[i]= FinalTrackPointCollection[i].fEdepX;
		HitEdepY[i]= FinalTrackPointCollection[i].fEdepY;	
		HitTimeX[i]= FinalTrackPointCollection[i].fTimeX;
		HitTimeY[i]= FinalTrackPointCollection[i].fTimeY;  
		//cout<<"HitX[i]="<<HitX[i]<<"  HitY[i]="<<HitY[i]<<"  HitZ[i]="<<HitZ[i]<<"  IsMuonHit[i]="<<IsMuonHit[i]<<"  HitPlaneId[i]="<<HitPlaneId[i]<<endl;
		//cout<<"HitEdepX[i]="<<HitEdepX[i]<<"  HitEdepY[i]="<<HitEdepY[i]<<"  HitTimeX[i]="<<HitTimeX[i]<<"  HitTimeY[i]="<<HitTimeY[i]<<endl;
	}
		
	//## Fill Track tree
	TrackTree->Fill();
		
	nMuonTracks++;

	//return true;
	
}//close G4MuonCounterReconstructor::FillSimParticleInfo()


void G4MuonCounterReconstructor::ProcessParticleHits(){

	//cout<<"G4MuonCounterReconstructor::ProcessParticleHits()"<<endl;
	//## Access to current un-selected scint hit collection
	std::vector<TScintHit> thisScintHitCollection;
	if(fProcessSelectedHits)
		thisScintHitCollection= currentSimParticle->fSelScintHitCollection;	
	else 
		thisScintHitCollection= currentSimParticle->fScintHitCollection;
	currentSimParticle->nStrips= nStrips;
	currentSimParticle->nPlanes= nPlanes;
	//cout<<"G4MuonCounterReconstructor::ProcessParticleHits(): StripSeparation="<<StripSeparation<<endl;
	currentSimParticle->StripSeparation= StripSeparation;

	currentSimParticle->fHitSigmaX= fHitSigmaX;
	currentSimParticle->fHitSigmaY= fHitSigmaY;
	currentSimParticle->fHitSigmaZ= fHitSigmaZ;

	std::vector<TScintHit> reselectedScintHitCollection;
	
	//## Resize track point collection
	(currentSimParticle->fTrackPointCollection).clear();
	(currentSimParticle->fTrackPointCollection).resize(0);		


	//## Re-select the hit collection 
	unsigned int nHitsForThisParticle= thisScintHitCollection.size();
	//cout<<"nHitsForThisParticle="<<nHitsForThisParticle<<"  fEnergyThreshold="<<fEnergyThreshold<<endl;
	
	for(unsigned int s=0;s<nHitsForThisParticle;s++){
		int StripId= (thisScintHitCollection.at(s)).StripId;
		int PlaneId= (thisScintHitCollection.at(s)).PlaneId;
		int SuperPlaneId= (thisScintHitCollection.at(s)).SuperPlaneId;
		int AbsPlaneId= 2*SuperPlaneId+PlaneId;
		double TotEnergyDep= (thisScintHitCollection.at(s)).Etot;
		//cout<<"TotEnergyDep="<<TotEnergyDep<<endl;

		//cout<<"G4MuonCounterReconstructor::ProcessParticleHits(): AbsPlaneId="<<AbsPlaneId<<"  StripId="<<StripId<<endl;

		if(TotEnergyDep<fEnergyThreshold/utl::MeV) continue;

		reselectedScintHitCollection.push_back(thisScintHitCollection.at(s));					
	}//end loop strip hits for this particle
	
	//## Assign the reselected hit to the current sim particle
	//## Re-call CreateTrackPoints
	currentSimParticle->fSelScintHitCollection= reselectedScintHitCollection;
	currentSimParticle->CreateTrackPoints();
	
	//cout<<"end"<<endl;

}//close ProcessParticleHits()


void G4MuonCounterReconstructor::InitStationHits(){

	cout<<"*** G4MuonCounterReconstructor::InitStationHits() ***"<<endl;

	currentStationSimData->nStrips= nStrips;
	currentStationSimData->nPlanes= nPlanes;
	currentStationSimData->fHitSigmaX= fHitSigmaX;
	currentStationSimData->fHitSigmaY= fHitSigmaY;
	currentStationSimData->fHitSigmaZ= fHitSigmaZ;

	//## Resize hit collection
	(currentStationSimData->fSelScintHitCollection).clear();
	(currentStationSimData->fSelScintHitCollection).resize(0);

	(currentStationSimData->fScintHitCollection).clear();
	(currentStationSimData->fScintHitCollection).resize(0);

	//## Resize track point collection
	(currentStationSimData->fTrackPointCollection).clear();
	(currentStationSimData->fTrackPointCollection).resize(0);	

	//## Resize list of hit strips
	(currentStationSimData->fListOfHitAbsStripId).clear();
	(currentStationSimData->fListOfHitAbsStripId).resize(0);

	(currentStationSimData->fSelListOfHitAbsStripId).clear();
	(currentStationSimData->fSelListOfHitAbsStripId).resize(0);		

}//close G4MuonCounterReconstructor::InitStationHits()


void G4MuonCounterReconstructor::MergeParticleHits(){

	cout<<"*** G4MuonCounterReconstructor::MergeParticleHits() ***"<<endl;

	//cout<<"PARTICLE: "<<currentSimParticle->fId<<"  Time[ns]="<<currentSimParticle->fTime<<endl;

	//## Access to current un-selected scint hit collection
	std::vector<TScintHit> thisScintHitCollection;
	TStationSimData::HitSearchMode whereToAdd;
	if(fProcessSelectedHits){
		thisScintHitCollection= currentSimParticle->fSelScintHitCollection;	
		whereToAdd= TStationSimData::eInSelectionCollection;
	}
	else{ 
		thisScintHitCollection= currentSimParticle->fScintHitCollection;
		whereToAdd= TStationSimData::eInFullCollection;
	}

	//######################################
	//##    MERGE ALL HITS IN STATION (not needed after bug correction!)
	//######################################
	//## Merging all hits for this station
	unsigned int nHitsForThisParticle= thisScintHitCollection.size();
	//cout<<"nHitsForThisParticle="<<nHitsForThisParticle<<endl;
	
	for(unsigned int s=0;s<nHitsForThisParticle;s++){
		int StripId= (thisScintHitCollection.at(s)).StripId;
		int PlaneId= (thisScintHitCollection.at(s)).PlaneId;
		int SuperPlaneId= (thisScintHitCollection.at(s)).SuperPlaneId;
		int AbsPlaneId= 2*SuperPlaneId+PlaneId;
		double TotEnergyDep= (thisScintHitCollection.at(s)).Etot;
		//cout<<"TotEnergyDep="<<TotEnergyDep<<endl;
		if(TotEnergyDep<fEnergyThreshold/utl::MeV) continue;

		unsigned int nHitsForThisStrip= ((thisScintHitCollection.at(s)).Time).size();
		for(unsigned int j=0;j<nHitsForThisStrip;j++){
			double Time= (thisScintHitCollection.at(s)).Time[j]; 
			double Edep= (thisScintHitCollection.at(s)).Edep[j]; 
			int TrackId= (thisScintHitCollection.at(s)).TrackId[j]; 
			int ParticleType= (thisScintHitCollection.at(s)).ParticleType[j]; 
			cout<<"SP="<<SuperPlaneId<<"  P="<<PlaneId<<"  S="<<StripId<<"  edep[MeV]="<<Edep<< "  Time[ns]="<<Time<<"  TrackId="<<TrackId<<"  PartType="<<ParticleType<<endl; 
		}

		currentStationSimData->AddScintHit(thisScintHitCollection.at(s), whereToAdd);		
	}//end loop strip hits for this particle


}//close G4MuonCounterReconstructor::MergeParticleHits()


void G4MuonCounterReconstructor::FindStationStartTime(){
		
	cout<<"G4MuonCounterReconstructor::FindStationStartTime()"<<endl;
	//######################################
	//##    FIND TRIGGER START TIME 
	//######################################
	
	//## Looping over all sim particles for this station	
	std::vector<TParticleSimData> currentParticleSimDataCollection= currentStationSimData->fParticleSimDataCollection;
	unsigned int nSimParticle= currentParticleSimDataCollection.size();
	
	cout<<"==> "<<nSimParticle<<"  particles for this station"<<endl;

	std::vector<double> listOfHitTimes;
	listOfHitTimes.clear();
	listOfHitTimes.resize(0);

	for(unsigned int j=0;j<nSimParticle;j++){
		currentSimParticle= &(currentParticleSimDataCollection.at(j));
		
		//access to current particle hit collection
		std::vector<TScintHit> thisScintHitCollection= currentSimParticle->fScintHitCollection;

		unsigned int nHitsForThisParticle= thisScintHitCollection.size();
		//cout<<"nHitsForThisParticle="<<nHitsForThisParticle<<endl;
	

		for(unsigned int s=0;s<nHitsForThisParticle;s++){
			//int StripId= (thisScintHitCollection.at(s)).StripId;
			//int PlaneId= (thisScintHitCollection.at(s)).PlaneId;
			//int SuperPlaneId= (thisScintHitCollection.at(s)).SuperPlaneId;
			//int AbsPlaneId= 2*SuperPlaneId+PlaneId;
			double TotEnergyDep= (thisScintHitCollection.at(s)).Etot;
		
			if(TotEnergyDep<fEnergyThreshold) continue;

			//fill hit time list
			unsigned int nHitsForThisStrip= ((thisScintHitCollection.at(s)).Time).size();
			for(unsigned int kk=0;kk<nHitsForThisStrip;kk++){
				listOfHitTimes.push_back( (thisScintHitCollection.at(s)).Time[kk] );			
			}
		
		}//end loop strip hits for this particle

	}//end loop sim particles
	cout<<endl;

	if( listOfHitTimes.empty() ){
		fStartTime= -999;
	}
	else{
		fStartTime= *( std::min_element( listOfHitTimes.begin(), listOfHitTimes.end() ) );
	}

	fNSamplingBins= (int)(fMergeTimeInterval/fSamplingTime);


}//close G4MuonCounterReconstructor::FindStationStartTime()



void G4MuonCounterReconstructor::ProcessStationHits(){
	
	//## Looping over all sim particles for this station	
	std::vector<TParticleSimData> currentParticleSimDataCollection= currentStationSimData->fParticleSimDataCollection;
	
	unsigned int nSimParticle= currentParticleSimDataCollection.size();

	TH1D* dummyTimeHisto= new TH1D("dummyTimeHisto","dummyTimeHisto",fNSamplingBins,fStartTime,fStartTime + fMergeTimeInterval);

	std::vector< std::vector<int> > ListHitAbsStripIdInTimeBin;
	std::vector< std::vector<TScintHit> > ScintHitInTimeBin;
	std::vector< std::vector<TParticleSimData> > ParticleSimDataCollectionInTimeBin;
	

	for(int k=0;k<fNSamplingBins;k++){
		ListHitAbsStripIdInTimeBin.push_back ( std::vector<int>() ); 
		ScintHitInTimeBin.push_back ( std::vector<TScintHit>() ); 
		ParticleSimDataCollectionInTimeBin.push_back ( std::vector<TParticleSimData>() ); 
	}

	int nHitsForThisStation= 0;

	for(unsigned int j=0;j<nSimParticle;j++){
		currentSimParticle= &(currentParticleSimDataCollection.at(j));
		
		//access to current particle hit collection
		std::vector<TScintHit> thisScintHitCollection= currentSimParticle->fScintHitCollection;

		unsigned int nHitsForThisParticle= thisScintHitCollection.size();
		//cout<<"nHitsForThisParticle="<<nHitsForThisParticle<<endl;
		
		std::vector<double> tmpTimeList;
		tmpTimeList.clear();
		tmpTimeList.resize(0);
		std::vector<double> tmpEdepList;
		tmpEdepList.clear();
		tmpEdepList.resize(0);

		//## Create a new TParticleSimData, reset hit collections
		TParticleSimData* aNewParticleSimData= new TParticleSimData;
		aNewParticleSimData->fId= currentSimParticle->fId;
		aNewParticleSimData->fName= currentSimParticle->fName;
		aNewParticleSimData->fGlobalPosition= currentSimParticle->fGlobalPosition;
		aNewParticleSimData->fPosition= currentSimParticle->fPosition;
		aNewParticleSimData->fRadius= currentSimParticle->fRadius;
		aNewParticleSimData->fPsi= currentSimParticle->fPsi;
		aNewParticleSimData->fRadiusGrd= currentSimParticle->fRadiusGrd;
		aNewParticleSimData->fPsiGrd= currentSimParticle->fPsiGrd;
		aNewParticleSimData->fDirection= currentSimParticle->fDirection;
		aNewParticleSimData->fMomentum= currentSimParticle->fMomentum;
		aNewParticleSimData->fEnergy= currentSimParticle->fEnergy;
		aNewParticleSimData->fTheta= currentSimParticle->fTheta;
		aNewParticleSimData->fPhi= currentSimParticle->fPhi;
		aNewParticleSimData->fTime= currentSimParticle->fTime;
		aNewParticleSimData->nStrips= nStrips;
	 	aNewParticleSimData->nPlanes= nPlanes;
		aNewParticleSimData->fHitSigmaX= fHitSigmaX;
		aNewParticleSimData->fHitSigmaY= fHitSigmaY;
		aNewParticleSimData->fHitSigmaZ= fHitSigmaZ;
		aNewParticleSimData->fPMTHitCollection= currentSimParticle->fPMTHitCollection;

		for(unsigned int s=0;s<nHitsForThisParticle;s++){
			int StripId= (thisScintHitCollection.at(s)).StripId;
			int PlaneId= (thisScintHitCollection.at(s)).PlaneId;
			int SuperPlaneId= (thisScintHitCollection.at(s)).SuperPlaneId;
			int AbsPlaneId= 2*SuperPlaneId+PlaneId;
			int AbsStripId= AbsPlaneId*nStrips+StripId;
			double TotEnergyDep= (thisScintHitCollection.at(s)).Etot;
			TVector3 StripPosition= (thisScintHitCollection.at(s)).StripPosition;
			std::vector<TVector3> PhotocathodeSurfPosList= (thisScintHitCollection.at(s)).PhotocathodeSurfacePosition;

		
			if(TotEnergyDep<fEnergyThreshold) continue;

			nHitsForThisStation++;

			TScintHit aNewHit;
			aNewHit.StripId= StripId;
			aNewHit.PlaneId= PlaneId;	
			aNewHit.SuperPlaneId= SuperPlaneId;	
			aNewHit.StripPosition= StripPosition;
			aNewHit.PhotocathodeSurfacePosition= PhotocathodeSurfPosList;	
		
	

			//fill hit time list
			unsigned int nHitsForThisStrip= ((thisScintHitCollection.at(s)).Time).size();
			
			int thisTimeBin= -1;

			for(unsigned int kk=0;kk<nHitsForThisStrip;kk++){
				double Time= (thisScintHitCollection.at(s)).Time[kk];
				double Edep= (thisScintHitCollection.at(s)).Edep[kk];
				int TrackId= (thisScintHitCollection.at(s)).TrackId[kk];
				int ParticleType= (thisScintHitCollection.at(s)).ParticleType[kk];
				TVector3 Position= (thisScintHitCollection.at(s)).Position[kk];	
				//TVector3 TrackDirection= (thisScintHitCollection.at(s)).TrackDirection[kk];
		
				thisTimeBin= dummyTimeHisto->FindBin(Time);
			
				#if ROOT_VERSION_CODE > ROOT_VERSION(5,22,0)
					if(dummyTimeHisto->IsBinOverflow(thisTimeBin) || dummyTimeHisto->IsBinUnderflow(thisTimeBin))
						continue;
				#else
					if(thisTimeBin==0 || thisTimeBin==dummyTimeHisto->GetNbinsX()+1){
						continue;	
				#endif

				tmpTimeList.push_back(Time);
				tmpEdepList.push_back(Edep);
				
				//add to the new hit
				(aNewHit.Edep).push_back(Edep);
				(aNewHit.Time).push_back(Time);
				(aNewHit.TrackId).push_back(TrackId);
				(aNewHit.ParticleType).push_back(ParticleType);
				(aNewHit.Position).push_back(Position);
				//(aNewHit.TrackDirection).push_back(TrackDirection);
				
				(aNewHit.Etot)+= Edep;	

				unsigned int nPhotonsForThisStrip= ((thisScintHitCollection.at(s)).PhotonTrackLength).size();
				for(unsigned int jj=0;jj<nPhotonsForThisStrip;jj++){
					double PhotonTrackLength= (thisScintHitCollection.at(s)).PhotonTrackLength[jj];
					double PhotonTrackDistance= (thisScintHitCollection.at(s)).PhotonTrackDistance[jj];
					double PhotonEmissionAngle= (thisScintHitCollection.at(s)).PhotonEmissionAngle[jj];
					double PhotonEmissionTime= (thisScintHitCollection.at(s)).PhotonEmissionTime[jj];
					double PhotonEmissionWavelength= (thisScintHitCollection.at(s)).PhotonEmissionWavelength[jj];
					int PhotonProcessType= (thisScintHitCollection.at(s)).PhotonProcessType[jj];

					//add for new hit
					(aNewHit.PhotonTrackLength).push_back(PhotonTrackLength);
					(aNewHit.PhotonTrackDistance).push_back(PhotonTrackDistance);	
					(aNewHit.PhotonEmissionAngle).push_back(PhotonEmissionAngle);
					(aNewHit.PhotonEmissionTime).push_back(PhotonEmissionTime);
					(aNewHit.PhotonEmissionWavelength).push_back(PhotonEmissionWavelength);
					(aNewHit.PhotonProcessType).push_back(PhotonProcessType);
				}//end loop photons for this strip

				//search if current abs strip id was already included
				std::vector<int>::const_iterator it= std::find(ListHitAbsStripIdInTimeBin[thisTimeBin-1].begin(), ListHitAbsStripIdInTimeBin[thisTimeBin-1].end(), AbsStripId);

				if(it == ListHitAbsStripIdInTimeBin[thisTimeBin-1].end()){
					//this strip does not exist yet in the list, add it!
					ScintHitInTimeBin[thisTimeBin-1].push_back(aNewHit);
					ListHitAbsStripIdInTimeBin[thisTimeBin-1].push_back(AbsStripId);
					
				}//close if create new hits
				else{
					//append to the list
					int collectionIndex= it-ListHitAbsStripIdInTimeBin[thisTimeBin-1].begin();
					ScintHitInTimeBin[thisTimeBin-1][collectionIndex].Append(aNewHit);
				}//close else append hits

			}//end loop hits for this strip

			aNewParticleSimData->AddScintHit(aNewHit);
			aNewParticleSimData->AddSelScintHit(aNewHit);

		}//end loop strip hits for this particle

		//find time bin of particle according to the minimum trigger time
		if(!tmpTimeList.empty()){
			int pos= std::min_element(tmpTimeList.begin(), tmpTimeList.end()) - tmpTimeList.begin();		
			double TimeOfMaxEdep= *(tmpTimeList.begin() + pos);
			int thisParticleTimeBin= dummyTimeHisto->FindBin(TimeOfMaxEdep);
			
			#if ROOT_VERSION_CODE > ROOT_VERSION(5,22,0)
				if(dummyTimeHisto->IsBinOverflow(thisParticleTimeBin) || dummyTimeHisto->IsBinUnderflow(thisParticleTimeBin)){
					cerr<<"G4MuonCounterReconstructor::ProcessStationHits(): Cannot assign this time ("<<TimeOfMaxEdep<<" ns) to any time bin...exit!"<<endl;
					exit(1);
				}
				else{
					//ParticleSimDataCollectionInTimeBin[thisParticleTimeBin-1].push_back(currentParticleSimDataCollection.at(j));
					aNewParticleSimData->CreateTrackPoints();
					ParticleSimDataCollectionInTimeBin[thisParticleTimeBin-1].push_back(*aNewParticleSimData);
				}
			#else
				if(thisParticleTimeBin==0 || thisParticleTimeBin==dummyTimeHisto->GetNbinsX()+1){
					cerr<<"G4MuonCounterReconstructor::ProcessStationHits(): Cannot assign this time ("<<TimeOfMaxEdep<<" ns) to any time bin...exit!"<<endl;
					exit(1);
				}
				else{
					//ParticleSimDataCollectionInTimeBin[thisParticleTimeBin-1].push_back(currentParticleSimDataCollection.at(j));
					aNewParticleSimData->CreateTrackPoints();
					ParticleSimDataCollectionInTimeBin[thisParticleTimeBin-1].push_back(*aNewParticleSimData);
				}
			#endif
		
		}//close if empty vect


	}//end loop sim particles
	cout<<endl;
	
	dummyTimeHisto->Delete();

	cout<<"****      ProcessStationHits()   ******"<<endl;
	cout<<nHitsForThisStation<<" hits for this station"<<endl;
	cout<<"*** HITS IN TIME BIN ***"<<endl;
	for(int k=0;k<fNSamplingBins-1;k++){
		cout<< ScintHitInTimeBin[k].size() <<" --- ";
	}
	cout<< ScintHitInTimeBin[fNSamplingBins-1].size() <<endl;
	cout<< endl;
	cout<<"*** PARTICLES IN TIME BIN ***"<<endl;
	for(int k=0;k<fNSamplingBins-1;k++){
		cout<< ParticleSimDataCollectionInTimeBin[k].size() <<" --- ";
	}
	cout<<"**************************************"<<endl;


	//## Loop over time bins 
	//## Do the tracking for each time bin
	for(int k=0;k<fNSamplingBins;k++){
		if(ScintHitInTimeBin[k].size()>0)
			cout<<"*** TIME BIN NO. "<<k+1<<" ***"<<endl;
		currentStationSimData->fScintHitCollectionInTimeBin= ScintHitInTimeBin[k];	
		currentStationSimData->fParticleSimDataCollectionInTimeBin= ParticleSimDataCollectionInTimeBin[k];
		currentStationSimData->CreateTrackPoints();

		FillSimStationInfo();
	}//end loop time bins

}//close ProcessStationHits()


void G4MuonCounterReconstructor::FillSimStationInfo(){
	
	cout<<"G4MuonCounterReconstructor::FillSimStationInfo()"<<endl;

	StationId= currentStationSimData->fId;
	RStat= currentStationSimData->fRadius;
	RStatGrd= currentStationSimData->fRadiusGrd;

	//std::vector<TParticleSimData> currentParticleSimDataCollection= currentStationSimData->fParticleSimDataCollection;
	std::vector<TParticleSimData> currentParticleSimDataCollection= currentStationSimData->fParticleSimDataCollectionInTimeBin;


	EventTag= currentStationSimData->fStationSimDataTag;
	nHitPlanes= (int)((currentStationSimData->fStripMultiplicityX).size());
	
	//cout<<"EventTag= "<<EventTag<<endl;
	
	std::vector<int> currentStripMultiplicity= currentStationSimData->fStripMultiplicity;
	std::vector<int> currentStripMultiplicityX= currentStationSimData->fStripMultiplicityX;
	std::vector<int> currentStripMultiplicityY= currentStationSimData->fStripMultiplicityY;
	std::vector<double> currentStripAverageXHit= currentStationSimData->fAverageXHit;
	std::vector<double> currentStripAverageYHit= currentStationSimData->fAverageYHit;
	std::vector<double> currentStripWeightRX= currentStationSimData->fRX;		
	std::vector<double> currentStripWeightRY= currentStationSimData->fRY;


	//## Fill average position, multiplicity, ...
	for(int i=0;i<nHitPlanes;i++){
		AverageXHit[i]= (currentStationSimData->fAverageXHit).at(i);
		AverageYHit[i]= (currentStationSimData->fAverageYHit).at(i);
		MultiplicityX[i]= (currentStationSimData->fStripMultiplicityX).at(i);
		MultiplicityY[i]= (currentStationSimData->fStripMultiplicityY).at(i);	
		RX[i]= (currentStationSimData->fRX).at(i);
		RY[i]= (currentStationSimData->fRY).at(i);
		//cout<<"plane no. "<<i<<"  MultiplicityX,Y="<<MultiplicityX[i]<<"  "<<MultiplicityY[i]<<endl;
	}



	//## Fill multiplicity histo
	for(unsigned t=0;t<currentStripMultiplicity.size();t++)
		fStripMultiplicityInPlaneHisto[t]->Fill(currentStripMultiplicity[t]);
			
	//## Fill RecTree
	RecTree->Fill();

	std::vector<double> currentMuonNumberPerPlane= currentStationSimData->fMuonNumber;
	std::vector<double> currentEmNumberPerPlane= currentStationSimData->fEmNumber; 
	std::vector<double> currentMuonNumberAlbedoPerPlane= currentStationSimData->fMuonAlbedoNumber;
	std::vector<double> currentEmNumberAlbedoPerPlane= currentStationSimData->fEmAlbedoNumber;
		
	for(unsigned int l=0;l<currentMuonNumberPerPlane.size();l++){
		double muDensity= currentMuonNumberPerPlane[l]/fDetectorAreaSize;
		double muDensityAlbedo= currentMuonNumberAlbedoPerPlane[l]/fDetectorAreaSize;
		double muDensityErr= sqrt(currentMuonNumberPerPlane[l])/fDetectorAreaSize;
		double muDensityAlbedoErr= sqrt(currentMuonNumberAlbedoPerPlane[l])/fDetectorAreaSize;
		double emDensity= currentEmNumberPerPlane[l]/fDetectorAreaSize;
		double emDensityAlbedo= currentEmNumberAlbedoPerPlane[l]/fDetectorAreaSize;
		double emDensityErr= sqrt(currentEmNumberPerPlane[l])/fDetectorAreaSize;
		double emDensityAlbedoErr= sqrt(currentEmNumberAlbedoPerPlane[l])/fDetectorAreaSize;
		if(muDensity!=0) {
			fMuonDensityInPlaneGraph[l]->SetPoint(nMuonDensityPointsCounter[l],RStatGrd,muDensity);
			fMuonDensityInPlaneGraph[l]->SetPointError(nMuonDensityPointsCounter[l],0.,muDensityErr);
			nMuonDensityPointsCounter[l]++;
		}
		if(emDensity!=0) {
			fEmDensityInPlaneGraph[l]->SetPoint(nEmDensityPointsCounter[l],RStatGrd,emDensity);
			fEmDensityInPlaneGraph[l]->SetPointError(nEmDensityPointsCounter[l],0.,emDensityErr);
			nEmDensityPointsCounter[l]++;
		}
		if(muDensityAlbedo!=0) {
			fMuonDensityAlbedoInPlaneGraph[l]->SetPoint(nMuonDensityAlbedoPointsCounter[l],RStatGrd,muDensityAlbedo);
			fMuonDensityAlbedoInPlaneGraph[l]->SetPointError(nMuonDensityAlbedoPointsCounter[l],0.,muDensityAlbedoErr);
			nMuonDensityAlbedoPointsCounter[l]++;
		}
		if(emDensityAlbedo!=0) {
			fEmDensityAlbedoInPlaneGraph[l]->SetPoint(nEmDensityAlbedoPointsCounter[l],RStatGrd,emDensityAlbedo);
			fEmDensityAlbedoInPlaneGraph[l]->SetPointError(nEmDensityAlbedoPointsCounter[l],0.,emDensityAlbedoErr);
			nEmDensityAlbedoPointsCounter[l]++;
		}
	}//end loop particles per plane	



	//## Track reconstructor (make sense for muon tracks & at least 3-fold events)
	ThetaRec= -999;
	PhiRec= -999;

	if(!fTrackMuons){
		ostringstream infoMess;
  	infoMess << "Tracking option not selected...skip tracking!";         
  	INFO(infoMess.str().c_str());
		//return true;
		return;
	}

	if(EventTag<TStationSimData::eThreeFold){
		ostringstream infoMess;
  	infoMess << "This station has no 3-Fold events...skip tracking!";         
  	INFO(infoMess.str().c_str());
		//return true;	
		return;
	}

	
	

	//## Do particle identification?
	bool IsMuonInPID= false;
	if(fUsePID){
		// get parameters for this particle
		double m1X= MultiplicityX[0];
		double m1Y= MultiplicityY[0];
		double m2X= MultiplicityX[1];
		double m2Y= MultiplicityY[1];
		double m3X= MultiplicityX[2];
		double m3Y= MultiplicityY[2];
		double r1X= RX[0];
		double r1Y= RY[0];
		double r2X= RX[1];
		double r2Y= RY[1];
		double r3X= RX[2];
		double r3Y= RY[2];

		//## normalize inputs
		double m_min= 0.;
		double m1X_max= 33.;
		double m1Y_max= 38.;
		double m2X_max= 31.;
		double m2Y_max= 35.;
		double m3X_max= 20.;
		double m3Y_max= 18.;
		double r_min= 0.;
		double r_max= 1.;

		// normalize in range [a',b'] from range [a,b]
		// formula: x' = ((x - a) / (b - a)) * (b' - a') + a'
		double a_new= -1.;
		double b_new= 1.;
		double m1X_norm= ((m1X-m_min)/(m1X_max-m_min)) * (b_new-a_new) + a_new;
		double m1Y_norm= ((m1Y-m_min)/(m1Y_max-m_min)) * (b_new-a_new) + a_new;
		double m2X_norm= ((m2X-m_min)/(m2X_max-m_min)) * (b_new-a_new) + a_new;
		double m2Y_norm= ((m2Y-m_min)/(m2Y_max-m_min)) * (b_new-a_new) + a_new;
		double m3X_norm= ((m3X-m_min)/(m3X_max-m_min)) * (b_new-a_new) + a_new;
		double m3Y_norm= ((m3Y-m_min)/(m3Y_max-m_min)) * (b_new-a_new) + a_new;

		double r1X_norm= ((r1X-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;
		double r1Y_norm= ((r1Y-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;	
		double r2X_norm= ((r2X-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;
		double r2Y_norm= ((r2Y-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;
		double r3X_norm= ((r3X-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;
		double r3Y_norm= ((r3Y-r_min)/(r_max-r_min)) * (b_new-a_new) + a_new;

		const int nParsNN= 12;
		std::vector<double> NNInputData(nParsNN,0);
		
		NNInputData[0]= m1X_norm;
		NNInputData[1]= m1Y_norm;
		NNInputData[2]= m2X_norm;
		NNInputData[3]= m2Y_norm;
		NNInputData[4]= m3X_norm;
		NNInputData[5]= m3Y_norm;
		NNInputData[6]= r1X_norm;
		NNInputData[7]= r1Y_norm;
		NNInputData[8]= r2X_norm;
		NNInputData[9]= r2Y_norm;
		NNInputData[10]= r3X_norm;
		NNInputData[11]= r3Y_norm;

		IsMuonInPID= IsMuonEvent(NNInputData);
	}
	else{
		//## Find if there is at least one muon in the list of particles
		int nSimParticle= currentParticleSimDataCollection.size();
		
		for(unsigned int j=0;j<nSimParticle;j++){
			TParticleSimData* thisSimParticle= &(currentParticleSimDataCollection.at(j));	
			Id= thisSimParticle->fId;
			if(abs(Id)==utl::Particle::eMuon) {
				IsMuonInPID= true;	
				break;
			}			
		}//end loop sim particles			
	}//end else

	IsMuon= -1;	
	nTrueTrack= 0;
	nTrueMuonTrack= 0;

	if(!IsMuonInPID){
		ostringstream infoMess;
  	infoMess << "No muons in this event...skip tracking!";         
  	INFO(infoMess.str().c_str());
		//return true;
		return;
	}


	//################################
	//##      TRACKING
	//###############################
	cout<<"*** Muon Track Reco ***"<<endl;
	
	std::vector<TrackPoint> currentTrackPointCollection= currentStationSimData->fTrackPointCollection;
	
	//## Get true muon track info
	//## Consider as true tracks those hitting all 3 planes
	unsigned int nSimParticle= currentParticleSimDataCollection.size();
	double E, Theta, Phi, Tx, Ty, X, Y, R, RGrd;		
	double MaxRelAngle;
	TVector3 particleTrueDirection;

	for(unsigned int j=0;j<nSimParticle;j++){
		TParticleSimData* thisSimParticle= &(currentParticleSimDataCollection.at(j));	
		int Id= thisSimParticle->fId;
		int ParticleTag= thisSimParticle->fParticleSimDataTag;

		if(abs(Id)!=utl::Particle::eMuon || ParticleTag!=TParticleSimData::eThreeFold) 
			continue;
				
		E= thisSimParticle->fEnergy;
		particleTrueDirection= TVector3(-thisSimParticle->fDirection.X(),-thisSimParticle->fDirection.Y(),-thisSimParticle->fDirection.Z());
		Theta= particleTrueDirection.Theta()*180./TMath::Pi();
		Phi= particleTrueDirection.Phi()*180./TMath::Pi();
		Tx= particleTrueDirection.X()/particleTrueDirection.Z();
		Ty= particleTrueDirection.Y()/particleTrueDirection.Z();
		X= (thisSimParticle->fPosition).X();
		Y= (thisSimParticle->fPosition).Y();
		R= thisSimParticle->fRadius;
		RGrd= thisSimParticle->fRadiusGrd;
		
		//Calculate max rel angle among hits
		MaxRelAngle= CalculateMaxRelAngle(thisSimParticle);

		GenE[nTrueTrack]= E;
		GenTheta[nTrueTrack]= Theta;
		GenPhi[nTrueTrack]= Phi;
		GenX[nTrueTrack]= X;
		GenY[nTrueTrack]= Y;
		GenTx[nTrueTrack]= Tx;
		GenTy[nTrueTrack]= Ty;
		GenR[nTrueTrack]= R;
		GenMaxRelAngle[nTrueTrack]= MaxRelAngle;
		
		nTrueTrack++;
		nTrueMuonTrack++;
	}//end loop sim particles
	cout<<"nTrueTrack="<<nTrueTrack<<endl;
	
				
	std::vector<Track*> RecTrackCollection;
	std::vector<TrackPoint> ClusterTrackPointCollection;
	std::vector<TrackPoint> FinalTrackPointCollection;
		
	TrackRecoStatus= -1;
	
	if(fTrackingAlgo == eLeastSquare){ 
		TrackFinder theTrackFinder;
		theTrackFinder.TrackWithMSCorrection(fTrackWithMSCorrection);
		theTrackFinder.UseEnergyLossCorrectionInMS(fUseEnergyLossCorrectionInMS);
		theTrackFinder.UseAverageEnergyInMS(fUseAverageEnergyInMS);
		theTrackFinder.SetAverageEnergyInMS(fAverageEnergyInMS/utl::MeV);
		theTrackFinder.UseTrueEnergyInMS(fUseTrueEnergyInMS);
		theTrackFinder.IncludeMomentumInTracking(fIncludeMomentumInTracking);
			
		theTrackFinder.UseClusterInTracking(fUseClusterInTracking);
		theTrackFinder.SplitClustersInTracking(fSplitClustersInTracking);
		theTrackFinder.SetSplitClusterThreshold(fSplitClusterThreshold);
		theTrackFinder.RemoveFakeHitsInTracking(fRemoveFakeHitsInTracking);
			
		theTrackFinder.SetTrackPoints(currentTrackPointCollection);
		theTrackFinder.SetHitSeparationInCluster(fHitSeparationInCluster);
		theTrackFinder.SetHitSigmaX(fHitSigmaX);
		theTrackFinder.SetHitSigmaY(fHitSigmaY);
		theTrackFinder.SetHitSigmaZ(fHitSigmaZ);
		theTrackFinder.SetSuperPlaneDepth(fSuperPlaneDepth);
		theTrackFinder.SetSuperPlaneSizeZ(fSuperPlaneSizeZ);
		theTrackFinder.SetVerbosity(1);
		theTrackFinder.SetGenEnergy(E);
		theTrackFinder.SetGenTheta(Theta);
		theTrackFinder.SetGenPhi(Phi);
		theTrackFinder.SetGenTx(Tx);
		theTrackFinder.SetGenTx(Ty);
		theTrackFinder.SetGenVertexX(X*100.);
		theTrackFinder.SetGenVertexY(Y*100.);
		theTrackFinder.SetRunNumber(nMuonTracks);
		TrackRecoStatus= theTrackFinder.TrackReconstructor();

		RecTrackCollection= theTrackFinder.GetRecTrackCollection();
		ClusterTrackPointCollection= theTrackFinder.GetClusterTrackPoints();
		FinalTrackPointCollection= theTrackFinder.GetTrackPoints();
	}
	else if(fTrackingAlgo == eKalmanFilter){
		KFTrackFinder theKFTrackFinder;
		theKFTrackFinder.TrackWithMSCorrection(fTrackWithMSCorrection);
		theKFTrackFinder.UseEnergyLossCorrectionInMS(fUseEnergyLossCorrectionInMS);
		theKFTrackFinder.UseAverageEnergyInMS(fUseAverageEnergyInMS);
		theKFTrackFinder.SetAverageEnergyInMS(fAverageEnergyInMS/utl::MeV);
		theKFTrackFinder.UseTrueEnergyInMS(fUseTrueEnergyInMS);
		theKFTrackFinder.IncludeMomentumInTracking(fIncludeMomentumInTracking);

		theKFTrackFinder.UseClusterInTracking(fUseClusterInTracking);
		theKFTrackFinder.SplitClustersInTracking(fSplitClustersInTracking);
		theKFTrackFinder.SetSplitClusterThreshold(fSplitClusterThreshold);
		theKFTrackFinder.RemoveFakeHitsInTracking(fRemoveFakeHitsInTracking);
			
		theKFTrackFinder.SetTrackPoints(currentTrackPointCollection);
		theKFTrackFinder.SetHitSeparationInCluster(fHitSeparationInCluster);
		theKFTrackFinder.SetHitSigmaX(fHitSigmaX);
		theKFTrackFinder.SetHitSigmaY(fHitSigmaY);
		theKFTrackFinder.SetHitSigmaZ(fHitSigmaZ);
		theKFTrackFinder.SetSuperPlaneDepth(fSuperPlaneDepth);
		theKFTrackFinder.SetSuperPlaneSizeZ(fSuperPlaneSizeZ);
		theKFTrackFinder.SetGenEnergy(E);
		theKFTrackFinder.SetGenTheta(Theta);
		theKFTrackFinder.SetGenPhi(Phi);
		theKFTrackFinder.SetGenTx(Tx);
		theKFTrackFinder.SetGenTx(Ty);
		theKFTrackFinder.SetGenVertexX(X*100.);
		theKFTrackFinder.SetGenVertexY(Y*100.);
		theKFTrackFinder.SetVerbosity(1);
		theKFTrackFinder.SetRunNumber(nMuonTracks);
		TrackRecoStatus= theKFTrackFinder.TrackReconstructor();

		RecTrackCollection= theKFTrackFinder.GetRecTrackCollection();
		ClusterTrackPointCollection= theKFTrackFinder.GetClusterTrackPoints();
		FinalTrackPointCollection= theKFTrackFinder.GetTrackPoints();
	}//close else if KalmanFilter
	else{
		cerr<<"G4MuonCounterReconstructor::FillSimParticleInfo(): Invalid tracking algorithm specified...exit"<<endl;
		exit(1);
	}


		
		
	//## Get rec muon track info and cluster points
	nTrack= (int)(RecTrackCollection.size());
	nClusterHit= (int)(ClusterTrackPointCollection.size());
	nTrackHit= 0;
	
	cout<<"G4MuonCounterReconstructor::FillSimParticleInfo(): nTrack true/rec="<<nTrueTrack<<"/"<<nTrack<<endl;
	
	for(int s=0;s<nTrack;s++){
		Track* currentRecTrack= RecTrackCollection[s];

		TVector3 RecDirection= currentRecTrack->GetDirection();
		RecOpeningAngle[s]= RecDirection.Angle(particleTrueDirection)*180./TMath::Pi();

		RecTheta[s]= currentRecTrack->GetTheta();
		RecPhi[s]= currentRecTrack->GetPhi();
		RecX[s]= (currentRecTrack->GetVertexPos()).X();
		RecY[s]= (currentRecTrack->GetVertexPos()).Y();
		RecTx[s]= currentRecTrack->GetTx();
		RecTy[s]= currentRecTrack->GetTy();
		RecThetaErr[s]= currentRecTrack->GetThetaErr();
		RecPhiErr[s]= currentRecTrack->GetPhiErr();
		RecXErr[s]= (currentRecTrack->GetVertexPosErr()).X();
		RecYErr[s]= (currentRecTrack->GetVertexPosErr()).Y();
		RecTxErr[s]= currentRecTrack->GetTxErr();
		RecTyErr[s]= currentRecTrack->GetTyErr();
		TrackFitStatus[s]= currentRecTrack->GetFitStatus();	
		TrackFitChi2[s]= currentRecTrack->GetFitChi2();

		RecThetaStart[s]= currentRecTrack->GetThetaStart();
		RecPhiStart[s]= currentRecTrack->GetPhiStart();
		RecXStart[s]= (currentRecTrack->GetVertexPosStart()).X();
		RecYStart[s]= (currentRecTrack->GetVertexPosStart()).Y();
		RecTxStart[s]= currentRecTrack->GetTxStart();
		RecTyStart[s]= currentRecTrack->GetTyStart();
		
		RecThetaStartErr[s]= currentRecTrack->GetThetaStartErr();
		RecPhiStartErr[s]= currentRecTrack->GetPhiStartErr();
		RecXStartErr[s]= (currentRecTrack->GetVertexPosStartErr()).X();
		RecYStartErr[s]= (currentRecTrack->GetVertexPosStartErr()).Y();
		RecTxStartErr[s]= currentRecTrack->GetTxStartErr();
		RecTyStartErr[s]= currentRecTrack->GetTyStartErr();
	
		double x,y,z,t;
		TrackPoint aTrackPoint;
		int id= currentRecTrack->GetId();
			
		for(int p=0;p<currentRecTrack->GetNpoints();p++){
			currentRecTrack->GetPoint(p,x,y,z,t);
			currentRecTrack->GetTrackPoint(p,aTrackPoint);
			TrackHitX[nTrackHit]= x;
			TrackHitY[nTrackHit]= y;
			TrackHitZ[nTrackHit]= z;
			TrackId[nTrackHit]= id;
			TrackHitEdepX[nTrackHit]= aTrackPoint.fEdepX;
			TrackHitEdepY[nTrackHit]= aTrackPoint.fEdepY;	
			TrackHitTimeX[nTrackHit]= aTrackPoint.fTimeX;
			TrackHitTimeY[nTrackHit]= aTrackPoint.fTimeY; 
				
			currentRecTrack->GetTrackHitChi2(p,TrackHitChi2[nTrackHit]);	
			currentRecTrack->GetTrackHitExpChi2(p,TrackHitExpChi2[nTrackHit]);
			nTrackHit++;
		}
	}//end loop rec tracks	

	cout<<"nClusterHit="<<nClusterHit<<endl;
	for(int s=0;s<nClusterHit;s++){
		TrackPoint currentClusterPoint= ClusterTrackPointCollection[s];	
		ClusterHitX[s]= (currentClusterPoint.fPosition).X();
		ClusterHitY[s]= (currentClusterPoint.fPosition).Y();
		ClusterHitZ[s]= (currentClusterPoint.fPosition).Z();
	}

	if(nTrack>0){
		ThetaRec= RecTrackCollection[0]->GetTheta();
		PhiRec= RecTrackCollection[0]->GetPhi();		
	}

	//nHit= (int)(currentTrackPointCollection.size());//## ALL SELECTED HITS
	nHit= (int)(FinalTrackPointCollection.size());//## ALL HITS PASSED TO TRACKING STEP

	cout<<"nHit="<<nHit<<endl;
	for(int i=0;i<nHit;i++){
		//HitX[i]= (currentTrackPointCollection[i].fPosition).X();
		//HitY[i]= (currentTrackPointCollection[i].fPosition).Y();
		//HitZ[i]= (currentTrackPointCollection[i].fPosition).Z();
		//IsMuonHit[i]= currentTrackPointCollection[i].fIsMuon;
		//HitPlaneId[i]= currentTrackPointCollection[i].fDetectorPlaneId;
		//HitEdepX[i]= currentTrackPointCollection[i].fEdepX;
		//HitEdepY[i]= currentTrackPointCollection[i].fEdepY;	
		//HitTimeX[i]= currentTrackPointCollection[i].fTimeX;
		//HitTimeY[i]= currentTrackPointCollection[i].fTimeY;

		HitX[i]= (FinalTrackPointCollection[i].fPosition).X();
		HitY[i]= (FinalTrackPointCollection[i].fPosition).Y();
		HitZ[i]= (FinalTrackPointCollection[i].fPosition).Z();
		IsMuonHit[i]= FinalTrackPointCollection[i].fIsMuon;
		HitPlaneId[i]= FinalTrackPointCollection[i].fDetectorPlaneId;
		HitEdepX[i]= FinalTrackPointCollection[i].fEdepX;
		HitEdepY[i]= FinalTrackPointCollection[i].fEdepY;	
		HitTimeX[i]= FinalTrackPointCollection[i].fTimeX;
		HitTimeY[i]= FinalTrackPointCollection[i].fTimeY;  
		//cout<<"HitX[i]="<<HitX[i]<<"  HitY[i]="<<HitY[i]<<"  HitZ[i]="<<HitZ[i]<<"  IsMuonHit[i]="<<IsMuonHit[i]<<"  HitPlaneId[i]="<<HitPlaneId[i]<<endl;
		//cout<<"HitEdepX[i]="<<HitEdepX[i]<<"  HitEdepY[i]="<<HitEdepY[i]<<"  HitTimeX[i]="<<HitTimeX[i]<<"  HitTimeY[i]="<<HitTimeY[i]<<endl;
	}
		
		
	//## Fill track tree
	TrackTree->Fill();
		
	nMuonTracks++;

	
	//return true;

}//close G4MuonCounterReconstructor::FillSimStationInfo()


bool G4MuonCounterReconstructor::IsMuonEvent(std::vector<double> NNInputPars){

	bool isMuon= false;	

	unsigned int nPars= NNInputPars.size();
	int nParsNN= 12;
	if(nPars!=nParsNN){
		cerr<<"G4MuonCounterReconstructor::IsMuonEvent(): ERROR: Passed array with size "<< nPars <<" (NN inputs must have size "<<nParsNN<<")...exit"<<endl;
		exit(1);
	}


	double m1x= NNInputPars[0];
	double m1y= NNInputPars[1];
	double m2x= NNInputPars[2];
	double m2y= NNInputPars[3];
	double m3x= NNInputPars[4];
	double m3y= NNInputPars[5];  
	double r1x= NNInputPars[6];
	double r1y= NNInputPars[7];
	double r2x= NNInputPars[8];
	double r2y= NNInputPars[9];
	double r3x= NNInputPars[10];
	double r3y= NNInputPars[11];

	double NNOutput;

	TString command="./GetNNResponse ";
	command+= Form("%f %f %f %f %f %f %f %f %f %f %f %f",m1x,m1y,m2x,m2y,m3x,m3y,r1x,r1y,r2x,r2y,r3x,r3y);	
	cout<<"G4MuonCounterReconstructor::IsMuonEvent(): Executing sistem command= "<< command<<endl;

	std::string cmdOutputString;
	cmdOutputString= Utilities::ExecSystemCommand(command);	
	//cout<<"NNOutput="<<cmdOutputString.c_str()<<endl;
	
	std::stringstream cmdOutputStream(cmdOutputString,ios_base::in);
	cmdOutputStream>> NNOutput;
	cout<<"NNOutput="<<NNOutput<<endl;

	if(NNOutput>fNNCutInPID)
		isMuon= true;
	else 
		isMuon= false;

	return isMuon;

}//close G4MuonCounterReconstructor::IsMuonEvent()


double G4MuonCounterReconstructor::CalculateMaxRelAngle(TParticleSimData* aSimParticle){
	
	//## Loop over all particle hits 
	std::vector<TrackPoint> thisTrackPointCollection= aSimParticle->fTrackPointCollection;	
	unsigned int nPointsForThisParticle= thisTrackPointCollection.size();
	
	std::vector<TVector3> TrackSegmentCollection;
	TrackSegmentCollection.clear();
	TrackSegmentCollection.resize(0);

	for(unsigned int s=0;s<nPointsForThisParticle;s++){
		int IsMuonPoint= thisTrackPointCollection[s].fIsMuon;
		
		if(IsMuonPoint==0)
			continue;
	
		int PlaneId= thisTrackPointCollection[s].fDetectorPlaneId;
		

		TVector3 Position= thisTrackPointCollection[s].fPosition;
		
		// Combine this point with other point from other planes
		// Store list of track segments
		for(unsigned int k=s+1;k<nPointsForThisParticle;k++){
			int IsNextMuonPoint= thisTrackPointCollection[k].fIsMuon;
			int NextPlaneId= thisTrackPointCollection[k].fDetectorPlaneId;
			

			if(IsNextMuonPoint==0 || NextPlaneId==PlaneId){
				continue;
			}

			TVector3 NextPosition= thisTrackPointCollection[k].fPosition;
			TVector3 thisTrackSegment= NextPosition-Position;
			TrackSegmentCollection.push_back(thisTrackSegment);
			
		}//end loop next points for this particle	
	}//end loop points for this particle

	
	// Calculate relative angle
	double MaxRelAngle= -1.e+60;		

	for(int s=0;s<TrackSegmentCollection.size();s++){
		for(int k=s+1;k<TrackSegmentCollection.size();k++){
			double thisRelAngle= TrackSegmentCollection[s].Angle(TrackSegmentCollection[k]) * 180./TMath::Pi();
			
			if(thisRelAngle>MaxRelAngle){
				MaxRelAngle= thisRelAngle;
			}
		}//end loop next track segments
	}//end loop track segments
	
	return MaxRelAngle;

}//close G4MuonCounterReconstructor::CalculateMaxRelAngle()


int G4MuonCounterReconstructor::GeomEventTagger(TParticleSimData* aSimParticle){

	//## Get particle vertex and direction
	double X0= (aSimParticle->fPosition).X();
 	double Y0= (aSimParticle->fPosition).Y();	

	TVector3 dir= TVector3(aSimParticle->fDirection.X(),aSimParticle->fDirection.Y(),-aSimParticle->fDirection.Z());
	double Theta= dir.Theta()* 180./TMath::Pi();
	double Phi= dir.Phi()* 180./TMath::Pi();

	//cout<<"G4MuonCounterReconstructor::GeomEventTagger(): (X0,Y0)=("<<X0<<","<<Y0<<"), (Theta,Phi)=("<<Theta<<","<<Phi<<")"<<endl;

	//## Get some detector info
	double XYPlaneDistance= fDetectorData->GetXYPlaneDistance()/CLHEP::m;
	double Xcoat= (fDetectorData->GetStripCoatingSize()).X()/CLHEP::m;
	double Ycoat= (fDetectorData->GetStripCoatingSize()).Y()/CLHEP::m;
	double Zcoat=	(fDetectorData->GetStripCoatingSize()).Z()/CLHEP::m;
	//cout<<"G4MuonCounterReconstructor::GeomEventTagger(): XYPlaneDistance="<<XYPlaneDistance<<"  (Xcoat,Ycoat,Zcoat)=("<<Xcoat<<","<<Ycoat<<","<<Zcoat<<")"<<endl;

	//## Calculate hit strips according to particle geometry
  //## HitCoord is y-coordinate for a x-plane and viceversa for a x-plane
  double HitCoord_geom_up= 0.;//hit in superior surface
  double HitCoord_geom_down= 0.;//hit in bottom surface
  double HitCoord_geom_center= 0.;//hit at center of the track

  int HitId_geom_up= -999;
  int HitId_geom_down= -999;
  int HitId_geom_center= -999;   

	std::vector<int> GeomHitStripId;
	GeomHitStripId.clear();
	GeomHitStripId.assign(2*nPlanes,-1);
	//for(unsigned int i=0;i<GeomHitStripId.size();i++) GeomHitStripId[i]= -1;
 
	std::vector<double> PlaneDepth;
	PlaneDepth.clear();
	PlaneDepth.resize(0);
	for(unsigned int i=0;i<fSuperPlaneDepth.size();i++)
		PlaneDepth.push_back(fSuperPlaneDepth[i]/100.);	

  //## Find geometrical intersections track-plane and store strip id 
	//cout<<"GeomEventTagger(): Expected hit strips from geometrical direction"<<endl;
  int n=0;
  for(unsigned int s=0;s<GeomHitStripId.size();s++){
  	
  	if(s%2==0){
    	//x plane (calculate y coord.)
      HitCoord_geom_up= (PlaneDepth[s-n]-Zcoat/2.)*tan(Theta*TMath::Pi()/180.)*sin(Phi*TMath::Pi()/180.) + Y0; 
      HitCoord_geom_down= (PlaneDepth[s-n]+Zcoat/2.)*tan(Theta*TMath::Pi()/180.)*sin(Phi*TMath::Pi()/180.) + Y0; 
			              
           
      //check if current track geom intersects more than one strip
      HitId_geom_up= (int)(ceil( (HitCoord_geom_up + Xcoat/2.)/Ycoat ) );        
      HitId_geom_down= (int)(ceil( (HitCoord_geom_down + Xcoat/2.)/Ycoat ) );        
      if(HitId_geom_up == HitId_geom_down){
      	GeomHitStripId[s]= HitId_geom_up;
      }
      else{
      	//find hottest strip according to greater track length        
      	HitCoord_geom_center= (PlaneDepth[s-n])*tan(Theta*TMath::Pi()/180.)*sin(Phi*TMath::Pi()/180.) + Y0; 
        HitId_geom_center= (int)(ceil( (HitCoord_geom_center+Xcoat/2.)/Ycoat ) );
        GeomHitStripId[s]= HitId_geom_center;        
      }
	
			//cout<<"plane "<<s<<" (X):  CoordUp="<<HitCoord_geom_up<<"  CoordDown="<<HitCoord_geom_down<<"  StripIdUp="<<HitId_geom_up<<"  StripIdDown="<<HitId_geom_down<<"  StripId="<<GeomHitStripId[s]<<"  PlaneDepth["<<s-n<<"]="<<PlaneDepth[s-n]<<endl;

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

			//cout<<"plane "<<s<<" (Y):  CoordUp="<<HitCoord_geom_up<<"  CoordDown="<<HitCoord_geom_down<<"  StripIdUp="<<HitId_geom_up<<"  StripIdDown="<<HitId_geom_down<<"  StripId="<<GeomHitStripId[s]<<"  PlaneDepth["<<s-(n+1)<<"]="<<PlaneDepth[s-(n+1)]<<endl;

      n++;
    }//close plane y

    //cout<<"Plane="<<s<<"  StripNo="<< GeomHitStripId[s] << "  coord="<<HitCoord_geom_up<<endl;     
		
	}//close plane s loop
	//cout<<endl;

	//## Loop over geometrical hit strip, check if they exists 
	//## Assign the IsInsideAcceptance flag
	
	int foldCounter= 0;
	HasGeomCrossedFirstPlane= 0;
	for(unsigned int s=0;s<GeomHitStripId.size();s+=2){
		if(GeomHitStripId[s]>=1 && GeomHitStripId[s]<=nStrips && 
			 GeomHitStripId[s+1]>=1 && GeomHitStripId[s+1]<=nStrips
		){
			foldCounter++;
			if(s==0) HasGeomCrossedFirstPlane= 1;
		}			
	}//end loop geom hit strips
	
	//int FOVTag= int(foldCounter/2.);	
	int FOVTag= foldCounter;
	//cout<<"foldCounter="<<foldCounter<<"  FOVTag="<<FOVTag<<endl;	

 	return FOVTag;

}//close G4MuonCounterReconstructor::GeomEventTagger()




VModule::ResultFlag
G4MuonCounterReconstructor::Finish()
{

	if(fSaveRecInfo && fOutputFile) {
		fOutputFile->cd();
	
		if(fSaveRecTree)
			RecTree->Write();
		if(fSaveTrackTree)	
			TrackTree->Write();
	
		if(fSaveHisto){
			for(int i=0;i<nPlanes*2;i++){
				fEnergyDepositInPlaneHisto[i]->Write();
				fElectronEnergyDepositInPlaneHisto[i]->Write();
				fMuonEnergyDepositInPlaneHisto[i]->Write();
				fGammaEnergyDepositInPlaneHisto[i]->Write();
				fProtonEnergyDepositInPlaneHisto[i]->Write();
				fStripMultiplicityInPlaneHisto[i]->Write();
				fMuonDensityInPlaneGraph[i]->Write();
				fEmDensityInPlaneGraph[i]->Write();
				fMuonDensityAlbedoInPlaneGraph[i]->Write();
				fEmDensityAlbedoInPlaneGraph[i]->Write();
			}
			fTrigger3FoldEnergyHisto->Write();
			fTrigger2FoldEnergyHisto->Write();
		}

		
		fOutputFile->Close();
	}

  return eSuccess;
}

