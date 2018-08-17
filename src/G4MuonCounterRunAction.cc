/**
* @file G4MuonCounterRunAction.cc
* @class G4MuonCounterRunAction
* @brief Handle all operations to be performed at run level
*
* Manage operations before the run starts (i.e. creation of output data structures, ...), actions to be performed for each event (access and storing of simulation information), operations after the run ends (i.e. saving of simulation outputs, ...).
* @author S. Riggi
* @date 05/04/2010
*/

#include "G4RunManager.hh"
#include "G4MuonCounterRunAction.hh"
#include "G4MuonCounterRunScorer.hh"
#include "RunMessenger.hh"
#include "G4MuonCounterRecorderBase.hh"

#include "G4MuonCounterPrimaryGenerator.hh"
#include "G4ScintillatorHit.hh"
#include "G4PMTHit.hh"
#include "UserEventInformation.hh"

#include "G4MuonCounterSimulator.hh"
#include "G4MuonCounterEventAction.hh"
#include "G4MuonCounterConstruction.hh"
#include "G4MuonCounterSteppingAction.hh"
#include "G4MuonCounterTrackingAction.hh"
#include "G4MuonCounterTrajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

#include "G4SDManager.hh"

#include "TParticleSimData.hh"
#include "TScintHit.hh"
#include "TPMTHit.hh"
#include "PMTSimulator.hh"

#include <evt/Event.h>
#include <evt/ShowerSimData.h>
#include <fwk/RunController.h>

#include <utl/CoordinateSystem.h>
#include <utl/CoordinateSystemPtr.h>
#include <utl/Point.h>
#include <utl/Particle.h>
#include <utl/ReferenceEllipsoid.h>
#include <utl/AugerCoordinateSystem.h>
#include <utl/AugerUnits.h>

#include <TBranch.h>
#include <TVector3.h>

#include <vector>

#include <iostream>

using namespace std;
using namespace fwk;
using namespace utl;
using namespace evt;
using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterRunAction::G4MuonCounterRunAction(G4MuonCounterRecorderBase* r,G4MuonCounterConstruction* det, G4MuonCounterPrimaryGenerator* gen)
  :recorder(r),Detector(det), generator(gen)
{
  //create a messenger for this class
  runMessenger = new RunMessenger(this);

	fStoreFullInfo= false;

	fG4MuonCounterSimulator =
    dynamic_cast<G4MuonCounterSimulator*>(&RunController::GetInstance().GetModule("G4MuonCounterSimulatorUSC"));

	//Prepare data member for G4MuonCounterRunScorer
  //vector represents a list of MultiFunctionalDetector names
  theSDName.push_back(G4String("PlaneModuleSD"));

}

G4MuonCounterRunAction::~G4MuonCounterRunAction()
{
	theSDName.clear();
	delete runMessenger;
}



G4Run* G4MuonCounterRunAction::GenerateRun(){

  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in RE02Run.hh/cc.
  return new G4MuonCounterRunScorer(theSDName);
}


void G4MuonCounterRunAction::BeginOfRunAction(const G4Run* aRun){

  //cout<<"G4MuonCounterRunAction::BeginOfRunAction()"<<endl;

  if(recorder)recorder->RecordBeginOfRun(aRun);

	Init();

}


void G4MuonCounterRunAction::Init(){

	//## Get detector info
	nStrips= Detector->GetNumberOfStripsInPlane();
	nPlanes= Detector->GetNumberOfPlanes();
	double StripSizeY= (Detector->GetScintSize()).y()/CLHEP::cm;
	double StripCoatHousingSizeY= (Detector->GetScintCoatingSize()).y()/CLHEP::cm;
	double StripCoatThickness= Detector->GetHousingThickness()/CLHEP::cm;
	double SuperPlaneModuleSizeZ= Detector->GetSuperPlaneModuleSizeZ()/CLHEP::cm;
	double StripModuleSizeY= (StripSizeY + 2.* StripCoatThickness);
	StripSeparation= StripModuleSizeY;

	fHitSigmaX= StripModuleSizeY/2.;
	fHitSigmaY= StripModuleSizeY/2.;
	fHitSigmaZ= SuperPlaneModuleSizeZ/2.;
	
	//cout<<"G4MuonCounterRunAction::Init(): StripSizeY="<<StripSizeY<<"  StripCoatHousingSizeY="<<StripCoatHousingSizeY<<"  StripSeparation="<<StripSeparation<<endl;

	//## Init vectors
	fMuonInNumber.clear();
	fMuonInNumber.resize(0);
	fMuonInNumber.assign(2*nPlanes,0);
	fMuonOutNumber.clear();
	fMuonOutNumber.resize(0);
	fMuonOutNumber.assign(2*nPlanes,0);
	fEmInNumber.clear();
	fEmInNumber.resize(0);
	fEmInNumber.assign(2*nPlanes,0);
	fEmOutNumber.clear();
	fEmOutNumber.resize(0);
	fEmOutNumber.assign(2*nPlanes,0);
	fHadronInNumber.clear();
	fHadronInNumber.resize(0);
	fHadronInNumber.assign(2*nPlanes,0);
	fHadronOutNumber.clear();
	fHadronOutNumber.resize(0);
	fHadronOutNumber.assign(2*nPlanes,0);

	

}//close GetDetectorInfo()


void G4MuonCounterRunAction::FillEventInfo(){

	//#### FILL CURRENT EVENT INFO IN FILE
	//cout<<"G4MuonCounterRunAction::FillEventInfo()"<<endl;

	//access to current event
	const G4Event* currentEvent= G4RunManager::GetRunManager()->GetCurrentEvent();
	UserEventInformation* eventInformation= (UserEventInformation*)currentEvent->GetUserInformation();	
  G4TrajectoryContainer* trajectoryContainer=currentEvent->GetTrajectoryContainer();
  

	if(!currentEvent){
		G4cerr<<"G4MuonCounterRunAction::FillEventInfo(): Cannot get current event!...exit!"<<endl;
		exit(1);
  }


	//Get current event info
	int EventId= currentEvent->GetEventID();
  
	//#######################################
	//###   CREATE STATION SIM DATA
	//#######################################
	//Get current event from simulator
	const evt::Event* currentAugerEvent= G4MuonCounterSimulator::GetCurrentEvent();

	//Get current TEventSimData from simulator
	currentEventSimData= G4MuonCounterSimulator::GetCurrentEventSimData();

	//Get coordinate systems
	const CoordinateSystemPtr siteCS = det::Detector::GetInstance().GetSiteCoordinateSystem();
	const ShowerSimData & simShower = currentAugerEvent->GetSimShower(); 
  const CoordinateSystemPtr showerCS = simShower.GetShowerCoordinateSystem();
	const CoordinateSystemPtr localCS = simShower.GetLocalCoordinateSystem();
	const ReferenceEllipsoid& e = ReferenceEllipsoid::GetWGS84();
  
	//Get current StationSimData
	//sevt::SEvent::StationIterator fCurrentStationIt = G4MuonCounterSimulator::GetCurrentEventStationIt();	
	const sdet::Station* fCurrentDetectorStation = G4MuonCounterSimulator::GetCurrentDetectorStation();
	utl::Point fCurrentStationPosition = fCurrentDetectorStation->GetPosition();
  CoordinateSystemPtr stationCS = AugerCoordinateSystem::Create(fCurrentStationPosition, e, siteCS);


	

	//Get corresponding TStationSimData
	TStationSimData* currentStationSimData= currentEventSimData->GetSimStation(fCurrentDetectorStation->GetId());
	if(!currentStationSimData){
		//cout<<"Current sim station not found in collection, creating new sim station"<<endl; 
		currentStationSimData= new TStationSimData;
		currentStationSimData->nStrips= nStrips;
		currentStationSimData->nPlanes= nPlanes;
		currentStationSimData->fHitSigmaX= fHitSigmaX;
		currentStationSimData->fHitSigmaY= fHitSigmaY;
		currentStationSimData->fHitSigmaZ= fHitSigmaZ; 
		
		currentStationSimData->fId= fCurrentDetectorStation->GetId();	
		currentStationSimData->fName= fCurrentDetectorStation->GetName();
		currentStationSimData->fPosition= TVector3(fCurrentDetectorStation->GetPosition().GetX(siteCS)/utl::m,
						 																   fCurrentDetectorStation->GetPosition().GetY(siteCS)/utl::m,
						 																   fCurrentDetectorStation->GetPosition().GetZ(siteCS)/utl::m);	
		currentStationSimData->fRadius= (fCurrentDetectorStation->GetPosition()).GetRho(showerCS)/utl::m;
		currentStationSimData->fPhi= (fCurrentDetectorStation->GetPosition()).GetPhi(showerCS)/utl::deg;
		currentStationSimData->fRadiusGrd= (fCurrentDetectorStation->GetPosition()).GetRho(localCS)/utl::m;
		currentStationSimData->fPhiGrd= (fCurrentDetectorStation->GetPosition()).GetPhi(localCS)/utl::deg;			
	}
	int StationId= currentStationSimData->fId;
	

	//Fill a corresponding TStationSimData
	fStationSimData= new TStationSimData;
	fStationSimData->nStrips= nStrips;
	fStationSimData->nPlanes= nPlanes;
	fStationSimData->fHitSigmaX= fHitSigmaX;
	fStationSimData->fHitSigmaY= fHitSigmaY;
	fStationSimData->fHitSigmaZ= fHitSigmaZ; 
	fStationSimData->fId= fCurrentDetectorStation->GetId();	
	fStationSimData->fName= fCurrentDetectorStation->GetName();
	fStationSimData->fPosition= TVector3(fCurrentDetectorStation->GetPosition().GetX(siteCS)/utl::m,
						 																fCurrentDetectorStation->GetPosition().GetY(siteCS)/utl::m,
						 																fCurrentDetectorStation->GetPosition().GetZ(siteCS)/utl::m);	
	fStationSimData->fRadius= (fCurrentDetectorStation->GetPosition()).GetRho(showerCS)/utl::m;
	fStationSimData->fPhi= (fCurrentDetectorStation->GetPosition()).GetPhi(showerCS)/utl::deg;
	fStationSimData->fRadiusGrd= (fCurrentDetectorStation->GetPosition()).GetRho(localCS)/utl::m;
	fStationSimData->fPhiGrd= (fCurrentDetectorStation->GetPosition()).GetPhi(localCS)/utl::deg;
	fStationSimData->fMuonNumber= currentStationSimData->fMuonNumber;
	fStationSimData->fEmNumber= currentStationSimData->fEmNumber;


	//#######################################
	//###   CREATE PARTICLE SIM DATA
	//#######################################
	//Get current particle
	sevt::StationSimData::ParticleIterator fCurrentParticleIt = G4MuonCounterSimulator::GetCurrentParticleIt();
	CoordinateSystemPtr csDir = fCurrentParticleIt->GetDirection().GetCoordinateSystem();

	//Fill a corresponding TParticleSimData
	TParticleSimData currentParticleSimData;
	currentParticleSimData.fId= fCurrentParticleIt->GetType();
	currentParticleSimData.fName= fCurrentParticleIt->GetName();
	//currentParticleSimData.fPDGCode= fCurrentParticleIt.Get;
	currentParticleSimData.fGlobalPosition= TVector3(fCurrentParticleIt->GetPosition().GetX(siteCS)/utl::m,
																						       fCurrentParticleIt->GetPosition().GetY(siteCS)/utl::m,
																						       fCurrentParticleIt->GetPosition().GetZ(siteCS)/utl::m);
	currentParticleSimData.fPosition= TVector3(fCurrentParticleIt->GetPosition().GetX(stationCS)/utl::m,
																						 fCurrentParticleIt->GetPosition().GetY(stationCS)/utl::m,
																						 fCurrentParticleIt->GetPosition().GetZ(stationCS)/utl::m);

	currentParticleSimData.fRadius= (fCurrentParticleIt->GetPosition()).GetRho(showerCS)/utl::m;
	currentParticleSimData.fPsi= (fCurrentParticleIt->GetPosition()).GetPhi(showerCS)/utl::deg;
	currentParticleSimData.fRadiusGrd= (fCurrentParticleIt->GetPosition()).GetRho(localCS)/utl::m;
	currentParticleSimData.fPsiGrd= (fCurrentParticleIt->GetPosition()).GetPhi(localCS)/utl::deg;

	currentParticleSimData.fDirection= TVector3(fCurrentParticleIt->GetDirection().GetX(csDir),
																						  fCurrentParticleIt->GetDirection().GetY(csDir),
																						  fCurrentParticleIt->GetDirection().GetZ(csDir));
	currentParticleSimData.fMomentum= TVector3(fCurrentParticleIt->GetMomentum().GetX(csDir)/utl::MeV,
																						 fCurrentParticleIt->GetMomentum().GetY(csDir)/utl::MeV,
																						 fCurrentParticleIt->GetMomentum().GetZ(csDir)/utl::MeV);
	currentParticleSimData.fEnergy= fCurrentParticleIt->GetKineticEnergy()/utl::MeV;
	currentParticleSimData.fTheta= fCurrentParticleIt->GetDirection().GetTheta(csDir)/utl::deg;
	currentParticleSimData.fPhi= fCurrentParticleIt->GetDirection().GetPhi(csDir)/utl::deg;
	currentParticleSimData.fTime= fCurrentParticleIt->GetTime().GetInterval()/utl::ns;

	currentParticleSimData.nStrips= nStrips;
	currentParticleSimData.nPlanes= nPlanes;
	currentParticleSimData.StripSeparation= StripSeparation;
	currentParticleSimData.fHitSigmaX= fHitSigmaX;
	currentParticleSimData.fHitSigmaY= fHitSigmaY;
	currentParticleSimData.fHitSigmaZ= fHitSigmaZ;

	currentParticleSimData.fDetectorData= fG4MuonCounterSimulator->GetCurrentMuonDetector();

	cout<<"G4MuonCounterRunAction(): INFO: Theta="<<fCurrentParticleIt->GetDirection().GetTheta(csDir)/utl::deg<<"  Phi="<<fCurrentParticleIt->GetDirection().GetPhi(csDir)/utl::deg<<endl;
	//cout<<"theta,phi (runAction)="<< currentParticleSimData.fTheta << "  "<<currentParticleSimData.fPhi <<endl;

	//TStationSimData& stationSimData = *(fG4MuonCounterSimulator->GetCurrentStationSimData());
	//int StationId= stationSimData.fId; 
	
	//########################################

  //Initialize hit structures
  int scintCollID= -1;
	int pmtCollID= -1;

  int StripId= -1;
	int PlaneId= -1;
	int SuperPlaneId= -1;
	double Edep= 0.;
  double TotEdep=0;
  int Nphotons= 0;
  int Nphotons_Scint= 0;
	int Nphotons_Cerenk= 0;
	int Nphotons_WLS= 0;
  int Nphotons_Abs= 0;
  int Nphotons_AbsBoundary= 0;
	
	double Time=0.;
	int TrackId=-1;
	int ParticleType=-1;
	TVector3 Position(0,0,0);
	TVector3 StripPosition(0,0,0);
	TVector3 TrackDirection(0,0,0);

	int PMTId=-1;
	int PhotonCounts=0;
	int PhotonCounts_scint=0;
  int PhotonCounts_cerenk=0;
	int PhotonCounts_wls=0;

  int pmtThreshold= 1;

  //get user event info
  std::vector<double> emissAngles= eventInformation->GetPhotonEmissionAngle();
  std::vector<double> emissTime= eventInformation->GetPhotonEmissionTime();
	std::vector<double> emissWavelength= eventInformation->GetPhotonEmissionWavelength();
	std::vector<int> processType= eventInformation->GetProcessType();
  std::vector<int> stripNo= eventInformation->GetStripNumber();
	std::vector<int> planeNo= eventInformation->GetPlaneNumber();
	std::vector<int> superplaneNo= eventInformation->GetSuperPlaneNumber();

  std::vector<double> trackLength= eventInformation->GetPhotonTrackLength();
	std::vector<double> trackHorizLength= eventInformation->GetPhotonHorizontalTrackLength();
  std::vector<int> stripNo2= eventInformation->GetStripNumber2();
	std::vector<int> planeNo2= eventInformation->GetPlaneNumber2();
	std::vector<int> superplaneNo2= eventInformation->GetSuperPlaneNumber2();
	std::vector<int> absFlag= eventInformation->GetAbsorptionFlag();


	//Get hit collections info
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(scintCollID<0) scintCollID=SDman->GetCollectionID("scintHitCollection");
  if(pmtCollID<0) pmtCollID=SDman->GetCollectionID("pmtHitCollection");


	ScintHitsCollection* SHC = 0;
  PMTHitsCollection* PHC = 0;
  G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();

  
  if(HCE){
    if(scintCollID>=0) SHC = (ScintHitsCollection*)(HCE->GetHC(scintCollID));
    if(pmtCollID>=0) PHC = (PMTHitsCollection*)(HCE->GetHC(pmtCollID));
  }



	std::vector<double> muInCounts= eventInformation->GetInMuonsAtPlaneSurface();
	std::vector<double> muOutCounts= eventInformation->GetOutMuonsAtPlaneSurface();
	std::vector<double> emInCounts= eventInformation->GetInEmAtPlaneSurface();
	std::vector<double> emOutCounts= eventInformation->GetOutEmAtPlaneSurface();
	std::vector<double> hadronInCounts= eventInformation->GetInHadronsAtPlaneSurface();
	std::vector<double> hadronOutCounts= eventInformation->GetOutHadronsAtPlaneSurface();

	fMuonInNumberInEvent.clear();
	fMuonInNumberInEvent.resize(0);
	fMuonInNumberInEvent.assign(2*nPlanes,0);
	fMuonOutNumberInEvent.clear();
	fMuonOutNumberInEvent.resize(0);
	fMuonOutNumberInEvent.assign(2*nPlanes,0);
	fEmInNumberInEvent.clear();
	fEmInNumberInEvent.resize(0);
	fEmInNumberInEvent.assign(2*nPlanes,0);
	fEmOutNumberInEvent.clear();
	fEmOutNumberInEvent.resize(0);
	fEmOutNumberInEvent.assign(2*nPlanes,0);
	fHadronInNumberInEvent.clear();
	fHadronInNumberInEvent.resize(0);
	fHadronInNumberInEvent.assign(2*nPlanes,0);
	fHadronOutNumberInEvent.clear();
	fHadronOutNumberInEvent.resize(0);
	fHadronOutNumberInEvent.assign(2*nPlanes,0);

	cout << "=============================================================" <<endl;
  cout << " Number of event processed for event : "<< currentEvent->GetEventID() << endl;
  cout << "=============================================================" <<endl;
	for(unsigned int i=0;i<2*nPlanes;i++){ 
		//cout << "Plane no. " << i << " : muons (in/out) "<< muInCounts[i] << "/" << muOutCounts[i] << "  em (in/out) "<< emInCounts[i] << "/" << emOutCounts[i] << "  hadrons (in/out) "<< hadronInCounts[i] << "/" << hadronOutCounts[i]<<endl;
		fMuonInNumber[i]+= muInCounts[i];
		fEmInNumber[i]+= emInCounts[i];
		fHadronInNumber[i]+= hadronInCounts[i];
		fMuonOutNumber[i]+= muOutCounts[i];
		fEmOutNumber[i]+= emOutCounts[i];
		fHadronOutNumber[i]+= hadronOutCounts[i];
	
		fMuonInNumberInEvent[i]= muInCounts[i];
		fEmInNumberInEvent[i]= emInCounts[i];
		fHadronInNumberInEvent[i]= hadronInCounts[i];
		fMuonOutNumberInEvent[i]= muOutCounts[i];
		fEmOutNumberInEvent[i]= emOutCounts[i];
		fHadronOutNumberInEvent[i]= hadronOutCounts[i];
	}//end loop planes
  //cout << "============================================="<<endl;


	//####### STORE SCINT HIT ############
	//first clean and resize hit collection for current event
	//add other vectors
	for(unsigned k=0;k<fScintHitCollection.size();k++){
		(fScintHitCollection[k].PhotonTrackLength).clear();
		(fScintHitCollection[k].PhotonTrackLength).resize(0);

		(fScintHitCollection[k].PhotonTrackDistance).clear();
		(fScintHitCollection[k].PhotonTrackDistance).resize(0);

		(fScintHitCollection[k].PhotonEmissionAngle).clear();
		(fScintHitCollection[k].PhotonEmissionAngle).resize(0);

  	(fScintHitCollection[k].PhotonEmissionTime).clear();
		(fScintHitCollection[k].PhotonEmissionTime).resize(0);

		(fScintHitCollection[k].PhotonEmissionWavelength).clear();
		(fScintHitCollection[k].PhotonEmissionWavelength).resize(0);

  	(fScintHitCollection[k].PhotonProcessType).clear();
		(fScintHitCollection[k].PhotonProcessType).resize(0);

		fScintHitCollection[k].Nphotons= 0;
    fScintHitCollection[k].Nphotons_scint= 0;
    fScintHitCollection[k].Nphotons_cerenk= 0;
		fScintHitCollection[k].Nphotons_wls= 0;
    fScintHitCollection[k].Nphotons_abs= 0;
		fScintHitCollection[k].Nphotons_absBoundary= 0;
	}//close for	

	fScintHitCollection.clear();
	fScintHitCollection.resize(0);

	
  if(SHC){
    int Nhits_scint = SHC->entries();
		//cout<<"Nhits in scintillator="<<Nhits_scint<<endl;

    G4ThreeVector eWeightPos(0.);     
    double EdepMax=0;

 		G4ThreeVector HitPos(0.);
		G4ThreeVector HitStripPos(0.);
		G4ThreeVector HitTrackDir(0.);
	
		int collection_index= -1;
		bool hitInCollection= false;
		int selcollection_index= -1;
		bool hitInSelectedCollection= false;
		bool isSelectedHit= false;

    for(int i=0;i<Nhits_scint;i++){ //gather info on hits in scintillator

			TScintHit currentHit;
			TScintHit currentSelectedHit;

			//init hit
			(currentHit.PhotonTrackLength).clear();
			(currentHit.PhotonTrackLength).resize(0);

			(currentHit.PhotonTrackDistance).clear();
			(currentHit.PhotonTrackDistance).resize(0);

			(currentHit.PhotonEmissionAngle).clear();
			(currentHit.PhotonEmissionAngle).resize(0);

  		(currentHit.PhotonEmissionTime).clear();
			(currentHit.PhotonEmissionTime).resize(0);

			(currentHit.PhotonEmissionWavelength).clear();
			(currentHit.PhotonEmissionWavelength).resize(0);

  		(currentHit.PhotonProcessType).clear();
			(currentHit.PhotonProcessType).resize(0);

			(currentHit.PhotocathodeSurfacePosition).clear();
			(currentHit.PhotocathodeSurfacePosition).resize(0);

			(currentHit.Edep).clear();
			(currentHit.Edep).resize(0);

			(currentHit.Time).clear();
			(currentHit.Time).resize(0);
	
			(currentHit.TrackId).clear();
			(currentHit.TrackId).resize(0);

			(currentHit.ParticleType).clear();
			(currentHit.ParticleType).resize(0);

			(currentHit.Position).clear();
			(currentHit.Position).resize(0);

			(currentHit.TrackDirection).clear();
			(currentHit.TrackDirection).resize(0);

			std::vector<TVector3> PhotocathodeSurfacePosTVector3;
			PhotocathodeSurfacePosTVector3.clear();
			PhotocathodeSurfacePosTVector3.resize(0);

			currentHit.StripId= -1;
			currentHit.PlaneId= -1;
			currentHit.SuperPlaneId= -1;
			currentHit.Etot= 0;
			(currentHit.StripPosition).SetXYZ(-1,-1,-1);
			currentHit.Nphotons= 0;
      currentHit.Nphotons_scint= 0;
      currentHit.Nphotons_cerenk= 0;
			currentHit.Nphotons_wls= 0;
      currentHit.Nphotons_abs= 0;
		  currentHit.Nphotons_absBoundary= 0;

			//get info for current hit
      Edep=(*SHC)[i]->GetEdep();
      StripId= (*SHC)[i]->GetSciCopyNo();	
			PlaneId= (*SHC)[i]->GetPlaneCopyNo();		
			SuperPlaneId= (*SHC)[i]->GetSuperPlaneCopyNo();	
			Time= (*SHC)[i]->GetTimeHit();
			TrackId= (*SHC)[i]->GetTrackId();
			ParticleType= (*SHC)[i]->GetTrackPartType();
  		HitPos= (*SHC)[i]->GetPos();
			Position= TVector3(HitPos.x()/CLHEP::cm,HitPos.y()/CLHEP::cm,HitPos.z()/CLHEP::cm);
			HitStripPos= (*SHC)[i]->GetStripPos();
			StripPosition= TVector3(HitStripPos.x()/CLHEP::cm,HitStripPos.y()/CLHEP::cm,HitStripPos.z()/CLHEP::cm);
			HitTrackDir= (*SHC)[i]->GetTrackDirection();
			TrackDirection= TVector3(HitTrackDir.x(),HitTrackDir.y(),HitTrackDir.z());

			std::vector<G4ThreeVector> PhotocathodeSurfacePos= (*SHC)[i]->GetPhotocathodeSurfacePos();
			for(unsigned int l=0;l<PhotocathodeSurfacePos.size();l++){
				PhotocathodeSurfacePosTVector3.push_back(TVector3(PhotocathodeSurfacePos[l].x()/CLHEP::cm,PhotocathodeSurfacePos[l].y()/CLHEP::cm,PhotocathodeSurfacePos[l].z()/CLHEP::cm));
			}

			cout<<"G4MuonCounterRunAction(): INFO: Hit no. "<<i+1<<"  (S,P,SP)=("<<StripId<<","<< PlaneId<<","<<SuperPlaneId<<")  Edep="<<Edep<<endl;

			/*
			cout<<"Get expPMTTime vector"<<endl;
			std::vector<double> ExpPMTHitTime= (*SHC)[i]->GetExpPMTTime();
			unsigned int nPMTForThisStrip= ExpPMTHitTime.size();
			cout<<"nPMTForThisStrip="<<nPMTForThisStrip<<endl;
			for(unsigned int g=0;g<nPMTForThisStrip;g++){
      	cout<<"TrackId="<<TrackId<<"  (SP,P,S)= ("<<SuperPlaneId<<","<<PlaneId<<","<<StripId<<")  Edep="<<Edep/CLHEP::MeV<<"  time="<<Time/CLHEP::ns<< "  expTime "<<g+1<<" ="<<ExpPMTHitTime[g]<<endl;
			}
			*/

      eventInformation->IncEDep(Edep); //sum up the edep
      eWeightPos += (*SHC)[i]->GetPos()*Edep;//calculate energy weighted pos
      if(Edep>EdepMax){
				EdepMax=Edep;//store max energy deposit
				G4ThreeVector posMax=(*SHC)[i]->GetPos();
				eventInformation->SetPosMax(posMax,Edep);
      }
    
			
			//search if the current strip has already been stored in collection
			hitInCollection= false;

			for(unsigned k=0;k<fScintHitCollection.size();k++){
				if((*SHC)[i]->GetSciCopyNo()== fScintHitCollection[k].StripId && (*SHC)[i]->GetPlaneCopyNo()== fScintHitCollection[k].PlaneId && (*SHC)[i]->GetSuperPlaneCopyNo()== fScintHitCollection[k].SuperPlaneId){
					//hit already stored in vector for this strip
					hitInCollection= true;
					collection_index= k;
					break;
				}
			}//close for

			if(hitInCollection){
				//add info to strip vector element
				(fScintHitCollection[collection_index].Edep).push_back(Edep);	
				(fScintHitCollection[collection_index].Time).push_back(Time);
				(fScintHitCollection[collection_index].TrackId).push_back(TrackId);
				(fScintHitCollection[collection_index].ParticleType).push_back(ParticleType);
				(fScintHitCollection[collection_index].Position).push_back(Position);
				(fScintHitCollection[collection_index].TrackDirection).push_back(TrackDirection);
				//cout<<"dir=("<<TrackDirection.X()<<","<<TrackDirection.Y()<<","<<TrackDirection.Z()<<")"<<endl;

				(fScintHitCollection[collection_index].Etot)+= Edep;

				/*
				cout<<"Appending to existing exp time entries"<<endl;		
				for(unsigned int hh=0;hh<nPMTForThisStrip;hh++) {
					(fScintHitCollection[collection_index].ExpPMTTime[hh]).push_back(ExpPMTHitTime[hh]);
				}
				cout<<"done"<<endl;
				*/
			}
			else{
				//new hit vector
				//set tree branch values
				(currentHit.Edep).push_back(Edep);	
				(currentHit.Time).push_back(Time);
				(currentHit.TrackId).push_back(TrackId);
				(currentHit.ParticleType).push_back(ParticleType);
				(currentHit.Position).push_back(Position);
				(currentHit.TrackDirection).push_back(TrackDirection);
				//cout<<"dir=("<<TrackDirection.X()<<","<<TrackDirection.Y()<<","<<TrackDirection.Z()<<")"<<endl;
				

				currentHit.StripId= StripId;
				currentHit.PlaneId= PlaneId;	
				currentHit.SuperPlaneId= SuperPlaneId;	
				currentHit.Etot= Edep;
				currentHit.StripPosition= StripPosition;
				currentHit.PhotocathodeSurfacePosition= PhotocathodeSurfacePosTVector3;

				/*
				cout<<"Adding new exp time entry"<<endl;
				std::vector< std::vector<double> > ExpHitTimeAtPMTAnode;
				for(int hh=0;hh<nPMTForThisStrip;hh++) {
					ExpHitTimeAtPMTAnode.push_back ( std::vector<double>() );
					ExpHitTimeAtPMTAnode[hh].push_back(ExpPMTHitTime[hh]);
				}
				currentHit.ExpPMTTime= ExpHitTimeAtPMTAnode;
				*/
				
				fScintHitCollection.push_back(currentHit);//add this strip in the collection
				
			}//close else
			
    }//close for hits

		
    if(eventInformation->GetEDep()==0.) {	
			//cout<<"No hits in scintillators for this event."<<endl;   	
		}
    else{
      //Finish calculation of energy weighted position
      eWeightPos/=eventInformation->GetEDep();
      eventInformation->SetEWeightPos(eWeightPos);
      //cout << "\t Energy weighted position of hits in Scintillator : "
	    //       << eWeightPos/CLHEP::mm << endl;
    }//close else
    
    //cout << "\t Total energy deposition in scintillator : "
	  // 			 << eventInformation->GetEDep()/CLHEP::MeV << " (MeV)" << endl;
  

		for(unsigned k=0;k<fScintHitCollection.size();k++){
			int StripId_collection= fScintHitCollection[k].StripId;
			int PlaneId_collection= fScintHitCollection[k].PlaneId;
			int SuperPlaneId_collection= fScintHitCollection[k].SuperPlaneId;

			for(unsigned int j=0;j<trackLength.size();j++){
				if(stripNo2[j]==StripId_collection && planeNo2[j]==PlaneId_collection && superplaneNo2[j]==SuperPlaneId_collection){	
					if(fStoreFullInfo){
						(fScintHitCollection[k].PhotonTrackLength).push_back(trackLength[j]);
						(fScintHitCollection[k].PhotonTrackDistance).push_back(trackHorizLength[j]);
					}

					if(absFlag[j]==G4MuonCounterSteppingAction::kAbsInBoundary) fScintHitCollection[k].Nphotons_absBoundary++;
					else fScintHitCollection[k].Nphotons_abs++;
				}//close if			
			}//close for j

			for(unsigned int j=0;j<emissAngles.size();j++){
				if(stripNo[j]==StripId_collection && planeNo[j]==PlaneId_collection && superplaneNo[j]==SuperPlaneId_collection) {
			
					if(fStoreFullInfo){
						(fScintHitCollection[k].PhotonEmissionAngle).push_back(emissAngles[j]);
						(fScintHitCollection[k].PhotonEmissionTime).push_back(emissTime[j]);
						(fScintHitCollection[k].PhotonEmissionWavelength).push_back(emissWavelength[j]);
						(fScintHitCollection[k].PhotonProcessType).push_back(processType[j]);
					}	

					if(processType[j]==G4MuonCounterTrackingAction::kScintillation) fScintHitCollection[k].Nphotons_scint++;
					else if(processType[j]==G4MuonCounterTrackingAction::kCerenkov) fScintHitCollection[k].Nphotons_cerenk++;
					else if(processType[j]==G4MuonCounterTrackingAction::kWLS) fScintHitCollection[k].Nphotons_wls++;

				}//close if
			}//close for

      (fScintHitCollection[k].Nphotons)= fScintHitCollection[k].Nphotons_scint+fScintHitCollection[k].Nphotons_cerenk+fScintHitCollection[k].Nphotons_wls;
				
		}//close loop collection

		
		//##################################
		//##    Fill StationSimData
		//##################################
		unsigned int nHitStrips= fScintHitCollection.size();
			
		//cout<<endl;
		//cout<<"*** Station no "<<StationId<<" ***"<<endl;
		//cout<< nHitStrips <<" hit strips found"<<endl;	
		//loop over all hit strips 	
		for(unsigned int kk=0;kk<nHitStrips;kk++){
			TScintHit currentScintHit= fScintHitCollection[kk];
			TScintHit currentSelectedScintHit;
			int StripId= currentScintHit.StripId;
			int PlaneId= currentScintHit.PlaneId;
			int SuperPlaneId= currentScintHit.SuperPlaneId;
			int AbsPlaneId= 2*SuperPlaneId+PlaneId;
			double TotEnergyDep= currentScintHit.Etot;

			
			bool isGoodHit= currentStationSimData->SelectScintHit(currentScintHit,currentSelectedScintHit);
			currentStationSimData->AddScintHit(currentScintHit,TStationSimData::eInFullCollection);
			currentParticleSimData.AddScintHit(currentScintHit);

			if(isGoodHit) {
				currentStationSimData->AddScintHit(currentSelectedScintHit,TStationSimData::eInSelectionCollection);
				currentParticleSimData.AddSelScintHit(currentSelectedScintHit);
			}

			cout<<"G4MuonCounterRunAction(): INFO:  (SP,P,S)= ("<<SuperPlaneId<<","<<PlaneId<<","<<StripId<<") AP="<<AbsPlaneId<<" ==> Edep="<<TotEnergyDep/CLHEP::MeV<<" IsGoodHit?"<<isGoodHit<<endl;
			
			
			unsigned int nHitsForThisStrip= (currentScintHit.Edep).size();
			//unsigned int nPMTForThisStrip= (currentScintHit.ExpPMTTime).size();
			for(unsigned int ss=0;ss<nHitsForThisStrip;ss++){	
				double Edep= currentScintHit.Edep[ss];
				double Time= currentScintHit.Time[ss];			
				int TrackId= currentScintHit.TrackId[ss];
				int ParticleType= currentScintHit.ParticleType[ss];
				TVector3 TrackDirection= currentScintHit.TrackDirection[ss];
				//cout<<"*** RunAction ***"<<endl;
				//TrackDirection.Print();
				//cout<<"******************"<<endl;
				
				//cout<<"Track "<<TrackId<<"  ParticleType="<<ParticleType<<"  (SP,P,S)= ("<<SuperPlaneId<<","<<PlaneId<<","<<StripId<<") AP="<<AbsPlaneId<<" ==> Edep="<<Edep/CLHEP::MeV<<"  t="<<Time/CLHEP::ns<<endl;
				/*
				for(unsigned int ww=0;ww<nPMTForThisStrip;ww++) {
					double expPMTTime= currentScintHit.ExpPMTTime[ww][ss];
					cout<<"Track "<<TrackId<<"  (SPlane,Plane,Strip)= ("<<SuperPlaneId<<","<<PlaneId<<","<<StripId<<") ==> Edep="<<Edep/CLHEP::MeV<<"  t="<<Time/CLHEP::ns<<"  texp "<<ww+1<<" ="<<expPMTTime/CLHEP::ns<<endl;
				}
				*/
				
			}//end loop nHits for this strip
		
		}//end loop nHits for this station	
		//cout<<endl;
		/*
		cout<<"*********************"<<endl;
		cout<<"*** Station no "<<StationId<<" ***"<<endl;
		cout<< "Hits (all/sel): "<<(currentParticleSimData.fSelScintHitCollection).size()<<"/"<<(currentParticleSimData.fSelScintHitCollection).size() <<endl;
		cout<<"*********************"<<endl;
		cout<<endl;
		*/
		currentParticleSimData.CreateTrackPoints();
			
		//## draw all scintillator hits
		SHC->DrawAllHits();

  }//close if ScintHits


	//####### STORE PMT HIT ############
	//first clean and resize hit collection for current event
	fPMTHitCollection.clear();
	fPMTHitCollection.resize(0);


	//Hits in PMTs
  if(PHC){
		int Nhits_pmt = PHC->entries();
		//cout<<"Nhits in PMTs="<<Nhits_pmt<<endl;

    G4ThreeVector reconPos(0,0,0);
		G4ThreeVector HitPos(0,0,0);
		int collection_index=-1;
		bool hitInCollection= false;

    for(int i=0;i<Nhits_pmt;i++){ //gather info on hits in scintillator

			TPMTHit currentHit;

			//clean and resize vectors for current Hit
			(currentHit.PhotonTime).clear();
			(currentHit.PhotonTime).resize(0);

			(currentHit.PhotonEnergy).clear();
			(currentHit.PhotonEnergy).resize(0);

			(currentHit.PhotonPosition).clear();
			(currentHit.PhotonPosition).resize(0);

			currentHit.StripId= -1;
			currentHit.PlaneId= -1;
			currentHit.SuperPlaneId= -1;
			currentHit.PMTId= -1;
			currentHit.PhotonCounts= 0;
      currentHit.PhotonCounts_scint= 0;
      currentHit.PhotonCounts_cerenk= 0;
			currentHit.PhotonCounts_wls= 0;
      currentHit.PhotonType= 0;
		  currentHit.Position(0);


			//get info for current hit
      StripId= (*PHC)[i]->GetStripNumber();
			PlaneId= (*PHC)[i]->GetPlaneNumber();
			SuperPlaneId= (*PHC)[i]->GetSuperPlaneNumber();
			PMTId= (*PHC)[i]->GetPMTNumber();
			PhotonCounts= (*PHC)[i]->GetPhotonCount();
			PhotonCounts_scint= (*PHC)[i]->GetPhotonCount_scint();
			PhotonCounts_cerenk= (*PHC)[i]->GetPhotonCount_cerenk();
			PhotonCounts_wls= (*PHC)[i]->GetPhotonCount_wls();
			
  		HitPos= (*PHC)[i]->GetPos();
			Position= TVector3(HitPos.x()/CLHEP::cm,HitPos.y()/CLHEP::cm,HitPos.z()/CLHEP::cm);
			
			std::vector<double> PhotonTime= (*PHC)[i]->GetPhotonTime();
			std::vector<double> PhotonEnergy= (*PHC)[i]->GetPhotonEnergy();
	    std::vector<G4ThreeVector> PhotonPosition= (*PHC)[i]->GetPhotonPosition();
			

      //add info of current hit to user event vectors 
			for(unsigned int j=0;j<PhotonTime.size();j++){
				//cout<<"PhotonTime["<<j<<"]= "<<PhotonTime[j]<<endl;
				eventInformation->InsertPMTPhotonTime(PhotonTime[j]);
        //cout<<"PhotonEnergy["<<j<<"]= "<<PhotonEnergy[j]<<endl;
				eventInformation->InsertPMTPhotonEnergy(PhotonEnergy[j]);
				eventInformation->InsertPMTPhotonPosition(PhotonPosition[j]);
				eventInformation->InsertPMTStripNumber(StripId);
				eventInformation->InsertPMTPlaneNumber(PlaneId);
				eventInformation->InsertPMTSuperPlaneNumber(SuperPlaneId);
				eventInformation->InsertPMTNumber(PMTId);
			}//close for time entries


			//fill vectors distiguishing among different PMTs
			eventInformation->InsertPMTId(PMTId);	
			eventInformation->InsertPMTStripId(StripId);	
			eventInformation->InsertPMTPlaneId(PlaneId);
			eventInformation->InsertPMTSuperPlaneId(SuperPlaneId);			
			eventInformation->InsertPMTCounts(PhotonCounts);					
			eventInformation->InsertPMTScintillationCounts(PhotonCounts_scint);
			eventInformation->InsertPMTCerenkovCounts(PhotonCounts_cerenk);
			eventInformation->InsertPMTWLSCounts(PhotonCounts_wls);


			//fill cumulative info for all PMTs
      eventInformation->IncHitCount(PhotonCounts);
			eventInformation->IncHitScintCount(PhotonCounts_scint);
			eventInformation->IncHitCerenkCount(PhotonCounts_cerenk);
			eventInformation->IncHitWLSCount(PhotonCounts_wls);



			//search if current strip has already been stored in collection
			for(unsigned k=0;k<fPMTHitCollection.size();k++){
				if((*PHC)[i]->GetStripNumber()== fPMTHitCollection[k].StripId && (*PHC)[i]->GetPMTNumber()== fPMTHitCollection[k].PMTId && (*PHC)[i]->GetPlaneNumber()== fPMTHitCollection[k].PlaneId && (*PHC)[i]->GetSuperPlaneNumber()== fPMTHitCollection[k].SuperPlaneId) {
					//hit already stored in vector for this strip
					hitInCollection= true;
					collection_index= k;
					break;
				}
			}//close for


			if(hitInCollection){
				//add info to strip vector element			
				(fPMTHitCollection[collection_index].PhotonCounts)+= PhotonCounts;
				(fPMTHitCollection[collection_index].PhotonCounts_scint)+= PhotonCounts_scint;
				(fPMTHitCollection[collection_index].PhotonCounts_cerenk)+= PhotonCounts_cerenk;
				(fPMTHitCollection[collection_index].PhotonCounts_wls)+= PhotonCounts_wls;
	
				for(unsigned k=0;k<PhotonTime.size();k++){
					(fPMTHitCollection[collection_index].PhotonTime).push_back(PhotonTime[k]);	
					(fPMTHitCollection[collection_index].PhotonEnergy).push_back(PhotonEnergy[k]);	
	
					TVector3 thisPosition= TVector3(PhotonPosition[k].x()/CLHEP::cm,PhotonPosition[k].y()/CLHEP::cm,PhotonPosition[k].z()/CLHEP::cm);
					(fPMTHitCollection[collection_index].PhotonPosition).push_back(thisPosition);
				}//close for

			}
			else{
				//new hit vector
				//set tree branch values
				currentHit.StripId= StripId;
				currentHit.PlaneId= PlaneId;	
				currentHit.SuperPlaneId= SuperPlaneId;	
				currentHit.PMTId= PMTId;
				currentHit.Position= Position;
				currentHit.PhotonCounts= PhotonCounts;
				currentHit.PhotonCounts_scint= PhotonCounts_scint;
				currentHit.PhotonCounts_cerenk= PhotonCounts_cerenk;
				currentHit.PhotonCounts_wls= PhotonCounts_wls;

				for(unsigned k=0;k<PhotonTime.size();k++){
					(currentHit.PhotonTime).push_back(PhotonTime[k]);	
					(currentHit.PhotonEnergy).push_back(PhotonEnergy[k]);	

					TVector3 thisPosition= TVector3(PhotonPosition[k].x()/CLHEP::cm,PhotonPosition[k].y()/CLHEP::cm,PhotonPosition[k].z()/CLHEP::cm);
					(currentHit.PhotonPosition).push_back(thisPosition);
				}//close for

				fPMTHitCollection.push_back(currentHit);//add this strip in the collection
			}//close else
	
      
      reconPos+=(*PHC)[i]->GetPMTPos()*(*PHC)[i]->GetPhotonCount();
      if((*PHC)[i]->GetPhotonCount()>=pmtThreshold){
				eventInformation->IncPMTSAboveThreshold();
      }
      else{//wasnt above the threshold, turn it back off
				(*PHC)[i]->SetDrawit(false);
      }
      
			cout<<"pmtNo="<< (*PHC)[i]->GetPMTNumber()<<"  stripNo="<< (*PHC)[i]->GetStripNumber()<<"  planeNo="<<(*PHC)[i]->GetPlaneNumber()<<"  superplaneNo="<<(*PHC)[i]->GetSuperPlaneNumber()<<"  has "<<(*PHC)[i]->GetPhotonCount()<<" counts"<<endl;
			cout<<"**  scintillation counts==> "<<(*PHC)[i]->GetPhotonCount_scint()<<endl;
	    cout<<"**  cerenkov counts==> "<<(*PHC)[i]->GetPhotonCount_cerenk()<<endl;
			cout<<"**  wls counts==> "<<(*PHC)[i]->GetPhotonCount_wls()<<endl;		
	
		}//close for hits

    if(eventInformation->GetHitCount()>0){//dont bother unless there were hits
      reconPos/=eventInformation->GetHitCount();
		  cout << "\tReconstructed position of hits in PMT : "
	       		<< reconPos/CLHEP::mm << endl;
      eventInformation->SetReconPos(reconPos);
    }
    PHC->DrawAllHits();
   
		//## FILL StationSimData
		for(unsigned int kk=0;kk<fPMTHitCollection.size();kk++){
			TPMTHit currentPMTHit= fPMTHitCollection[kk];
			//stationSimData.AddPMTHit(currentPMTHit);
			currentParticleSimData.AddPMTHit(currentPMTHit);
		}//end loop scint hits

	}//close if PMTHits
	

	//#######################################
	//##  ADD STATION SIM DATA TO EVENT
	//#######################################
	unsigned int nHitForThisParticle= (currentParticleSimData.fScintHitCollection).size();
	
	currentParticleSimData.AddMuonCount(fMuonInNumberInEvent,std::string("in"));
  currentParticleSimData.AddEmCount(fEmInNumberInEvent,std::string("in"));
	currentParticleSimData.AddHadronCount(fHadronInNumberInEvent,std::string("in"));
	currentParticleSimData.AddMuonCount(fMuonOutNumberInEvent,std::string("out"));
  currentParticleSimData.AddEmCount(fEmOutNumberInEvent,std::string("out"));
	currentParticleSimData.AddHadronCount(fHadronOutNumberInEvent,std::string("out"));

	fStationSimData->AddSimParticle(currentParticleSimData);
	//fStationSimData->AddSimParticleInTimeBin(currentParticleSimData);

	//cout<<"G4MuonCounterRunAction(): INFO: fStationSimData nParticles="<<(fStationSimData->fParticleSimDataCollection).size()<<endl;
	currentEventSimData->AddSimStation(*fStationSimData);
	
}//close function


void G4MuonCounterRunAction::HandlePrimitiveScorers(const G4Run* aRun){

	// Set G4MuonCounterRunScorer object.
 	G4MuonCounterRunScorer* runScorer = (G4MuonCounterRunScorer*)aRun;
  // Dump all socred quantities involved in G4MuonCounterRunScorer
  //runScorer->DumpAllScorer();
  
  //---------------------------------------------
  // Dump accumulated quantities for this RUN.
  //  (Display only central region of x-y plane)
  //---------------------------------------------
  G4THitsMap<double>* muonSurfaceInCurrent = runScorer->GetHitsMap("PlaneModuleSD/muonSurfaceInCurrent");
  G4THitsMap<double>* emSurfaceInCurrent = runScorer->GetHitsMap("PlaneModuleSD/emSurfaceInCurrent");
	G4THitsMap<double>* muonSurfaceOutCurrent = runScorer->GetHitsMap("PlaneModuleSD/muonSurfaceOutCurrent");
  G4THitsMap<double>* emSurfaceOutCurrent = runScorer->GetHitsMap("PlaneModuleSD/emSurfaceOutCurrent");
	

  cout << "=============================================================" <<endl;
  cout << " Number of event processed : "<< aRun->GetNumberOfEvent() << endl;
  cout << "=============================================================" <<endl;
	
	for(int iz=0;iz<nPlanes;iz++){ 
		double* muonSurfaceInCounts= (*muonSurfaceInCurrent)[iz];
		double* emSurfaceInCounts= (*emSurfaceInCurrent)[iz];
		double* muonSurfaceOutCounts= (*muonSurfaceOutCurrent)[iz];
		double* emSurfaceOutCounts= (*emSurfaceOutCurrent)[iz];
		
		if(!muonSurfaceInCounts) muonSurfaceInCounts= new double(0.0);
		if(!emSurfaceInCounts) emSurfaceInCounts= new double(0.0);
		if(!muonSurfaceOutCounts) muonSurfaceOutCounts= new double(0.0);
		if(!emSurfaceOutCounts) emSurfaceOutCounts= new double(0.0);
		
		fMuonInNumber[iz]= (*muonSurfaceInCounts);
		fEmInNumber[iz]= (*emSurfaceInCounts);
		fMuonOutNumber[iz]= (*muonSurfaceOutCounts);
		fEmOutNumber[iz]= (*emSurfaceOutCounts);

		cout << "Plane no. " << iz << " : muons (in/out) "<< (*muonSurfaceInCounts)  << "/" << (*muonSurfaceOutCounts) << "  em (in/out) "<< (*emSurfaceInCounts) << "/" << (*muonSurfaceOutCounts) << endl;
	}//end loop planes
  cout << "============================================="<<endl;
	
	
}//close HandlePrimitiveScorers()


void G4MuonCounterRunAction::EndOfRunAction(const G4Run* aRun){

	//cout<<"G4MuonCounterRunAction::EndOfRunAction"<<endl;
	//## Handle primitive scorers
	//HandlePrimitiveScorers(aRun);
	
	/*
	cout << "=============================================================" <<endl;
  cout << " Number of event processed : "<< aRun->GetNumberOfEvent() << endl;
  cout << "=============================================================" <<endl;
	for(unsigned int i=0;i<fMuonInNumber.size();i++){ 
		cout << "Plane no. " << i << " : muons (in) "<< fMuonInNumber[i] << "  em (in) "<< fEmInNumber[i] << "  hadr (in) "<< fHadronInNumber[i]<<endl;
	}//end loop planes
  cout << "============================================="<<endl;
	*/
	
	//#######################################
	//##  ADD STATION SIM DATA TO EVENT
	//#######################################
	//## Add current station sim data to event sim data
	//## internal check if this station has been already included
	//## If so just append particle info to stored station
	//Retrieve current TStationSimData
	TStationSimData* currentStationSimData= currentEventSimData->GetSimStation(G4MuonCounterSimulator::GetCurrentDetectorStation()->GetId());
	currentStationSimData->AddMuonCount(fMuonInNumber,std::string("in"));
  currentStationSimData->AddEmCount(fEmInNumber,std::string("in"));
	currentStationSimData->AddHadronCount(fHadronInNumber,std::string("in"));
	currentStationSimData->AddMuonCount(fMuonOutNumber,std::string("out"));
  currentStationSimData->AddEmCount(fEmOutNumber,std::string("out"));
	currentStationSimData->AddHadronCount(fHadronOutNumber,std::string("out"));

  if(recorder) recorder->RecordEndOfRun(aRun);

}
