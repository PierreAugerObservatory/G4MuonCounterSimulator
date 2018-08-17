/**
* @file G4MuonCounterEventAction.cc
* @class G4MuonCounterEventAction
* @brief Handle all operations to be performed at event level
*
* Manage operations before the event starts (i.e. initialize hit structures, ...), operations after the event ends (i.e. access to simulation hits and tracks, hit manipulation,...).
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterEventAction.hh"
#include "G4ScintillatorHit.hh"
#include "G4PMTHit.hh"
#include "UserEventInformation.hh"
#include "G4MuonCounterTrajectory.hh"
#include "G4MuonCounterRecorderBase.hh"
#include "G4MuonCounterRunAction.hh"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"

#include "TMath.h"

#include <cstddef>
#include <iostream>
#include <sstream>

using namespace std;
using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterEventAction::G4MuonCounterEventAction(G4MuonCounterRecorderBase* r,G4MuonCounterRunAction* run)
  :recorder(r),Run(run),saveThreshold(0),scintCollID(-1),pmtCollID(-1),verbose(0),
   pmtThreshold(1),forcedrawphotons(false),forcenophotons(false)
{
  
  eventMessenger= new EventMessenger(this);
  
  
}
 

G4MuonCounterEventAction::~G4MuonCounterEventAction(){
	//delete [] pmtTimeStruct.pmtPhotonTime;
  //pmtTimeStruct.pmtPhotonTime = NULL;
}


void G4MuonCounterEventAction::BeginOfEventAction(const G4Event* anEvent){
  
	//cout<<"G4MuonCounterEventAction::BeginOfEventAction()"<<endl;	


  //##########################################################
  //##  Actions to be done at the beginning of event sim: 
  //##     - create a user event info class,
  //##     - initialize hit collections
  //##     - create/append to a TTree/TH1 the event output? done in the RunAction
  //##########################################################
  
  
  //New event, add the user information object
  G4EventManager::GetEventManager()->SetUserInformation(new UserEventInformation);
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(scintCollID<0) scintCollID=SDman->GetCollectionID("scintHitCollection");
  if(pmtCollID<0) pmtCollID=SDman->GetCollectionID("pmtHitCollection");
  if(recorder)recorder->RecordBeginOfEvent(anEvent);

 
}//close G4MuonCounterEventAction::BeginOfEventAction()
 

void G4MuonCounterEventAction::EndOfEventAction(const G4Event* anEvent){
  
  //cout<<"G4MuonCounterEventAction::EndOfEventAction()"<<endl;
  
  //##########################################################
  //##  Actions to be done at the end of event sim: 
  //##     - store sim data into ROOT file output ==> done in RunAction class
  //##     - access to the user event info class,
  //##     - access to the created tracks and draw them
  //##     - access to the created hits
  //##########################################################

  UserEventInformation* eventInformation= (UserEventInformation*)anEvent->GetUserInformation();
  G4TrajectoryContainer* trajectoryContainer=anEvent->GetTrajectoryContainer();
  
  int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // extract the trajectories and draw them
  if (G4VVisManager::GetConcreteInstance()){
    for (int i=0; i<n_trajectories; i++){ 
      G4MuonCounterTrajectory* trj = (G4MuonCounterTrajectory*)((*(anEvent->GetTrajectoryContainer()))[i]);
      if(trj->GetParticleName()=="opticalphoton"){
				trj->SetForceDrawTrajectory(forcedrawphotons);
				trj->SetForceNoDrawTrajectory(forcenophotons);
      }
      
      //if(i>0) trj->SetOptPhotonTrackColor(G4Colour(0.,0.,1.));//set to blue instead of dafault green ### ONLY TO DISTINGUISH FEW TRACKS
      trj->DrawTrajectory(50);
    }
  }
  
	//Fill event info via run
	Run->FillEventInfo();


 	/*
  //End of event output. later to be controlled by a verbose level
  cout << "\tNumber of photons that hit PMTs in this event : "
	       << eventInformation->GetHitCount() << endl;
  cout << "\tNumber of photons that hit PMTs in this event (scint): "
	       << eventInformation->GetHitScintCount() << endl;
	cout << "\tNumber of photons that hit PMTs in this event (cerenk): "
	       << eventInformation->GetHitCerenkCount() << endl;
  cout << "\tNumber of photons that hit PMTs in this event (wls): "
	       << eventInformation->GetHitWLSCount() << endl;
  cout << "\tNumber of PMTs above threshold("<<pmtThreshold<<") : "
	       << eventInformation->GetPMTSAboveThreshold() << endl;
  cout << "\tNumber of photons produced by scintillation in this event : "
	       << eventInformation->GetPhotonCount_Scint() << endl;
  cout << "\tNumber of photons produced by cerenkov in this event : "
	       << eventInformation->GetPhotonCount_Ceren() << endl;
  cout << "\tNumber of photons produced by wls in this event : "
	       << eventInformation->GetPhotonCount_WLS() << endl;
  cout << "\tNumber of photons absorbed (OpAbsorption) in this event : "
	       << eventInformation->GetAbsorptionCount() << endl;
  cout << "\tNumber of photons absorbed at boundaries (OpBoundary) in "
	       << "this event : " << eventInformation->GetBoundaryAbsorptionCount() 
	       << endl;
  cout << "\tNumber of photons at boundary step too small in this event : "
	       << eventInformation->GetBoundaryStepTooSmallCount() << endl;
	cout << "\tNumber of photons killed : "
	       << eventInformation->GetPhotonCount_killed() << endl;
  cout << "Unacounted for photons in this event : " 
	       << (eventInformation->GetPhotonCount_Scint() + 
	           eventInformation->GetPhotonCount_Ceren() -
						 eventInformation->GetHitCount() - 
	           eventInformation->GetAbsorptionCount() -     
	           eventInformation->GetBoundaryAbsorptionCount()-
						 eventInformation->GetPhotonCount_killed())
	       << endl;
	*/
		
  //If we have set the flag to save 'special' events, save here
  if(saveThreshold&&eventInformation->GetPhotonCount() <= saveThreshold) G4RunManager::GetRunManager()->rndmSaveThisEvent();

  if(recorder)recorder->RecordEndOfEvent(anEvent);
  
}

void G4MuonCounterEventAction::SetSaveThreshold(int save){
/*Sets the save threshold for the random number seed. If the number of photons
generated in an event is lower than this, then save the seed for this event
in a file called run###evt###.rndm
*/
	//cout<<"G4MuonCounterEventAction saveThreshold"<<endl;

  saveThreshold=save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/random/setSavingFlag 1");
}




