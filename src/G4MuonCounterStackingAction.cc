/**
* @file G4MuonCounterStackingAction.cc
* @class G4MuonCounterStackingAction
* @brief Handle all operations to be performed at stack level
*
* Allows to classify each new track and manipulate its status (kill,track urgent, ...), count tracks of a given type (scintillator,cerenkov,...).
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterStackingAction.hh"
#include "UserEventInformation.hh"
#include "G4MuonCounterSteppingAction.hh"
#include "G4MuonCounterRecorderBase.hh"

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

#include "TVector3.h"
#include "TMath.h"

using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterStackingAction::G4MuonCounterStackingAction(G4MuonCounterRecorderBase* r)
:recorder(r)
{
	
  G4cout<<"Create G4MuonCounterStackingAction"<<G4endl;

	Nscint=0;
	Ncerenk=0;

}

G4MuonCounterStackingAction::~G4MuonCounterStackingAction(){
	

}

G4ClassificationOfNewTrack G4MuonCounterStackingAction::ClassifyNewTrack(const G4Track * aTrack){
  
  //G4cout<<"Call ClassifyNewTrack"<<G4endl;	

  UserEventInformation* eventInformation= (UserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();



  //Count what process generated the optical photons
	// particle is optical photon
  if(aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){ 
    

		// current particle is secondary and its parent is a primary track
		// current particle is inside the scintillator
    if( aTrack->GetParentID()==1){
			
			/*
			//get here the direction of the parent track (this must be a primary!) of this secondary
			// to be implemented...see this tips:
			//    http://geant4.slac.stanford.edu/Tips/event/1.html
			//    http://geant4.slac.stanford.edu/Tips/event/3.html
			std::vector<G4ThreeVector> primaryTrackDirectionVect= eventInformation->GetPrimaryParticleDirection();
			G4ThreeVector primaryTrackDirection= primaryTrackDirectionVect[0];//TO BE REPLACED!!!
			TVector3 primaryTrackDirection_ROOT;
			primaryTrackDirection_ROOT.SetXYZ(primaryTrackDirection.x(),primaryTrackDirection.y(),primaryTrackDirection.z());		


			//get the direction of the current secondary photon
			G4ThreeVector currentTrackDirection= aTrack->GetMomentumDirection();
			TVector3 currentTrackDirection_ROOT;
			currentTrackDirection_ROOT.SetXYZ(currentTrackDirection.x(),currentTrackDirection.y(),currentTrackDirection.z());			

			//calculate and set emission angle (angle between primary and secondary directions)
			G4double EmissionAngle= currentTrackDirection_ROOT.Angle(primaryTrackDirection_ROOT);
			EmissionAngle*= 180./TMath::Pi();
			eventInformation->InsertPhotonEmissionAngle(EmissionAngle);

			//set emission time
			G4double EmissionTime= aTrack->GetGlobalTime();
			eventInformation->InsertPhotonEmissionTime(EmissionTime/ns);
			*/

	
			
			//set process creator
			if(aTrack->GetCreatorProcess()->GetProcessName()=="Scintillation") {
				//eventInformation->IncPhotonCount_Scint();
				//eventInformation->InsertProcessType(1);
				Nscint++;
			}
			else if(aTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov"){
				//eventInformation->IncPhotonCount_Ceren();
				//eventInformation->InsertProcessType(2);
				Ncerenk++;
			}
			else if(aTrack->GetCreatorProcess()->GetProcessName()=="OpWLS") {
				//eventInformation->IncPhotonCount_WLS();			
				//eventInformation->InsertProcessType(3);
			}//close else if

			
			// get module volume of the current track
			//G4VPhysicalVolume* currentModuleVolume = aTrack->GetStep()->GetPreStepPoint()->GetTouchableHandle()->GetVolume(3);
  		//G4int currentModuleCopyNo= currentModuleVolume->GetCopyNo();
			//G4cout<<"G4MuonCounterStackingAction: "<<"  currentModuleCopyNo="<< currentModuleCopyNo<<G4endl;
			
			//insert current strip number
			//eventInformation->InsertStripNumber(currentModuleCopyNo);
			//G4cout<<"fill StripNumber vector"<<G4endl;	
			
    }//close if secondary
		
		//G4cout<<"TRACKID="<<aTrack->GetTrackID()<<"  Nscint="<<Nscint<<"  Ncerenk="<<Ncerenk<<G4endl;

  }//close if optical photon
  		
	

  return fUrgent;
}

void G4MuonCounterStackingAction::NewStage(){
}

void G4MuonCounterStackingAction::PrepareNewEvent(){ 
	
	Nscint=0;
	Ncerenk=0;
}

