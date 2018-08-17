/**
* @file UserTrackInformation.cc
* @class UserTrackInformation
* @brief User-defined container of track information produced in the event
* @author Dr. Simone Riggi
* @date 05/04/2010
*/


#include "UserTrackInformation.hh"

using namespace G4MuonCounterSimulatorUSC;

G4Allocator<UserTrackInformation> aTrackInformationAllocator;


UserTrackInformation::UserTrackInformation()
  :status(active),reflections(0),forcedraw(false)
{
	originalTrackID = 0;
  particleDefinition = 0;
  originalPosition = G4ThreeVector(0.,0.,0.);
  originalMomentum = G4ThreeVector(0.,0.,0.);
  originalEnergy = 0.;
  originalTime = 0.;

}//close constructor

UserTrackInformation::UserTrackInformation(const G4Track* aTrack){
	originalTrackID = aTrack->GetTrackID();
  particleDefinition = aTrack->GetDefinition();
  originalPosition = aTrack->GetPosition();
  originalMomentum = aTrack->GetMomentum();
  originalEnergy = aTrack->GetTotalEnergy();
  originalTime = aTrack->GetGlobalTime();
}

UserTrackInformation::UserTrackInformation(const UserTrackInformation* aTrackInfo){
    
	originalTrackID = aTrackInfo->originalTrackID;
  particleDefinition = aTrackInfo->particleDefinition;
  originalPosition = aTrackInfo->originalPosition;
  originalMomentum = aTrackInfo->originalMomentum;
  originalEnergy = aTrackInfo->originalEnergy;
  originalTime = aTrackInfo->originalTime;
}


UserTrackInformation::~UserTrackInformation(){

}//close destructor


void UserTrackInformation::AddTrackStatusFlag(int s)
{
  if(s&active) //track is now active
    status&=~inactive; //remove any flags indicating it is inactive 
  else if(s&inactive) //track is now inactive
    status&=~active; //remove any flags indicating it is active
  status|=s; //add new flags
}


void UserTrackInformation::Print() const{
    G4cout << "Original track ID " << originalTrackID 
           << " at " << originalPosition << G4endl;
}


