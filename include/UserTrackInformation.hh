/**
* @file UserTrackInformation.hh
* @class UserTrackInformation
* @brief User-defined container of track information produced in the event
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_UserTrackInformation_h
#define _G4MuonCounterSimulatorUSC_UserTrackInformation_h 1

#include "G4VUserTrackInformation.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Allocator.hh"


namespace G4MuonCounterSimulatorUSC {

enum TrackStatus { active=1, hitPMT=2, internalReflected=3,absorbed=4, rayleighScattered=5, boundaryAbsorbed=8, inactive=14};

/*TrackStatus:
  active: still being tracked
  hitPMT: stopped by being detected in a PMT
  absorbed: stopped by being absorbed with G4OpAbsorbtion
  boundaryAbsorbed: stopped by being aborbed with G4OpAbsorbtion
  hitSphere: track hit the sphere at some point
  inactive: track is stopped for some reason
   -This is the sum of all stopped flags so can be used to remove stopped flags
  
 */

class UserTrackInformation : public G4VUserTrackInformation
{
public:

	/** 
	\brief Class constructor: initialize structures.
 	*/
  UserTrackInformation();
	/** 
	\brief Overloaded class constructor: initialize structures.
 	*/
  UserTrackInformation(const G4Track* aTrack);
	/** 
	\brief Overloaded class constructor: initialize structures.
 	*/
  UserTrackInformation(const UserTrackInformation* aTrackInfo);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~UserTrackInformation();

	/**
	* \brief Allocate new class
	*/
	inline void *operator new(size_t);
	/**
	* \brief Delete allocated class
	*/
  inline void operator delete(void *aTrackInfo);
	/**
	* \brief Tells if two tracks are equals
	*/
  inline int operator ==(const UserTrackInformation& right) const {return (this==&right);}
	/**
	* \brief Get track ID
	*/
	inline int GetOriginalTrackID() const {return originalTrackID;}
	/**
	* \brief Get particle type
	*/
  inline G4ParticleDefinition* GetOriginalParticle() const {return particleDefinition;}
	/**
	* \brief Get particle position
	*/
  inline G4ThreeVector GetOriginalPosition() const {return originalPosition;}
	/**
	* \brief Get particle momentum
	*/
  inline G4ThreeVector GetOriginalMomentum() const {return originalMomentum;}
	/**
	* \brief Get particle energy
	*/
  inline double GetOriginalEnergy() const {return originalEnergy;}
	/**
	* \brief Get particle time
	*/
  inline double GetOriginalTime() const {return originalTime;}


  
  //Sets the track status to s (does not check validity of flags)
	/**
	* \brief Set the track status flag (active,absorbed,hitPMT,...)
	*/
  void SetTrackStatusFlags(int s){status=s;}
  //Does a smart add of track status flags (disabling old flags that conflict)
  //If s conflicts with itself it will not be detected
	/**
	* \brief Add the track status flag (active,absorbed,hitPMT,...)
	*/
  void AddTrackStatusFlag(int s);
  /**
	* \brief Get the track status flag (active,absorbed,hitPMT,...)
	*/
  int GetTrackStatus()const {return status;}
  /**
	* \brief Increment number of tracks being reflected on the walls (+1)
	*/
  void IncReflections(){reflections++;}
	/**
	* \brief Get number of tracks being reflected on the walls
	*/
  int GetReflectionCount()const {return reflections;}
	/**
	* \brief Force on/off the drawing of trajectories
	*/
  void SetForceDrawTrajectory(bool b){forcedraw=b;}
	/**
	* \brief Get the force drawing flag
	*/
  bool GetForceDrawTrajectory(){return forcedraw;}
	/**
	* \brief Print user information
	*/
  inline void Print()const;


private:
  int status;
  int reflections;
  bool forcedraw;

	int                 originalTrackID;
  G4ParticleDefinition* particleDefinition;
  G4ThreeVector         originalPosition;
  G4ThreeVector         originalMomentum;
  double              originalEnergy;
  double              originalTime;


};

}//close namespace


extern G4Allocator<G4MuonCounterSimulatorUSC::UserTrackInformation> aTrackInformationAllocator;


namespace G4MuonCounterSimulatorUSC {

inline void* UserTrackInformation::operator new(size_t){ 
	void* aTrackInfo;
  aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void UserTrackInformation::operator delete(void *aTrackInfo){ 
	aTrackInformationAllocator.FreeSingle((UserTrackInformation*)aTrackInfo);
}

}//close namespace



#endif
