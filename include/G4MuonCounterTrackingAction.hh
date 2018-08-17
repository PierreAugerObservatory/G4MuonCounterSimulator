/**
* @file G4MuonCounterTrackingAction.hh
* @class G4MuonCounterTrackingAction
* @brief Handle all operations to be performed at track level
*
* Manage operations before the track starts (i.e. access to primary and secondary track information, count given tracks,...), operations when the track stops.
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterTrackingAction_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterRecorderBase;
class G4MuonCounterPrimaryGenerator;

class G4MuonCounterTrackingAction : public G4UserTrackingAction {

public:  
	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4MuonCounterTrackingAction(G4MuonCounterRecorderBase*,G4MuonCounterPrimaryGenerator*);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~G4MuonCounterTrackingAction() {};


	/**
	* \brief Which process generate the opt photon
	*/
	enum processFlag { kScintillation= 1, kCerenkov= 2, kWLS= 3 };


  /**
	* \brief Operations to be done at the beginning of the track
	*/
  void PreUserTrackingAction(const G4Track*);
	/**
	* \brief Operations to be done at the end of the track
	*/
  void PostUserTrackingAction(const G4Track*);

  
private:
  G4MuonCounterRecorderBase* recorder;
  G4MuonCounterPrimaryGenerator* generator;

	//#DEBUG
	int Nscint;
	int Ncerenk;
	
	bool recordEmissionInfo;
  
};

}//close namespace
#endif
