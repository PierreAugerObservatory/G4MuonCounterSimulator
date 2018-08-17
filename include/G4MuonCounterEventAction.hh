/**
* @file G4MuonCounterEventAction.hh
* @class G4MuonCounterEventAction
* @brief Handle all operations to be performed at event level
*
* Manage operations before the event starts (i.e. initialize hit structures, ...), operations after the event ends (i.e. access to simulation hits and tracks, hit manipulation,...).
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterEventAction_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterEventAction_h 1

#include "EventMessenger.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

class G4Event;

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterRecorderBase;
class G4MuonCounterRunAction;

class G4MuonCounterEventAction : public G4UserEventAction
{
public:
	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4MuonCounterEventAction(G4MuonCounterRecorderBase*,G4MuonCounterRunAction*);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~G4MuonCounterEventAction();
  
public:

	/**
	* \brief Operations to be done before event starts
	*/
  void BeginOfEventAction(const G4Event*);
	/**
	* \brief Operations to be done when event ends
	*/
  void EndOfEventAction(const G4Event*);
  /**
	* \brief Sets the save threshold for the random number seed. If the number of photons generated in an event is lower than this, then save the seed for this event in a file called runXXXevtXXX.rndm
	*/
  void SetSaveThreshold(int save);
	/**
	* \brief Set verbosity of the events
	*/
  void SetEventVerbose(int v){verbose=v;}
	/**
	* \brief Set a count threshold for the PMTs
	*/
  void SetPMTThreshold(int t){pmtThreshold=t;}
	/**
	* \brief Force the drawing of the produced photons
	*/
  void SetForceDrawPhotons(bool b){forcedrawphotons=b;}
	/**
	* \brief Exclude the drawing of the produced photons
	*/
  void SetForceDrawNoPhotons(bool b){forcenophotons=b;}

private:
	/**
	* \brief Recorder base for this class
	*/
  G4MuonCounterRecorderBase* recorder;
	/**
	* \brief Event messenger for this class
	*/
  EventMessenger* eventMessenger;
	/**
	* \brief Run action for this class
	*/
  G4MuonCounterRunAction* Run;
	/**
	* \brief Save threshold flag
	*/
  int saveThreshold;

	/**
	* \brief Scintillator hit collection flag
	*/
	int scintCollID;
	/**
	* \brief PMT hit collection flag
	*/
  int pmtCollID;
  /**
	* \brief Verbosity of events
	*/
  int verbose;	
	/**
	* \brief PMT count threshold
	*/
  int pmtThreshold;
	/**
	* \brief Flag to force the drawing of produced photons
	*/
  bool forcedrawphotons;
	/**
	* \brief Flag to exclude the drawing of produced photons
	*/
  bool forcenophotons;
  
  

};

}//close namespace

#endif

    
