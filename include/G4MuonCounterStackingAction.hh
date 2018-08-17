/**
* @file G4MuonCounterStackingAction.hh
* @class G4MuonCounterStackingAction
* @brief Handle all operations to be performed at stack level
*
* Allows to classify each new track and manipulate its status (kill,track urgent, ...), count tracks of a given type (scintillator,cerenkov,...).
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterStackingAction_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterRecorderBase;

class G4MuonCounterStackingAction : public G4UserStackingAction
{
public:
	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4MuonCounterStackingAction(G4MuonCounterRecorderBase*);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~G4MuonCounterStackingAction();
  /**
	* \brief Classify new tracks being created 
	*/
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
	/**
	* \brief ### To be documented ! ###
	*/
  virtual void NewStage();
	/**
	* \brief Prepare operations to be done for new events
	*/
  virtual void PrepareNewEvent();
  
private:

  G4MuonCounterRecorderBase* recorder;

	//## DEBUG
	int Nscint;
	int Ncerenk;


};

}//close namespace

#endif
