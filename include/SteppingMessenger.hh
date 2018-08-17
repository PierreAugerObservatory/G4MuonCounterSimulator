/**
* @file SteppingMessenger.hh
* @class SteppingMessenger
* @brief Create the UI commands to handle the stepping action parameters, i.e. from a configuration file.  
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_SteppingMessenger_h
#define _G4MuonCounterSimulatorUSC_SteppingMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIcmdWithABool;


namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterSteppingAction;

class SteppingMessenger: public G4UImessenger
{
public:
	/** 
	\brief Class constructor: create UI commands
 	*/
  SteppingMessenger(G4MuonCounterSteppingAction*);
	/**
	* \brief Class destructor: delete allocated UI commands
	*/
  ~SteppingMessenger();
  /**
	* \brief Set (on-fly) the values parsed by the command to the detector geometry, by calling members of the detector construction class
	*/
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  G4MuonCounterSteppingAction*        stepping;
  G4UIcmdWithABool*  oneStepPrimariesCmd;
  
};

}//close namespace

#endif

