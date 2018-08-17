/**
* @file RunMessenger.hh
* @class RunMessenger
* @brief Create the UI commands to handle the run action parameters, i.e. from a configuration file  
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_RunMessenger_h
#define _G4MuonCounterSimulatorUSC_RunMessenger_h 1

#include "G4UImessenger.hh"
#include "G4RunMessenger.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterRunAction;

class RunMessenger: public G4UImessenger
{
  public:
		/** 
		\brief Class constructor: create UI commands
 		*/
    RunMessenger(G4MuonCounterRunAction*);
		/**
		* \brief Class destructor: delete allocated UI commands
		*/
   ~RunMessenger();
    /**
		* \brief Set (on-fly) the values parsed by the command to the detector geometry, by calling members of the detector construction class
		*/
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    G4MuonCounterRunAction* 									 theRun;
    G4UIdirectory*               runDir; 
		G4UIcmdWithABool*     storeInfoCmd;
		
};

}//close namespace

#endif

