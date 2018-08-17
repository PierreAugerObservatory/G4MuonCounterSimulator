/**
* @file EventMessenger.hh
* @class EventMessenger
* @brief Create the UI commands to handle the event action parameters, i.e. from a configuration file 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_EventMessenger_h
#define _G4MuonCounterSimulatorUSC_EventMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"


class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterEventAction;

class EventMessenger: public G4UImessenger
{
public:
	/** 
	\brief Class constructor: create UI commands
 	*/
  EventMessenger(G4MuonCounterEventAction*);
	/**
	* \brief Class destructor: delete allocated UI commands
	*/
  ~EventMessenger();
  /**
	* \brief Set (on-fly) the values parsed by the command to the detector geometry, by calling members of the detector construction class
	*/
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  G4MuonCounterEventAction*        Event;
  G4UIcmdWithAnInteger*  saveThresholdCmd;
  G4UIcmdWithAnInteger*  verboseCmd;
  G4UIcmdWithAnInteger*  pmtThresholdCmd;
  G4UIcmdWithABool*      forceDrawPhotonsCmd;
  G4UIcmdWithABool*      forceDrawNoPhotonsCmd;
};

}//close namespace

#endif

