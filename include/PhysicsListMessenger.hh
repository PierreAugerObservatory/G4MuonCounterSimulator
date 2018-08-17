/**
* @file PhysicsListMessenger.hh
* @class PhysicsListMessenger
* @brief Create the UI commands to handle the physics process parameters, i.e. from a configuration file 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef PhysicsListMessenger_h
#define PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4MuonCounterPhysicsList;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;


class PhysicsListMessenger: public G4UImessenger
{
  public:

	  /** 
		\brief Class constructor: create UI commands
 		*/
    PhysicsListMessenger(G4MuonCounterPhysicsList* );
		/**
	  * \brief Class destructor: delete allocated UI commands
	  */
   ~PhysicsListMessenger();
    /**
		* \brief Set (on-fly) the values parsed by the command to the physics list, by calling members of the physics list class
		*/
    void SetNewValue(G4UIcommand*, G4String);
    
  private:  
    G4MuonCounterPhysicsList*     pPhysicsList;
    
    G4UIdirectory*        physDir;

/*
    G4UIcmdWithAnInteger* verboseCmd;
    G4UIcmdWithAnInteger* cerenkovCmd;
		G4UIcmdWithoutParameter*     activateOpticalProcessesCmd;
		G4UIcmdWithoutParameter*     inactivateOpticalProcessesCmd;
*/
};


#endif

