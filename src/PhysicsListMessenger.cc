/**
* @file PhysicsListMessenger.cc
* @class PhysicsListMessenger
* @brief Create the UI commands to handle the physics process parameters, i.e. from a configuration file 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#include "PhysicsListMessenger.hh"

#include "G4MuonCounterPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithoutParameter.hh"


PhysicsListMessenger::PhysicsListMessenger(G4MuonCounterPhysicsList* pPhys)
:pPhysicsList(pPhys)
{
  
  physDir = new G4UIdirectory("/StripSimulation/phys/");
  physDir->SetGuidance("PhysicsList control");
 /*
  verboseCmd = new G4UIcmdWithAnInteger("/StripSimulation/phys/verbose",this);  
  verboseCmd->SetGuidance("set verbose for physics processes");
  verboseCmd->SetParameterName("verbose",true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("verbose>=0");
  verboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cerenkovCmd = new G4UIcmdWithAnInteger("/StripSimulation/phys/cerenkovMaxPhotons",this);  
  cerenkovCmd->SetGuidance("set max nb of photons per step");
  cerenkovCmd->SetParameterName("MaxNumber",false);
  cerenkovCmd->SetRange("MaxNumber>=0");
  cerenkovCmd->AvailableForStates(G4State_Idle);


	activateOpticalProcessesCmd = new G4UIcmdWithoutParameter("/StripSimulation/phys/activateOpticalProcesses",this);
  activateOpticalProcessesCmd->SetGuidance("Activate optical physics processes.");
	activateOpticalProcessesCmd->AvailableForStates(G4State_PreInit,G4State_Idle, G4State_GeomClosed, G4State_EventProc);  

	inactivateOpticalProcessesCmd = new G4UIcmdWithoutParameter("/StripSimulation/phys/inactivateOpticalProcesses",this);
  inactivateOpticalProcessesCmd->SetGuidance("Inactivate optical physics processes.");
	inactivateOpticalProcessesCmd->AvailableForStates(G4State_PreInit,G4State_Idle, G4State_GeomClosed, G4State_EventProc); 
*/
}


PhysicsListMessenger::~PhysicsListMessenger(){
	delete physDir;
/*
  delete verboseCmd;
  delete cerenkovCmd;
	delete activateOpticalProcessesCmd;
	delete inactivateOpticalProcessesCmd;
*/
}


void PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 
      /*
  if(command == verboseCmd){ 
		pPhysicsList->SetVerbose(verboseCmd->GetNewIntValue(newValue));
	}
   
  if(command == cerenkovCmd){
		pPhysicsList->SetNbOfPhotonsCerenkov(cerenkovCmd->GetNewIntValue(newValue));
	}
	if(command==activateOpticalProcessesCmd){
	  pPhysicsList->SetOpticalPhysics(true);
	}
	if(command==inactivateOpticalProcessesCmd){
		pPhysicsList->SetOpticalPhysics(false);
	}
*/

}


