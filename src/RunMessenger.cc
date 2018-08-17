/**
* @file RunMessenger.cc
* @class RunMessenger
* @brief Create the UI commands to handle the run action parameters, i.e. from a configuration file  
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "RunMessenger.hh"

#include "G4MuonCounterRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

using namespace G4MuonCounterSimulatorUSC;

RunMessenger::RunMessenger(G4MuonCounterRunAction* Run)
:theRun(Run)
{
  runDir = new G4UIdirectory("/StripSimulation/run/");
  runDir->SetGuidance("RunAction control");
   

	storeInfoCmd = new G4UIcmdWithABool("/StripSimulation/run/FullStorage",this);
  storeInfoCmd->SetGuidance("Enable/Disable the full storage for simulation output");
  storeInfoCmd->SetParameterName("FullStorage",true,true);
	storeInfoCmd->SetDefaultValue(true);
  
}


RunMessenger::~RunMessenger(){

	delete storeInfoCmd;  
  delete runDir;
	
	
}


void RunMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
 
	if (command == storeInfoCmd){
    theRun->StoreFullInfo(storeInfoCmd->GetNewBoolValue(newValue));
  }

}


