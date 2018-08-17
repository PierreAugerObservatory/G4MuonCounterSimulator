/**
* @file SteppingMessenger.cc
* @class SteppingMessenger
* @brief Create the UI commands to handle the stepping action parameters, i.e. from a configuration file.   
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#include "SteppingMessenger.hh"
#include "G4MuonCounterSteppingAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

using namespace G4MuonCounterSimulatorUSC;

SteppingMessenger::SteppingMessenger(G4MuonCounterSteppingAction* step)
:stepping(step)
{
  oneStepPrimariesCmd = new G4UIcmdWithABool("/StripSim/oneStepPrimaries",this);
  oneStepPrimariesCmd->SetGuidance("Only allows primaries to go one step in the scintillator volume before being killed.");
}


SteppingMessenger::~SteppingMessenger(){
  delete oneStepPrimariesCmd;
}


void SteppingMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 
  if( command == oneStepPrimariesCmd ){ 
    stepping->SetOneStepPrimaries(oneStepPrimariesCmd->GetNewBoolValue(newValue));
  }
}


