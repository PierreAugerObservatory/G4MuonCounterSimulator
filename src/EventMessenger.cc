/**
* @file EventMessenger.cc
* @class EventMessenger
* @brief Create the UI commands to handle the event action parameters, i.e. from a configuration file  
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#include "EventMessenger.hh"
#include "G4MuonCounterEventAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

using namespace G4MuonCounterSimulatorUSC;

EventMessenger::EventMessenger(G4MuonCounterEventAction* event)
:Event(event)
{
  saveThresholdCmd = new G4UIcmdWithAnInteger("/StripSim/saveThreshold",this);
  saveThresholdCmd->SetGuidance("Set the photon count threshold for saving the random number seed for an event.");
  saveThresholdCmd->SetParameterName("photons",true);
  saveThresholdCmd->SetDefaultValue(4500);
  saveThresholdCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  verboseCmd = new G4UIcmdWithAnInteger("/StripSim/eventVerbose",this);
  verboseCmd->SetGuidance("Set the verbosity of event data.");
  verboseCmd->SetParameterName("verbose",true);
  verboseCmd->SetDefaultValue(1);

  pmtThresholdCmd = new G4UIcmdWithAnInteger("/StripSim/pmtThreshold",this);
  pmtThresholdCmd->SetGuidance("Set the pmtThreshold (in # of photons)");

  forceDrawPhotonsCmd=new G4UIcmdWithABool("/StripSim/forceDrawPhotons",this);
  forceDrawPhotonsCmd->SetGuidance("Force drawing of photons.");
  forceDrawPhotonsCmd
    ->SetGuidance("(Higher priority than /StripSim/forceDrawNoPhotons)");

  forceDrawNoPhotonsCmd=new G4UIcmdWithABool("/StripSim/forceDrawNoPhotons",this);
  forceDrawNoPhotonsCmd->SetGuidance("Force no drawing of photons.");
  forceDrawNoPhotonsCmd
    ->SetGuidance("(Lower priority than /StripSim/forceDrawPhotons)");
}


EventMessenger::~EventMessenger(){
  delete saveThresholdCmd;
  delete verboseCmd;
  delete pmtThresholdCmd;
  delete forceDrawPhotonsCmd;
  delete forceDrawNoPhotonsCmd;
}


void EventMessenger::SetNewValue(G4UIcommand* command, G4String newValue){ 
  if( command == saveThresholdCmd ){ 
    Event->SetSaveThreshold(saveThresholdCmd->GetNewIntValue(newValue));
  }
  else if( command == verboseCmd ){
    Event->SetEventVerbose(verboseCmd->GetNewIntValue(newValue));
  }
  else if( command == pmtThresholdCmd ){
    Event->SetPMTThreshold(pmtThresholdCmd->GetNewIntValue(newValue));
  }
  else if(command == forceDrawPhotonsCmd){
    Event->SetForceDrawPhotons(forceDrawPhotonsCmd
				  ->GetNewBoolValue(newValue));
  }
  else if(command == forceDrawNoPhotonsCmd){
    Event->SetForceDrawNoPhotons(forceDrawNoPhotonsCmd
				  ->GetNewBoolValue(newValue));
    G4cout<<"TEST"<<G4endl;
  }
}


