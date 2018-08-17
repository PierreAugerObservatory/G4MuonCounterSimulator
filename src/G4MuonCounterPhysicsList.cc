/**
* @file G4MuonCounterPhysicsList.cc
* @class G4MuonCounterPhysicsList
* @brief Build the physical processes and particles involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#include "G4MuonCounterPhysicsList.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"

#include "G4MuonCounterDecayPhysics.hh"
#include "G4MuonCounterBosonPhysics.hh"
#include "G4MuonCounterLeptonPhysics.hh"
#include "G4MuonCounterNeutronPhysics.hh"
#include "G4MuonCounterHadronPhysics.hh"
#include "G4MuonCounterIonPhysics.hh"
#include "G4MuonCounterOpticalPhysics.hh"

#include "G4MuonCounterSimulator.hh"

using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterPhysicsList::G4MuonCounterPhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*CLHEP::mm;
  // SetVerboseLevel(1);

	fIncludeOpticalPhysics= false;

  // Particle decays
  RegisterPhysics( new G4MuonCounterDecayPhysics("decay"));

  // Bosons (gamma + geantinos)
  RegisterPhysics( new G4MuonCounterBosonPhysics("boson"));

  // Leptons
  RegisterPhysics( new G4MuonCounterLeptonPhysics("lepton"));

  // Neutron Physics
  RegisterPhysics( new G4MuonCounterNeutronPhysics("neutron"));

  // Hadron Physics
  RegisterPhysics( new G4MuonCounterHadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new G4MuonCounterIonPhysics("ion"));

	// Optical Physics (cerenkov,scintillation, wls, ...)
  if(G4MuonCounterSimulator::OpticalOn()) RegisterPhysics( new G4MuonCounterOpticalPhysics("optical"));
}

G4MuonCounterPhysicsList::~G4MuonCounterPhysicsList()
{;}

void G4MuonCounterPhysicsList::SetCuts()
{
  // Use default cut values gamma and e processes
  SetCutsWithDefault(); 
	//DumpCutValuesTable();  
}

void G4MuonCounterPhysicsList::SetOpticalPhysics(bool value)
{
  fIncludeOpticalPhysics= value;  
}

