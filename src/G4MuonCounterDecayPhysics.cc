/**
* @file DecayPhysics.cc
* @class DecayPhysics
* @brief Build decay physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#include "G4MuonCounterDecayPhysics.hh"
#include "G4ProcessManager.hh"
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterDecayPhysics::G4MuonCounterDecayPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{;}

G4MuonCounterDecayPhysics::~G4MuonCounterDecayPhysics()
{;}

void G4MuonCounterDecayPhysics::ConstructParticle()
{;}

void G4MuonCounterDecayPhysics::ConstructProcess()
{
  // Add Decay Process
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (fDecayProcess.IsApplicable(*particle)) { 
      pmanager->AddProcess(&fDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(&fDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(&fDecayProcess, idxAtRest);
    }
  }



	// Declare radioactive decay to the GenericIon in the IonTable.
  const G4IonTable *theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();

  for (G4int i=0; i<theIonTable->Entries(); i++) {
  	G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
    G4String particleType = theIonTable->GetParticle(i)->GetParticleType();
      
    if (particleName == "GenericIon") {
	  	G4ProcessManager* pmanager = 
	    theIonTable->GetParticle(i)->GetProcessManager();
	  	//pmanager->SetVerboseLevel(VerboseLevel);
	  	pmanager ->AddProcess(theRadioactiveDecay);
	  	pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
	  	pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
		} 
	}


}//close ConstructProcess()


