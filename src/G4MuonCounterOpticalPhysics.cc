/**
* @file G4MuonCounterOpticalPhysics.cc
* @class G4MuonCounterOpticalPhysics
* @brief Build the optical physical processes and optical particles involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#include "G4MuonCounterOpticalPhysics.hh"
#include "G4MuonCounterSimulator.hh"

#include "G4LossTableManager.hh"
//#include "G4EmSaturation.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4OpticalPhoton.hh"

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4MuonCounterSimulator.hh"

using namespace G4MuonCounterSimulatorUSC;


G4MuonCounterOpticalPhysics::G4MuonCounterOpticalPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name)
{
  //  new OpticalPhysicsMessenger(this);
}

G4MuonCounterOpticalPhysics::~G4MuonCounterOpticalPhysics(){
}


void G4MuonCounterOpticalPhysics::ConstructParticle(){
  G4OpticalPhoton::OpticalPhotonDefinition();
}


void G4MuonCounterOpticalPhysics::ConstructProcess(){

	G4cout<<"*** Building optical processes ***"<<G4endl;
 
  theScintillationProcess = new G4Scintillation();
  theCerenkovProcess=new G4Cerenkov();
  theAbsorptionProcess=new G4OpAbsorption();
  theRayleighScatteringProcess=new G4OpRayleigh();
  theBoundaryProcess=new G4OpBoundaryProcess();
  theWLSProcess=new G4OpWLS();

	//theCerenkovProcess->DumpPhysicsTable();
  //theScintillationProcess->DumpPhysicsTable();
  //theScintillationProcess->DumpInfo();
  //theAbsorptionProcess->DumpPhysicsTable();
  //theRayleighScatteringProcess->DumpPhysicsTable();

  //###############################
  //##  CHERENKOV PROCESS    ######
  //###############################
  theCerenkovProcess->SetMaxNumPhotonsPerStep(300);//300 original value
  theCerenkovProcess->SetTrackSecondariesFirst(true);

	//#########################
  //##  WLS PROCESS    ######
  //#########################
  theWLSProcess->UseTimeProfile("delta");
  //theWLSProcess->UseTimeProfile("exponential");
 
	//###################################
  //##  SCINTILLATION PROCESS    ######
  //###################################
  theScintillationProcess->SetScintillationYieldFactor(1);//1.0
  theScintillationProcess->SetTrackSecondariesFirst(true);
  //theScintillationProcess->SetScintillationExcitationRatio(0.0);
	
	//Use Birks Correction in the Scintillation process
	//G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  //theScintillationProcess->AddSaturation(emSaturation);

	//##############################
  //##  BOUNDARY PROCESS    ######
  //##############################
  G4OpticalSurfaceModel themodel = unified;
  theBoundaryProcess->SetModel(themodel);


	//##############################
  //## ADD PROCESSES        ######
  //##############################
  G4ProcessManager * pManager = 0;
  
  pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

	//add processes
  //if(G4MuonCounterSimulator::AbsorptionOn()) pManager->AddDiscreteProcess(theAbsorptionProcess);
  //if(G4MuonCounterSimulator::RayleighOn()) pManager->AddDiscreteProcess(theRayleighScatteringProcess);
  //if(G4MuonCounterSimulator::BoundaryOn()) pManager->AddDiscreteProcess(theBoundaryProcess);
  //if(G4MuonCounterSimulator::WLSOn()) pManager->AddDiscreteProcess(theWLSProcess);
	pManager->AddDiscreteProcess(theAbsorptionProcess);
  pManager->AddDiscreteProcess(theRayleighScatteringProcess);
  pManager->AddDiscreteProcess(theBoundaryProcess);
  pManager->AddDiscreteProcess(theWLSProcess);

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    pManager = particle->GetProcessManager();
    //if(theCerenkovProcess->IsApplicable(*particle) && G4MuonCounterSimulator::CherenkovOn()){
		if(theCerenkovProcess->IsApplicable(*particle)){
      pManager->AddProcess(theCerenkovProcess);
      pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
    }
    //if(theScintillationProcess->IsApplicable(*particle) && G4MuonCounterSimulator::ScintillationOn()){
		if(theScintillationProcess->IsApplicable(*particle)){
      pManager->AddProcess(theScintillationProcess);
      pManager->SetProcessOrderingToLast(theScintillationProcess,idxAtRest);
      pManager->SetProcessOrderingToLast(theScintillationProcess,idxPostStep);
    }
  }

}//close function


void G4MuonCounterOpticalPhysics::SetVerbose(G4int verbose){

  theCerenkovProcess->SetVerboseLevel(verbose);
  theScintillationProcess->SetVerboseLevel(verbose);
  theAbsorptionProcess->SetVerboseLevel(verbose);
  theRayleighScatteringProcess->SetVerboseLevel(verbose);
  theBoundaryProcess->SetVerboseLevel(verbose);  
  theWLSProcess->SetVerboseLevel(verbose);
}


void G4MuonCounterOpticalPhysics::SetScintYieldFactor(G4double yf){
  if(theScintillationProcess) theScintillationProcess->SetScintillationYieldFactor(yf);
}

void G4MuonCounterOpticalPhysics::SetNbOfPhotonsCerenkov(G4int MaxNumber){  
  theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumber);
}


