/**
* @file G4MuonCounterIonPhysics.cc
* @class G4MuonCounterIonPhysics
* @brief Build ion particles and their physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterIonPhysics.hh"
#include "G4Version.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonConstructor.hh"

// processes
#if (G4VERSION_NUMBER>=930)
	#include "G4hMultipleScattering.hh"
#else
	#include "G4MultipleScattering.hh"
#endif

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// models
#include "G4LElastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"

#include "G4MuonCounterSimulator.hh"

using namespace G4MuonCounterSimulatorUSC;


G4MuonCounterIonPhysics::G4MuonCounterIonPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{;}


G4MuonCounterIonPhysics::~G4MuonCounterIonPhysics()
{;}


void G4MuonCounterIonPhysics::ConstructParticle()
{ 
  // Construct light ions (d, t, 3He, alpha, and generic ion)
  G4IonConstructor ionConstruct;
  ionConstruct.ConstructParticle();
}


void G4MuonCounterIonPhysics::ConstructProcess()
{
  // Hadronic Elastic Process and Model (for all ions except generic ion)

  G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();
  G4LElastic* elasticModel = new G4LElastic();
  elasticProcess->RegisterMe(elasticModel);

  // Hadronic inelastic models

  G4ProcessManager * pManager = 0;
  
  //#########################
  //##  Deuteron           ##
  //#########################

  pManager = G4Deuteron::Deuteron()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif

  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
	//pManager->AddProcess(new G4ionIonisation(), -1, 2, 2);//added

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4DeuteronInelasticProcess* dinelProc = new G4DeuteronInelasticProcess();
  dinelProc->RegisterMe( new G4LEDeuteronInelastic() );
  pManager->AddDiscreteProcess(dinelProc);


  //#########################
  //##  Triton             ##
  //#########################

  pManager = G4Triton::Triton()->GetProcessManager(); 

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif

  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
	//pManager->AddProcess(new G4ionIonisation(), -1, 2, 2);//added

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4TritonInelasticProcess* tinelProc = new G4TritonInelasticProcess();
  tinelProc->RegisterMe( new G4LETritonInelastic() );
  pManager->AddDiscreteProcess(tinelProc);

  //#########################
  //##  He3                ##
  //#########################
  pManager = G4He3::He3()->GetProcessManager(); 

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif

  //pManager->AddProcess(new G4hIonisation(), -1, 2, 2);
	pManager->AddProcess(new G4ionIonisation(), -1, 2, 2);//added

  // hadron inelastic

  //#########################
  //##  Alpha              ##
  //#########################

  pManager = G4Alpha::Alpha()->GetProcessManager(); 

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  //pManager->AddProcess(new G4hIonisation(), -1, 2, 2);
  pManager->AddProcess(new G4ionIonisation(), -1, 2, 2);//added

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AlphaInelasticProcess* ainelProc = new G4AlphaInelasticProcess();
  ainelProc->RegisterMe( new G4LEAlphaInelastic() );
  pManager->AddDiscreteProcess(ainelProc);

  //#########################
  //##  Generic ion        ##
  //#########################

  pManager = G4GenericIon::GenericIon()->GetProcessManager();

  // EM processes for generic ion
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
	//pManager->AddProcess(new G4hIonisation(),      -1, 2, 2);
  pManager->AddProcess(new G4ionIonisation(),      -1, 2, 2);
}
