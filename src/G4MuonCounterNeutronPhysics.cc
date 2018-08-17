/**
* @file G4MuonCounterNeutronPhysics.cc
* @class G4MuonCounterNeutronPhysics
* @brief Build neutron particles and their physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterNeutronPhysics.hh"

#include "G4ParticleTable.hh"

// processes
#include "G4ProcessManager.hh"
#include "G4HadronElasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

// cross sections

// models
#include "G4LElastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

#include "G4MuonCounterSimulator.hh"

using namespace G4MuonCounterSimulatorUSC;


G4MuonCounterNeutronPhysics::G4MuonCounterNeutronPhysics(const G4String& name)
                          :G4VPhysicsConstructor(name)
{;}

G4MuonCounterNeutronPhysics::~G4MuonCounterNeutronPhysics()
{;}


void G4MuonCounterNeutronPhysics::ConstructParticle()
{
  //  No need to construct anything here as long as 
  //  HadronPhysics is constructed
}



void G4MuonCounterNeutronPhysics::ConstructProcess()
{

  G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();

  //## Neutron elastic process, models and cross sections
  G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();

  G4LElastic* theLElasticModel = new G4LElastic();
	//theLElasticModel->SetMinEnergy(19*MeV);//check!
  elasticProcess->RegisterMe(theLElasticModel);

/*
	G4NeutronHPElastic* neutronHPElProcModel= new G4NeutronHPElastic;//check!
	elasticProcess->RegisterMe(neutronHPElProcModel);//check!
	G4NeutronHPElasticData* neutronHPElProcData = new G4NeutronHPElasticData;//check!
	elasticProcess->AddDataSet(neutronHPElProcData);//check!
*/

  pManager->AddDiscreteProcess(elasticProcess);


  //## Neutron inelastic process, models and cross sections
  //## Use Quark-Gluon String Model between 15 GeV and 100 TeV
  G4NeutronInelasticProcess* ninelProc = new G4NeutronInelasticProcess();

  // Use LEP model between 9.5 and 25 GeV
  G4LENeutronInelastic* LEPnModel = new G4LENeutronInelastic();
	//LEPnModel->SetMinEnergy(19*MeV);//check!
  ninelProc->RegisterMe(LEPnModel);

/*
	G4NeutronHPInelastic* neutronHPInelProcModel = new G4NeutronHPInelastic;//check!
	ninelProc->RegisterMe(neutronHPInelProcModel);//check!
	G4NeutronHPInelasticData* neutronHPInelProcData = new G4NeutronHPInelasticData;//check!
	ninelProc->AddDataSet(neutronHPInelProcData);//check!
*/

  //## Neutron-induced fission process, models and cross sections
  G4HadronFissionProcess* neutronFission = new G4HadronFissionProcess();
  G4LFission* theLFissionModel = new G4LFission();
  theLFissionModel->SetMaxEnergy(20.*TeV);
  neutronFission->RegisterMe(theLFissionModel);
  pManager->AddDiscreteProcess(neutronFission);

  //## Neutron capture process, models and cross sections
  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
  G4LCapture* theLCaptureModel = new G4LCapture();
	//theLCaptureModel->SetMinEnergy(19*MeV);
  theLCaptureModel->SetMaxEnergy(20.*TeV);
  neutronCapture->RegisterMe(theLCaptureModel);

/*
	G4NeutronHPCapture* neutronHPCaptureProcModel= new G4NeutronHPCapture;//check!
	neutronCapture->RegisterMe(neutronHPCaptureProcModel);//check!
	G4NeutronHPCaptureData* neutronHPCaptureProcData = new G4NeutronHPCaptureData;//check!
	neutronCapture->AddDataSet(neutronHPCaptureProcData);//check!
*/	
  pManager->AddDiscreteProcess(neutronCapture);


}
