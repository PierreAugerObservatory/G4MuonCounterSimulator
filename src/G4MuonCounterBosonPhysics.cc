/**
* @file G4MuonCounterBosonPhysics.cc
* @class G4MuonCounterBosonPhysics
* @brief Build boson particles and their physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterBosonPhysics.hh"
#include "G4ProcessManager.hh"

// processes
#include "G4GammaConversion.hh"
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4PhotoNuclearProcess.hh"
#include "G4GammaConversionToMuons.hh"

//## optical
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"
#include "G4Decay.hh"

// models
#include "G4PhotoNuclearCrossSection.hh"
#include "G4GammaNuclearReaction.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSModel.hh"

#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"

using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterBosonPhysics::G4MuonCounterBosonPhysics(const G4String& name)
:  G4VPhysicsConstructor(name)
{;}


G4MuonCounterBosonPhysics::~G4MuonCounterBosonPhysics()
{;}


void G4MuonCounterBosonPhysics::ConstructParticle()
{
	// pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();
}


void G4MuonCounterBosonPhysics::ConstructProcess()
{
  //
  // First define the gamma-nuclear models
  //

  // low energy part
  G4GammaNuclearReaction* lowEGammaModel = new G4GammaNuclearReaction();
  lowEGammaModel->SetMaxEnergy(3.5*GeV);
 
  // high energy part
  G4TheoFSGenerator* highEGammaModel = new G4TheoFSGenerator();
  G4GeneratorPrecompoundInterface* preComModel = new G4GeneratorPrecompoundInterface();
  highEGammaModel->SetTransport(preComModel);
 
  G4QGSModel<G4GammaParticipants>* theStringModel = new G4QGSModel<G4GammaParticipants>;
  G4QGSMFragmentation* fragModel = new G4QGSMFragmentation();
  G4ExcitedStringDecay* stringDecay = new G4ExcitedStringDecay(fragModel);
  theStringModel->SetFragmentationModel(stringDecay);
 
  highEGammaModel->SetHighEnergyGenerator(theStringModel);
  highEGammaModel->SetMinEnergy(3.*GeV);
  highEGammaModel->SetMaxEnergy(100.*TeV); 

	// data cross sections
	G4PhotoNuclearCrossSection* photoCrossSection = new G4PhotoNuclearCrossSection();


  // Now add the processes to the gamma, including e+e- pair creation, 
  // Compton scattering and photo-electric effect

  G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();
  pManager->AddDiscreteProcess(new G4GammaConversion());
  pManager->AddDiscreteProcess(new G4ComptonScattering());
  pManager->AddDiscreteProcess(new G4PhotoElectricEffect());
	pManager->AddDiscreteProcess(new G4GammaConversionToMuons());	

  G4PhotoNuclearProcess* thePhotoNuclearProcess = new G4PhotoNuclearProcess();
  thePhotoNuclearProcess->RegisterMe(lowEGammaModel);//model LE
  thePhotoNuclearProcess->RegisterMe(highEGammaModel);//model HE
	//thePhotoNuclearProcess->AddDataSet(photoCrossSection);//cross section

  pManager->AddDiscreteProcess(thePhotoNuclearProcess);//using G4PhotoNuclearCrossSection as default


}



