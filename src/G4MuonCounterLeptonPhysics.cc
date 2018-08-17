/**
* @file LeptonPhysics.cc
* @class LeptonPhysics
* @brief Build lepton particles and their physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#include "G4MuonCounterLeptonPhysics.hh"
#include "G4Version.hh"


#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// processes
#if (G4VERSION_NUMBER>=930)
	#include "G4eMultipleScattering.hh"
#else
	#include "G4MultipleScattering.hh"
#endif

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4AnnihiToMuPair.hh"
#include "G4eeToHadrons.hh"

#if (G4VERSION_NUMBER>=930)
	#include "G4MuMultipleScattering.hh"
  #include "G4WentzelVIModel.hh"
  #include "G4CoulombScattering.hh"
#endif
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

#if (G4VERSION_NUMBER>=930)
	#include "G4hMultipleScattering.hh"
#endif

#include "G4hIonisation.hh"

#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"

// models
#include "G4ElectroNuclearReaction.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"

#include "G4EmProcessOptions.hh"

#include "G4MuonCounterSimulator.hh"

using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterLeptonPhysics::G4MuonCounterLeptonPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{;}


G4MuonCounterLeptonPhysics::~G4MuonCounterLeptonPhysics()
{;}


void G4MuonCounterLeptonPhysics::ConstructParticle()
{ 
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();
}


void G4MuonCounterLeptonPhysics::ConstructProcess()
{
  G4ProcessManager* pManager = 0;

	//## Options for nuclear reactions
  // 1) Model for e+/e- nuclear reactions
  G4ElectroNuclearReaction* theElectronReaction = new G4ElectroNuclearReaction();
	// 2) Photonuclear cross sections
	G4PhotoNuclearCrossSection* photoCrossSection = new G4PhotoNuclearCrossSection();

	//#########################
  //##  Electron physics   ##
  //#########################
  pManager = G4Electron::Electron()->GetProcessManager();
	#if (G4VERSION_NUMBER>=930)
		pManager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4eIonisation(),        -1, 2, 2);
  //pManager->AddProcess(new G4eBremsstrahlung(),    -1, 3, 3); //OLD
	pManager->AddProcess(new G4eBremsstrahlung(),    -1, -3, 3);  

  G4ElectronNuclearProcess* theElectronNuclearProcess = new G4ElectronNuclearProcess();	
  theElectronNuclearProcess->RegisterMe(theElectronReaction);//model
	//theElectronNuclearProcess->AddDataSet(photoCrossSection);//cross sections
	
  pManager->AddDiscreteProcess(theElectronNuclearProcess);//using G4PhotoNuclearCrossSection as default

  //#########################
  //##  Positron physics   ##
  //#########################
  pManager = G4Positron::Positron()->GetProcessManager();
	#if (G4VERSION_NUMBER>=930) 
		pManager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4eIonisation(),        -1, 2, 2);
  //pManager->AddProcess(new G4eBremsstrahlung(),    -1, 3, 3);//OLD 
	pManager->AddProcess(new G4eBremsstrahlung(),    -1, -3, 3);  
  pManager->AddProcess(new G4eplusAnnihilation(),   0,-1, 4);

  G4PositronNuclearProcess* thePositronNuclearProcess = new G4PositronNuclearProcess();
  thePositronNuclearProcess->RegisterMe(theElectronReaction);//model
	//thePositronNuclearProcess->AddDataSet(photoCrossSection);//cross sections
	
  pManager->AddDiscreteProcess(thePositronNuclearProcess);//using G4PhotoNuclearCrossSection as default

	pManager->AddDiscreteProcess(new G4AnnihiToMuPair());//added with respect to standard EM physics
  pManager->AddDiscreteProcess(new G4eeToHadrons());//added with respect to standard EM physics


  //#########################
  //##  Muon- physics      ##
  //#########################
  pManager = G4MuonMinus::MuonMinus()->GetProcessManager(); 
	#if (G4VERSION_NUMBER>=930) 
		G4MuMultipleScattering* msc = new G4MuMultipleScattering();
    msc->AddEmModel(0, new G4WentzelVIModel());
    pManager->AddProcess(msc,-1, 1, 1);

		pManager->AddDiscreteProcess(new G4CoulombScattering());
	#else
		pManager->AddProcess(new G4MultipleScattering(),    -1, 1, 1);
	#endif
  pManager->AddProcess(new G4MuIonisation(),          -1, 2, 2);
  //pManager->AddProcess(new G4MuBremsstrahlung(),      -1, 3, 3); //OLD
	pManager->AddProcess(new G4MuBremsstrahlung(),      -1, -3, 3);  
  //pManager->AddProcess(new G4MuPairProduction(),      -1, 4, 4);//OLD
	pManager->AddProcess(new G4MuPairProduction(),      -1, -4, 4);
	
	G4MuonMinusCaptureAtRest* fMuMinusCaptureAtRest = new G4MuonMinusCaptureAtRest();
	pManager->AddRestProcess(fMuMinusCaptureAtRest);
	

  //#########################
  //##  Muon+ physics      ##
  //#########################
  pManager = G4MuonPlus::MuonPlus()->GetProcessManager(); 
	#if (G4VERSION_NUMBER>=930) 
    pManager->AddProcess(msc,-1, 1, 1);

		pManager->AddDiscreteProcess(new G4CoulombScattering());
	#else
		pManager->AddProcess(new G4MultipleScattering(),    -1, 1, 1);
	#endif
  pManager->AddProcess(new G4MuIonisation(),          -1, 2, 2);
  //pManager->AddProcess(new G4MuBremsstrahlung(),      -1, 3, 3); //OLD
	pManager->AddProcess(new G4MuBremsstrahlung(),      -1, -3, 3);  
  //pManager->AddProcess(new G4MuPairProduction(),      -1, 4, 4);//OLD
	pManager->AddProcess(new G4MuPairProduction(),      -1, -4, 4);


	//#########################
  //##  Tau- physics      ##
  //#########################
  pManager = G4TauMinus::TauMinus()->GetProcessManager();
	#if (G4VERSION_NUMBER>=930) 
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif

  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
 
  //#########################
  //##  Tau+ physics      ##
  //#########################
	pManager = G4TauPlus::TauPlus()->GetProcessManager();
	
  #if (G4VERSION_NUMBER>=930) 
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

	#if (G4VERSION_NUMBER>=930) 
		G4EmProcessOptions opt;
  	//opt.SetVerbose(verbose);
 	 	opt.SetPolarAngleLimit(0.2);
	#endif
}
