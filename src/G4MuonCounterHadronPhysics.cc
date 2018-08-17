/**
* @file HadronPhysics.cc
* @class HadronPhysics
* @brief Build hadron particles and their physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#include "G4MuonCounterHadronPhysics.hh"
#include "G4MuonCounterSimulator.hh"

#include "G4Version.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// processes
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#if (G4VERSION_NUMBER>=930)
#include "G4hMultipleScattering.hh"
#else
#include "G4MultipleScattering.hh"
#endif

#include "G4hIonisation.hh"
#include "G4HadronElasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4HEKaonZeroInelastic.hh"

#include "G4ProtonInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"

// Stopping processes
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorptionAtRest.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

// cross sections
#include "G4PiNuclearCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4BGGNucleonInelasticXS.hh"

// models
#include "G4LElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEKaonZeroLongInelastic.hh"
#include "G4HEKaonZeroShortInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"


using namespace G4MuonCounterSimulatorUSC;


G4MuonCounterHadronPhysics::G4MuonCounterHadronPhysics(const G4String& name)
                          :G4VPhysicsConstructor(name)
{;}

G4MuonCounterHadronPhysics::~G4MuonCounterHadronPhysics()
{;}


void G4MuonCounterHadronPhysics::ConstructParticle()
{
  //  Construct all mesons
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  //  Construct all baryons
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}



void G4MuonCounterHadronPhysics::ConstructProcess()
{
  //## Hadronic Elastic Process and Model (the same for all hadrons)

  G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();
  G4LElastic* elasticModel = new G4LElastic();
  elasticProcess->RegisterMe(elasticModel);

  G4ProcessManager * pManager = 0;

  //####################
  //##  pi+ physics   ##
  //####################
  pManager = G4PionPlus::PionPlus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
		pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif

  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
 
  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4PionPlusInelasticProcess* pionPlusInelProc = new G4PionPlusInelasticProcess();
  G4PiNuclearCrossSection* pion_XC = new G4PiNuclearCrossSection();// check!
  pionPlusInelProc->AddDataSet(pion_XC);//check!

  G4LEPionPlusInelastic* pionPlusLEInelProcModel= new G4LEPionPlusInelastic();
  pionPlusInelProc->RegisterMe(pionPlusLEInelProcModel);
  G4HEPionPlusInelastic* pionPlusHEInelProcModel= new G4HEPionPlusInelastic();
  pionPlusInelProc->RegisterMe(pionPlusHEInelProcModel);

  pManager->AddDiscreteProcess(pionPlusInelProc);

  //####################
  //##  pi- physics   ##
  //####################
  pManager = G4PionMinus::PionMinus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);  
	#endif

	pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4PionMinusInelasticProcess* pionMinusInelProc= new G4PionMinusInelasticProcess();
  pionMinusInelProc->AddDataSet(pion_XC);//check!

  G4LEPionMinusInelastic* pionMinusLEInelProcModel= new G4LEPionMinusInelastic();
  pionMinusInelProc->RegisterMe(pionMinusLEInelProcModel);
  G4HEPionMinusInelastic* pionMinusHEInelProcModel= new G4HEPionMinusInelastic();
  pionMinusInelProc->RegisterMe(pionMinusHEInelProcModel);

  pManager->AddDiscreteProcess(pionMinusInelProc);

  // pi- absorption at rest
  //G4PionMinusAbsorptionAtRest* pionMinusAbsorbProc = new G4PionMinusAbsorptionAtRest();//check!
  //pManager->AddRestProcess(pionMinusAbsorbProc);//check!
	pManager->AddRestProcess(new G4PiMinusAbsorptionAtRest, ordDefault);
   
  //####################
  //##  K+ physics   ##
  //####################
  pManager = G4KaonPlus::KaonPlus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4KaonPlusInelasticProcess* kaonPlusInelProc= new G4KaonPlusInelasticProcess();
  G4LEKaonPlusInelastic* kaonPlusLEInelProcModel= new G4LEKaonPlusInelastic();
  kaonPlusInelProc->RegisterMe(kaonPlusLEInelProcModel);
  G4HEKaonPlusInelastic* kaonPlusHEInelProcModel= new G4HEKaonPlusInelastic();
  kaonPlusInelProc->RegisterMe(kaonPlusHEInelProcModel);
  pManager->AddDiscreteProcess(kaonPlusInelProc);

  //####################
  //##  K- physics   ##
  //####################
  pManager = G4KaonMinus::KaonMinus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4KaonMinusInelasticProcess* kaonMinusInelProc = new G4KaonMinusInelasticProcess();

  G4LEKaonMinusInelastic* kaonMinusLEInelProcModel = new G4LEKaonMinusInelastic();
  kaonMinusInelProc->RegisterMe(kaonMinusLEInelProcModel);
  G4HEKaonMinusInelastic* kaonMinusHEInelProcModel = new G4HEKaonMinusInelastic();
  kaonMinusInelProc->RegisterMe(kaonMinusHEInelProcModel);

  pManager->AddDiscreteProcess(kaonMinusInelProc);

  // K- absorption at rest
  //G4KaonMinusAbsorption* kaonMinusAbsorbProc = new G4KaonMinusAbsorption();//check!
  //pManager->AddRestProcess(kaonMinusAbsorbProc);//check!
	pManager->AddRestProcess(new G4KaonMinusAbsorptionAtRest, ordDefault);

  
	//####################
  //##  K0L physics   ##
  //####################
  pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4KaonZeroLInelasticProcess* kaon0LongInelProc = new G4KaonZeroLInelasticProcess();

	G4LEKaonZeroLInelastic* kaon0LongLEInelProcModel = new G4LEKaonZeroLInelastic;
	kaon0LongInelProc->RegisterMe(kaon0LongLEInelProcModel);

	G4HEKaonZeroInelastic* kaon0LongHEInelProcModel = new G4HEKaonZeroInelastic;
	kaon0LongInelProc->RegisterMe(kaon0LongHEInelProcModel);

	//G4LEKaonZeroInelastic* LEPk0LModel = new G4LEKaonZeroInelastic();//check!
  //kaon0LInelProc->RegisterMe(LEPk0LModel);//check!
  //G4HEKaonZeroLongInelastic* HEPk0LModel = new G4HEKaonZeroLongInelastic();//check!
  //kaon0LInelProc->RegisterMe(HEPk0LModel);//check!

  pManager->AddDiscreteProcess(kaon0LongInelProc);


  //####################
  //##  K0S physics   ##
  //####################
  pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4KaonZeroSInelasticProcess* kaon0ShortInelProc = new G4KaonZeroSInelasticProcess();

	//G4LEKaonZeroInelastic* kaon0ShortLEInelProcModel = new G4LEKaonZeroInelastic();//check!
  //kaon0ShortInelProc->RegisterMe(kaon0ShortLEInelProcModel);//check!
  //G4HEKaonZeroShortInelastic* kaon0ShortHEInelProcModel = new G4HEKaonZeroShortInelastic();//check!
  //kaon0ShortInelProc->RegisterMe(kaon0ShortHEInelProcModel);//check!

	G4LEKaonZeroSInelastic* kaon0ShortLEInelProcModel = new G4LEKaonZeroSInelastic;
	kaon0ShortInelProc->RegisterMe(kaon0ShortLEInelProcModel);
	G4HEKaonZeroInelastic* kaon0ShortHEInelProcModel = new G4HEKaonZeroInelastic;
	kaon0ShortInelProc->RegisterMe(kaon0ShortHEInelProcModel);

  pManager->AddDiscreteProcess(kaon0ShortInelProc );

  //#######################
  //##  proton physics   ##
  //#######################
  pManager = G4Proton::Proton()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4ProtonInelasticProcess* protonInelProc= new G4ProtonInelasticProcess();

  G4LEProtonInelastic* protonLEInelProcModel = new G4LEProtonInelastic();
  protonInelProc->RegisterMe(protonLEInelProcModel);
  G4HEProtonInelastic* protonHEInelProcModel = new G4HEProtonInelastic();
  protonInelProc->RegisterMe(protonHEInelProcModel);

  pManager->AddDiscreteProcess(protonInelProc);


  //############################
  //##  anti-proton physics   ##
  //############################
  pManager = G4AntiProton::AntiProton()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiProtonInelasticProcess* antiprotonInelProc = new G4AntiProtonInelasticProcess();

  G4LEAntiProtonInelastic* antiprotonLEInelProcModel = new G4LEAntiProtonInelastic(); 
  antiprotonInelProc->RegisterMe(antiprotonLEInelProcModel);
  G4HEAntiProtonInelastic* antiprotonHEInelProcModel = new G4HEAntiProtonInelastic(); 
  antiprotonInelProc->RegisterMe(antiprotonHEInelProcModel);
  pManager->AddDiscreteProcess(antiprotonInelProc);

  // anti-proton annihilation at rest
  G4AntiProtonAnnihilationAtRest* antiprotonAnnihilation = new G4AntiProtonAnnihilationAtRest();
  pManager->AddRestProcess(antiprotonAnnihilation);
  


  //#################################
  //##  anti-neutron physics       ##
  //#################################
  pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiNeutronInelasticProcess* antineutronInelProc = new G4AntiNeutronInelasticProcess();

  G4LEAntiNeutronInelastic* antineutronLEInelProcModel= new G4LEAntiNeutronInelastic(); 
  antineutronInelProc->RegisterMe(antineutronLEInelProcModel);

  G4HEAntiNeutronInelastic* antineutronHEInelProcModel= new G4HEAntiNeutronInelastic(); 
  antineutronInelProc->RegisterMe(antineutronHEInelProcModel);

  pManager->AddDiscreteProcess(antineutronInelProc);

  // anti-neutron annihilation at rest
  G4AntiNeutronAnnihilationAtRest* antineutronAnnihilProc = new G4AntiNeutronAnnihilationAtRest();
  pManager->AddRestProcess(antineutronAnnihilProc);


  //#################################
  //##  Lambda physics             ##
  //#################################
  pManager = G4Lambda::Lambda()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4LambdaInelasticProcess* linelProc = new G4LambdaInelasticProcess();

  G4LELambdaInelastic* LEPlModel = new G4LELambdaInelastic(); 
  linelProc->RegisterMe(LEPlModel);
  
  G4HELambdaInelastic* HEPlModel = new G4HELambdaInelastic(); 
  linelProc->RegisterMe(HEPlModel);
  
  pManager->AddDiscreteProcess(linelProc);

  //#################################
  //##  Anti-Lambda physics        ##
  //#################################
  pManager = G4AntiLambda::AntiLambda()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiLambdaInelasticProcess* alinelProc = new G4AntiLambdaInelasticProcess();
  G4LEAntiLambdaInelastic* LEPalModel = new G4LEAntiLambdaInelastic(); 
  alinelProc->RegisterMe(LEPalModel);

  G4HEAntiLambdaInelastic* HEPalModel = new G4HEAntiLambdaInelastic(); 
  alinelProc->RegisterMe(HEPalModel);

  pManager->AddDiscreteProcess(alinelProc);

  //#################################
  //##  Sigma- physics             ##
  //#################################
  pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4SigmaMinusInelasticProcess* sminelProc = new G4SigmaMinusInelasticProcess();
  G4LESigmaMinusInelastic* LEPsmModel = new G4LESigmaMinusInelastic(); 
  sminelProc->RegisterMe(LEPsmModel);
  G4HESigmaMinusInelastic* HEPsmModel = new G4HESigmaMinusInelastic(); 
  sminelProc->RegisterMe(HEPsmModel);
  pManager->AddDiscreteProcess(sminelProc);

  //#################################
  //##  AntiSigma- physics         ##
  //#################################
  pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiSigmaMinusInelasticProcess* asminelProc = new G4AntiSigmaMinusInelasticProcess();
  G4LEAntiSigmaMinusInelastic* LEPasmModel = new G4LEAntiSigmaMinusInelastic(); 
  asminelProc->RegisterMe(LEPasmModel);
  G4HEAntiSigmaMinusInelastic* HEPasmModel = new G4HEAntiSigmaMinusInelastic(); 
  asminelProc->RegisterMe(HEPasmModel);
  pManager->AddDiscreteProcess(asminelProc);

  //#################################
  //##  Sigma+ physics             ##
  //#################################
  pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4SigmaPlusInelasticProcess* spinelProc = new G4SigmaPlusInelasticProcess();
  G4LESigmaPlusInelastic* LEPspModel = new G4LESigmaPlusInelastic(); 
  spinelProc->RegisterMe(LEPspModel);
  G4HESigmaPlusInelastic* HEPspModel = new G4HESigmaPlusInelastic(); 
  spinelProc->RegisterMe(HEPspModel);
  pManager->AddDiscreteProcess(spinelProc);

  //#################################
  //##  AntiSigma+ physics         ##
  //#################################
  pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiSigmaPlusInelasticProcess* aspinelProc = new G4AntiSigmaPlusInelasticProcess();
  G4LEAntiSigmaPlusInelastic* LEPaspModel = new G4LEAntiSigmaPlusInelastic(); 
  aspinelProc->RegisterMe(LEPaspModel);
  G4HEAntiSigmaPlusInelastic* HEPaspModel = new G4HEAntiSigmaPlusInelastic(); 
  aspinelProc->RegisterMe(HEPaspModel);
  pManager->AddDiscreteProcess(aspinelProc);

  //#################################
  //##  Xi- physics                ##
  //#################################
  pManager = G4XiMinus::XiMinus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4XiMinusInelasticProcess* xminelProc = new G4XiMinusInelasticProcess();
  G4LEXiMinusInelastic* LEPxmModel = new G4LEXiMinusInelastic(); 
  xminelProc->RegisterMe(LEPxmModel);
  G4HEXiMinusInelastic* HEPxmModel = new G4HEXiMinusInelastic(); 
  xminelProc->RegisterMe(HEPxmModel);
  pManager->AddDiscreteProcess(xminelProc);

  //#################################
  //##  AntiXi- physics            ##
  //#################################
  pManager = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);	
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);		
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiXiMinusInelasticProcess* axminelProc = new G4AntiXiMinusInelasticProcess();
  G4LEAntiXiMinusInelastic* LEPaxmModel = new G4LEAntiXiMinusInelastic(); 
  axminelProc->RegisterMe(LEPaxmModel);
  G4HEAntiXiMinusInelastic* HEPaxmModel = new G4HEAntiXiMinusInelastic(); 
  axminelProc->RegisterMe(HEPaxmModel);
  pManager->AddDiscreteProcess(axminelProc);

  //#################################
  //##  Xi0 physics                ##
  //#################################
  pManager = G4XiZero::XiZero()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4XiZeroInelasticProcess* x0inelProc = new G4XiZeroInelasticProcess();
  G4LEXiZeroInelastic* LEPx0Model = new G4LEXiZeroInelastic(); 
  x0inelProc->RegisterMe(LEPx0Model);
  G4HEXiZeroInelastic* HEPx0Model = new G4HEXiZeroInelastic(); 
  x0inelProc->RegisterMe(HEPx0Model);
  pManager->AddDiscreteProcess(x0inelProc);

  //#################################
  //##  AntiXi0 physics            ##
  //#################################
  pManager = G4AntiXiZero::AntiXiZero()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiXiZeroInelasticProcess* ax0inelProc = new G4AntiXiZeroInelasticProcess();
  G4LEAntiXiZeroInelastic* LEPax0Model = new G4LEAntiXiZeroInelastic(); 
  ax0inelProc->RegisterMe(LEPax0Model);
  G4HEAntiXiZeroInelastic* HEPax0Model = new G4HEAntiXiZeroInelastic(); 
  ax0inelProc->RegisterMe(HEPax0Model);
  pManager->AddDiscreteProcess(ax0inelProc);

  //#################################
  //##  Omega- physics             ##
  //#################################
  pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4OmegaMinusInelasticProcess* ominelProc = new G4OmegaMinusInelasticProcess();
  G4LEOmegaMinusInelastic* LEPomModel = new G4LEOmegaMinusInelastic(); 
  ominelProc->RegisterMe(LEPomModel);
  G4HEOmegaMinusInelastic* HEPomModel = new G4HEOmegaMinusInelastic(); 
  ominelProc->RegisterMe(HEPomModel);
  pManager->AddDiscreteProcess(ominelProc);

  //#################################
  //##  AntiOmega- physics         ##
  //#################################
  pManager = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();

  // EM processes
	#if (G4VERSION_NUMBER>=930)
  	pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	#else
		pManager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
	#endif
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiOmegaMinusInelasticProcess* aominelProc = new G4AntiOmegaMinusInelasticProcess();
  G4LEAntiOmegaMinusInelastic* LEPaomModel = new G4LEAntiOmegaMinusInelastic(); 
  aominelProc->RegisterMe(LEPaomModel);
  G4HEAntiOmegaMinusInelastic* HEPaomModel = new G4HEAntiOmegaMinusInelastic(); 
  aominelProc->RegisterMe(HEPaomModel);
  pManager->AddDiscreteProcess(aominelProc);

}





