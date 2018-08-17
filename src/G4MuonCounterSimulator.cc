/**
* @mainpage StripSimulation tool
* <img src="SimSampleScreeshot.jpg" alt="Screenshot">
* <b>Description of the project </b> <br>
* This project implements a detailed GEANT4 simulation (including generation and tracking of optical photons) of particles injected in planes of scintillator modules. <br>
A strip module is made up of a scintillator bar, surrounded by a reflective coating. A groove is done in the bar to host a WLS fiber. Two PMTs are placed at the end of the scintillator bar. <br>
The detector plane module is made up of a series of strip modules, one aside to the other. The detector geometry can be customized by the user (i.e. by setting the number of the number of planes, number of strips in a plane, ...).
The simulation output is recorded in ROOT format.
* 
* @author S. Riggi, E. Trovato <br> Dept. of Physics and Astronomy, University of Catania <br> INFN, Section of Catania  <br> Italy
*/

/**
* @file G4MuonCounterSimulator.cc
* @class module
* @brief main class used to handle the entire project
*
* @author S. Riggi
* @date 05/04/2010
*/

#include <G4MuonCounterSimulator.hh>
#include <G4MuonCounterPhysicsList.hh>
#include <G4MuonCounterConstruction.hh>
#include <G4MuonCounterPrimaryGenerator.hh>
#include <G4MuonCounterTrackingAction.hh>
#include <G4MuonCounterStackingAction.hh>
#include <G4MuonCounterSteppingAction.hh>
#include <G4MuonCounterEventAction.hh>
#include <G4MuonCounterRunAction.hh>
#include <G4MuonCounterRecorderBase.hh>
#include <G4MuonCounterSteppingVerbose.hh>
//#include <G4MuonCounterVisManager.hh>

#include <G4MuonCounterSteppingVerbose.hh>
#include <G4MuonCounterMCReader.hh>
#include <Utilities.hh>

#include <TStationSimData.hh>
#include <MuonDetector.hh>

#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4VisManager.hh>

#include <G4UIterminal.hh>
#include <G4UItcsh.hh>

#include <G4ios.hh>
#include <Randomize.hh>

//#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
//#endif

//## offline headers
//#include "MyDetector.h"
//#include "Detector.h"
//#include <MDetector.h>
//#include <MStation.h>

#include <fwk/CentralConfig.h>
#include <fwk/RandomEngineRegistry.h>

#include <det/Detector.h>

#include <evt/Event.h>
#include <evt/ShowerSimData.h>

#include <sdet/SDetector.h>
#include <sdet/Station.h>

#include <sevt/PMT.h>
#include <sevt/PMTSimData.h>
#include <sevt/SEvent.h>
#include <sevt/Station.h>
#include <sevt/StationSimData.h>

#include <utl/ErrorLogger.h>
#include <utl/Reader.h>
#include <utl/Particle.h>
#include <utl/ShowerParticleIterator.h>
#include <utl/TimeDistribution.h>
#include <utl/TimeDistributionAlgorithm.h>

#include <CLHEP/Random/Random.h>

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TSystem.h>
#include <Rtypes.h>
#include <TObject.h>
#include <TRint.h>
#include <TMath.h>
#include "TScintHit.hh"
#include "TPMTHit.hh"

#include <vector>
#include <cstddef>
#include <iostream>
#include <sstream>


using namespace utl;
using namespace fwk;
using namespace std;
using namespace sevt;
using namespace sdet;
using namespace evt;
using namespace det;
//using namespace mdet;
using namespace G4MuonCounterSimulatorUSC;

const sdet::Station* G4MuonCounterSimulator::fCurrentDetectorStation = 0;
SEvent::StationIterator G4MuonCounterSimulator::fCurrentEventStationIt;
StationSimData::ParticleIterator G4MuonCounterSimulator::fCurrentParticleIt;
std::vector<StationSimData::ParticleIterator> G4MuonCounterSimulator::fCurrentParticleItList;

const evt::Event* G4MuonCounterSimulator::fCurrentEvent= 0;

TStationSimData* G4MuonCounterSimulator::fCurrentStationSimData;
std::vector<TStationSimData> G4MuonCounterSimulator::fStationSimDataCollection;
TEventSimData* G4MuonCounterSimulator::fCurrentEventSimData;
std::vector<TEventSimData> G4MuonCounterSimulator::fEventSimDataCollection;
MuonDetector* G4MuonCounterSimulator::fMuonDetector; 
G4MuonCounterConstruction* G4MuonCounterSimulator::fMuonCounterConstruction; 

bool G4MuonCounterSimulator::fgAllProcess;
bool G4MuonCounterSimulator::fgOptical;
bool G4MuonCounterSimulator::fgMuCapture;
bool G4MuonCounterSimulator::fgScintillation;
bool G4MuonCounterSimulator::fgWLS;
bool G4MuonCounterSimulator::fgCherenkov;
bool G4MuonCounterSimulator::fgAbsorption;
bool G4MuonCounterSimulator::fgRayleigh;
bool G4MuonCounterSimulator::fgBoundary;
bool G4MuonCounterSimulator::fOneParticlePerStation;

/**
* main 
*/
G4MuonCounterSimulator::G4MuonCounterSimulator() :
  fRunManager(0),
  fUImanager(0),
  fVisManager(0),
	fRecorder(0),
  fGeoVisOn(0),
  fTrajVisOn(0),
	fVisDriverType(0),
	fRunVerbosity(0),
	fEventVerbosity(0),
	fTrackingVerbosity(0),
  fDetectorConstructed(0),
  fFastMode(false),
	fUseG4StandAloneConfiguration(false),	
	fSaveSimInfo(0),
	fSaveFullInfo(0),
  fEventId("")
{
}


G4MuonCounterSimulator::~G4MuonCounterSimulator()
{
}


VModule::ResultFlag
G4MuonCounterSimulator::Init()
{

	INFO("Geant4 muon counter sim (beta-version)");
	//##################################
	//## get some XML config
	//##################################
  Branch topB = CentralConfig::GetInstance()->GetTopBranch("G4MuonCounterSimulatorUSC");
	if (!topB) {
    ERROR("No XML configuration found for G4MuonCounterSimulatorUSC!");
    return eFailure;
  }

	map<string, string> useAtt;
  useAtt["use"] = string("yes");

	Branch saveSimInfoB = topB.GetChild("SaveSimInfo", useAtt);
  if (saveSimInfoB) {
    fSaveSimInfo = true;
    saveSimInfoB.GetChild("OutFileName").GetData(fOutputFileName);
		saveSimInfoB.GetChild("saveFullInfo").GetData(fSaveFullInfo);		
  }

	if(fSaveSimInfo) {
		cout<<"Saving simulation output in "<< fOutputFileName.c_str()<< endl;
		fOutputFile= new TFile(fOutputFileName.c_str(),"RECREATE");
		fSimInfo= new TTree("SimInfo","SimInfo");	
		fSimInfo->Branch("EventSimData","TEventSimData",&fCurrentEventSimData);
		fSimInfo->Branch("StationSimDataCollection",&fStationSimDataCollection);

		fDetInfo= new TTree("DetInfo","DetInfo");	
		fDetInfo->Branch("Detector","MuonDetector",&fMuonDetector);
	}

	//## verbosity
  Branch verbosityB = topB.GetChild("verbosity");
	verbosityB.GetChild("runVerbosity").GetData(fRunVerbosity);
  verbosityB.GetChild("eventVerbosity").GetData(fEventVerbosity);
  verbosityB.GetChild("trackingVerbosity").GetData(fTrackingVerbosity);
	cout<<"Verbosity mode: run="<<fRunVerbosity<<" , event="<<fEventVerbosity<<" , tracking="<<fTrackingVerbosity<<endl;  


	//## visualization
  Branch visB = topB.GetChild("visualization");
	visB.GetChild("visDriver").GetData(fVisDriverType);
  visB.GetChild("geometry").GetData(fGeoVisOn);
  visB.GetChild("trajectories").GetData(fTrajVisOn);
	cout<<"Visualization mode: geometry="<<fGeoVisOn<<" , trajectories="<<fTrajVisOn<<" , driver="<<fVisDriverType<<endl;  

	//## run mode
	Branch runModeB= topB.GetChild("run");	
	Branch g4macroB= runModeB.GetChild("G4MacroFile");	
	if (!g4macroB){
		INFO("No g4 macro given, running with XML config");
		fUseG4StandAloneConfiguration= false;	
	}
	else{
		g4macroB.GetData(fG4MacroFile);
		if (fG4MacroFile != "") {
			INFO("Steering simulation with Geant4 macro file.");
			fUseG4StandAloneConfiguration= true;
		}
		else {
			INFO("Empty g4 macro string, running with XML config.");
			fUseG4StandAloneConfiguration= false;
		}
	}

  Branch fastModeB = runModeB.GetChild("fastMode");
  fastModeB.GetData(fFastMode);
	Branch fullModeB = runModeB.GetChild("fullMode");
  fullModeB.GetData(fFullMode);
	Branch oneParticlePerStationModeB = runModeB.GetChild("oneParticlePerStation");
  oneParticlePerStationModeB.GetData(fOneParticlePerStation);

	
	//## physics processes
  Branch physicsProcessesB = topB.GetChild("physicsProcesses");
	physicsProcessesB.GetChild("All").GetData(fgAllProcess);
	physicsProcessesB.GetChild("Optical").GetData(fgOptical);
	physicsProcessesB.GetChild("Scintillation").GetData(fgScintillation);
	physicsProcessesB.GetChild("Cherenkov").GetData(fgCherenkov);
	physicsProcessesB.GetChild("WLS").GetData(fgWLS);
	physicsProcessesB.GetChild("Absorption").GetData(fgAbsorption);
	physicsProcessesB.GetChild("Rayleigh").GetData(fgRayleigh);
	physicsProcessesB.GetChild("Boundary").GetData(fgBoundary);
  physicsProcessesB.GetChild("muCapture").GetData(fgMuCapture);

	//#ifdef G4VIS_USE
  if ( (fGeoVisOn || fTrajVisOn) && !fVisManager){
		cout<<"Creating G4 visualization manager"<<endl;
		fVisManager = new G4VisExecutive;		
    //fVisManager = new G4MuonCounterVisManager;
	}
  //#endif	


	G4VSteppingVerbose::SetInstance(new G4MuonCounterSteppingVerbose);

	if (!fRunManager) fRunManager = new G4RunManager;

 
	//## init muon counter detector
	G4MuonCounterConstruction* muonCounterDetector = new G4MuonCounterConstruction();
	if(fUseG4StandAloneConfiguration) muonCounterDetector->SetDetectorParametersFromXML(false);
	else muonCounterDetector->SetDetectorParametersFromXML(true);
  fRunManager->SetUserInitialization(muonCounterDetector);
	fDetectorConstructed = true;

	//## init muon counter physics
	G4MuonCounterPhysicsList* physicsList= new G4MuonCounterPhysicsList();
  fRunManager->SetUserInitialization(physicsList);	
	
	//## init recorder
	fRecorder = NULL;//No recording is done

	//## init stacking action
	fStackingAction = new G4MuonCounterStackingAction(fRecorder);
  fRunManager->SetUserAction(fStackingAction);

	//## init primary generator
	G4MuonCounterPrimaryGenerator* generator = new G4MuonCounterPrimaryGenerator(muonCounterDetector);
  fRunManager->SetUserAction(generator);
  
	//## init tracking action
  fRunManager->SetUserAction(new G4MuonCounterTrackingAction(fRecorder,generator));
	//## init stepping action
  fRunManager->SetUserAction(new G4MuonCounterSteppingAction(fRecorder));
  
	//## init run action
	G4MuonCounterRunAction* run= new G4MuonCounterRunAction(fRecorder,muonCounterDetector,generator);
	fRunManager->SetUserAction(run);
  
	//## init event action
  fRunManager->SetUserAction(new G4MuonCounterEventAction(fRecorder,run));
 	
	//## Initialize G4 kernel
	fRunManager->Initialize();

	//## Get pointer to the UI Manager
	fUImanager = G4UImanager::GetUIpointer();
  
	if(!fUseG4StandAloneConfiguration){
		//set verbosity of run
  	ostringstream runVerbosityCommand;
  	runVerbosityCommand << "/run/verbose " << fRunVerbosity; 
		ostringstream eventVerbosityCommand;
  	eventVerbosityCommand << "/event/verbose " << fEventVerbosity; 
		ostringstream trackingVerbosityCommand;
  	trackingVerbosityCommand << "/tracking/verbose " << fTrackingVerbosity;            
		fUImanager->ApplyCommand(runVerbosityCommand.str().c_str());
  	fUImanager->ApplyCommand(eventVerbosityCommand.str().c_str());
  	fUImanager->ApplyCommand(trackingVerbosityCommand.str().c_str());

		//dump all processes
		fUImanager->ApplyCommand("/process/list");


		//deselect all process
		if(!fgAllProcess){
			fUImanager->ApplyCommand("/process/inactivate Decay");
			fUImanager->ApplyCommand("/process/inactivate RadioactiveDecay");
			fUImanager->ApplyCommand("/process/inactivate conv");
			fUImanager->ApplyCommand("/process/inactivate compt");
			fUImanager->ApplyCommand("/process/inactivate phot");
			fUImanager->ApplyCommand("/process/inactivate GammaToMuPair");
			fUImanager->ApplyCommand("/process/inactivate PhotonInelastic");
			fUImanager->ApplyCommand("/process/inactivate msc");
			fUImanager->ApplyCommand("/process/inactivate eIoni");
			fUImanager->ApplyCommand("/process/inactivate eBrem");
			fUImanager->ApplyCommand("/process/inactivate ElectroNuclear");
			fUImanager->ApplyCommand("/process/inactivate annihil");
			fUImanager->ApplyCommand("/process/inactivate PositronNuclear");
			fUImanager->ApplyCommand("/process/inactivate AnnihiToMuPair");
			fUImanager->ApplyCommand("/process/inactivate ee2hadr");
			fUImanager->ApplyCommand("/process/inactivate muMsc");
			fUImanager->ApplyCommand("/process/inactivate CoulombScat");
			fUImanager->ApplyCommand("/process/inactivate muIoni");
			fUImanager->ApplyCommand("/process/inactivate muBrems");
			fUImanager->ApplyCommand("/process/inactivate muPairProd");
			fUImanager->ApplyCommand("/process/inactivate muMinusCaptureAtRest");
			fUImanager->ApplyCommand("/process/inactivate hIoni");
			fUImanager->ApplyCommand("/process/inactivate HadronElastic");
			fUImanager->ApplyCommand("/process/inactivate nFission");
			fUImanager->ApplyCommand("/process/inactivate nCapture");
			fUImanager->ApplyCommand("/process/inactivate PionPlusInelastic");
			fUImanager->ApplyCommand("/process/inactivate PionMinusInelastic");
			fUImanager->ApplyCommand("/process/inactivate PiMinusAbsorptionAtRest");
			fUImanager->ApplyCommand("/process/inactivate KaonPlusInelastic");
			fUImanager->ApplyCommand("/process/inactivate KaonMinusInelastic");
			fUImanager->ApplyCommand("/process/inactivate KaonMinusAbsorptionAtRest");
			fUImanager->ApplyCommand("/process/inactivate KaonZeroLInelastic");
			fUImanager->ApplyCommand("/process/inactivate KaonZeroSInelasti");
			fUImanager->ApplyCommand("/process/inactivate ProtonInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiProtonInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiProtonAnnihilationAtRest");
			fUImanager->ApplyCommand("/process/inactivate AntiNeutronInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiNeutronAnnihilationAtRest");
			fUImanager->ApplyCommand("/process/inactivate LambdaInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiLambdaInelastic");
			fUImanager->ApplyCommand("/process/inactivate SigmaMinusInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiSigmaMinusInelastic");
			fUImanager->ApplyCommand("/process/inactivate SigmaPlusInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiSigmaPlusInelastic");
			fUImanager->ApplyCommand("/process/inactivate XiMinusInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiXiMinusInelastic");
			fUImanager->ApplyCommand("/process/inactivate XiZeroInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiXiZeroInelasti");
			fUImanager->ApplyCommand("/process/inactivate OmegaMinusInelastic");
			fUImanager->ApplyCommand("/process/inactivate AntiOmegaMinusInelastic");
			fUImanager->ApplyCommand("/process/inactivate DeuteronInelastic");
			fUImanager->ApplyCommand("/process/inactivate TritonInelastic");
			fUImanager->ApplyCommand("/process/inactivate ionIoni");
		}//close deselect all process

	

		//deselect optical physics processes
		if(!fgScintillation && fgOptical) 
			fUImanager->ApplyCommand("/process/inactivate Scintillation");
		if(!fgCherenkov && fgOptical) 
			fUImanager->ApplyCommand("/process/inactivate Cerenkov");
		if(!fgWLS && fgOptical) 
			fUImanager->ApplyCommand("/process/inactivate OpWLS");
		if(!fgAbsorption && fgOptical) 
			fUImanager->ApplyCommand("/process/inactivate OpAbsorption");
		if(!fgRayleigh && fgOptical) 
			fUImanager->ApplyCommand("/process/inactivate OpRayleigh");
		if(!fgBoundary && fgOptical) 
			fUImanager->ApplyCommand("/process/inactivate OpBoundary");


		//set visualization 
		//#ifdef G4VIS_USE
		if (fGeoVisOn || fTrajVisOn) {
								
      fVisManager->Initialize();
		
			fUImanager->ApplyCommand("/vis/scene/create");
			if(fVisDriverType == 1){
      	fUImanager->ApplyCommand("/vis/open OGLIX"); //immediate OpenGL no Motif control
			}
			else if(fVisDriverType == 2){
				fUImanager->ApplyCommand("/vis/open OGLIXm"); //immediate OpenGL Motif control
			}
			else if(fVisDriverType == 3){
				fUImanager->ApplyCommand("/vis/open OGLSX"); //stored OpenGL Motif control
			}					
			else if(fVisDriverType == 4){
				fUImanager->ApplyCommand("/vis/open OGLSXm"); //stored OpenGL Motif control
			}
			else if(fVisDriverType == 5){
				fUImanager->ApplyCommand("/vis/open DAWNFILE"); //DAWN
				//fUImanager->ApplyCommand("/vis/scene/create");
			}
			else if(fVisDriverType == 6){
				fUImanager->ApplyCommand("/vis/open VRML2FILE");
				//fUImanager->ApplyCommand("/vis/scene/create");
			}
			else{
				cerr<<"Invalid visualization driver given in XML config...exit!"<<endl;	
				exit(1);
			}
   	
      fUImanager->ApplyCommand("/vis/sceneHandler/attach");
     	fUImanager->ApplyCommand("/vis/scene/add/volume");
     	fUImanager->ApplyCommand("/vis/viewer/set/style/wireframe");
			//fUImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 90 90 deg");//0-90 (X-Y plane), 90-90 (X-Z plane), 90-180 (Y-Z plane)
			fUImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 30 60 deg");//0-90 (X-Y plane), 90-90 (X-Z plane), 90-180 (Y-Z plane)

			fUImanager->ApplyCommand("/vis/viewer/set/background white");
			fUImanager->ApplyCommand("/vis/viewer/zoom 1");
			//fUImanager->ApplyCommand("/vis/viewer/flush");
      fUImanager->ApplyCommand("/vis/drawVolume");
     	fUImanager->ApplyCommand("/vis/scene/notifyHandlers");
			fUImanager->ApplyCommand("/vis/viewer/update");				
    }
    if (fTrajVisOn) {
    	fUImanager->ApplyCommand("/tracking/storeTrajectory 1");
      fUImanager->ApplyCommand("/vis/scene/add/trajectories"); 
			fUImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByCharge"); 
			fUImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true"); 
			fUImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2");
			fUImanager->ApplyCommand("/vis/scene/add/hits");
			fUImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2");	     	
    }
	
		fUImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
		fUImanager->ApplyCommand("/vis/viewer/refresh");

		//#endif
	}//if not stand-alone G4 config


	//initialize station sim data container
	fCurrentStationSimData= new TStationSimData;
	fStationSimDataCollection.clear();
	fStationSimDataCollection.resize(0);

	fCurrentEventSimData= new TEventSimData;
	fEventSimDataCollection.clear();
	fEventSimDataCollection.resize(0);

	//initialize and fill muon detector
	fMuonDetector= new MuonDetector;
	fMuonDetector->SetNumberOfStrips( muonCounterDetector->GetNumberOfStripsInPlane() );
	fMuonDetector->SetNumberOfPlanes( muonCounterDetector->GetNumberOfPlanes() );
	fMuonDetector->SetSuperPlaneThickness( muonCounterDetector->GetSuperPlaneModuleThickness() );
	fMuonDetector->SetSuperPlaneHeight( muonCounterDetector->GetSuperPlaneModuleSizeZ() );
	fMuonDetector->SetSuperPlaneDistance( muonCounterDetector->GetDistanceAmongPlanes() );
	fMuonDetector->CalculateSuperPlaneDepth();
	fMuonDetector->SetXYPlaneDistance( muonCounterDetector->GetDistanceAmongXYPlanes() );
	fMuonDetector->SetXYPlaneTiltAngle( muonCounterDetector->GetPlaneTiltAngle() );
	fMuonDetector->SetStripSize( TVector3(muonCounterDetector->GetScintSize().x(),muonCounterDetector->GetScintSize().y(),muonCounterDetector->GetScintSize().z()) );
	fMuonDetector->SetStripGrooveSize( TVector3(muonCounterDetector->GetGrooveSize().x(),muonCounterDetector->GetGrooveSize().y(),muonCounterDetector->GetGrooveSize().z()) );
	fMuonDetector->SetStripCoatingSize( TVector3(muonCounterDetector->GetScintCoatingSize().x(),muonCounterDetector->GetScintCoatingSize().y(),muonCounterDetector->GetScintCoatingSize().z()) );
	fMuonDetector->SetStripCoatingThickness( muonCounterDetector->GetHousingThickness());
	fMuonDetector->SetFiberLength( muonCounterDetector->GetFiberSizeZ() );
	fMuonDetector->SetFiberRadius( muonCounterDetector->GetFiberRadius() );
	fMuonDetector->SetPMTSize( TVector3(muonCounterDetector->GetPMTSize().x(),muonCounterDetector->GetPMTSize().y(),muonCounterDetector->GetPMTSize().z()) );
	fMuonDetector->SetPhotocathodeSize( TVector3(muonCounterDetector->GetPhotocathodeSize().x(),muonCounterDetector->GetPhotocathodeSize().y(),muonCounterDetector->GetPhotocathodeSize().z()) );
	fMuonDetector->SetOptCouplingSize( TVector3(muonCounterDetector->GetPMTOpticalCouplingSize().x(),muonCounterDetector->GetPMTOpticalCouplingSize().y(),muonCounterDetector->GetPMTOpticalCouplingSize().z()) );

	if(fSaveSimInfo)
		fDetInfo->Fill();
	
  return eSuccess;
}


VModule::ResultFlag
G4MuonCounterSimulator::Run(Event& event)
{
	
	INFO("Running muon counter simulation");
	
	//## Steering with G4 macro
	if (fUseG4StandAloneConfiguration) {
    std::string command = "/control/execute ";
    fUImanager->ApplyCommand(command+fG4MacroFile);
  }


  return fFastMode ? RunFast(event) : RunFull(event);
}


VModule::ResultFlag
G4MuonCounterSimulator::RunFull(Event& event)
{
  INFO("Full muon counter simulation");
  
	fCurrentEvent= &event;

	// Get the SEvent
  SEvent& sEvent = event.GetSEvent();

  // Get the SDetector
  const SDetector& sDetector = Detector::GetInstance().GetSDetector();

  for (SEvent::StationIterator sIt = sEvent.StationsBegin(); sIt != sEvent.StationsEnd(); ++sIt) {

    if (!fDetectorConstructed) {
      ERROR("The detector hasn't been initialized when looping over stations");
      return eFailure;
    }

		int injectedParticles= 0;

    if (sIt->HasSimData()) {

      StationSimData& simData = sIt->GetSimData();

      // full-blown tank sim (no optimization of photon tracking)
      //simData.SetSimulatorSignature("G4TankSimulatorFullOG");

      const int stationId = sIt->GetId();

      const unsigned long numParticles = simData.ParticlesEnd() - simData.ParticlesBegin();

      if (numParticles)
        cerr << stationId << ':' << numParticles << ' ' << flush;
      else
        continue;

      fCurrentEventStationIt = sIt;
      fCurrentDetectorStation = &sDetector.GetStation(stationId);
      
			fCurrentParticleItList.clear();
      fCurrentParticleItList.resize(0);

			//## loop over shower particles
			//## fill station sim data in RunAction class
			if(fOneParticlePerStation){
				cout<<"Station "<<stationId<< ": injecting "<< numParticles<< " particles ...";
      	for (StationSimData::ParticleIterator pIt = simData.ParticlesBegin(); pIt != simData.ParticlesEnd(); ++pIt) {
				
        	fCurrentParticleIt = pIt;
        	fRunManager->BeamOn(1);

					injectedParticles++;
					if(injectedParticles%200 == 0) cout<<injectedParticles<<"...";
      	}//end loop injected particles
				cout<<"done"<<endl;		
			}//close if fOneParticlePerStation
			else{
				cout<<"Station "<<stationId<< ": injecting "<< numParticles<< " particles ...";
      	for (StationSimData::ParticleIterator pIt = simData.ParticlesBegin(); pIt != simData.ParticlesEnd(); ++pIt) {
					fCurrentParticleIt = pIt;
        	fCurrentParticleItList.push_back(pIt);
        	
					injectedParticles++;
      	}//end loop injected particles
		
				cout<<"Injecting "<<injectedParticles<<" particles...";
				fRunManager->BeamOn(1);
				cout<<"done"<<endl;
			}//close else
			
			

    }//if has sim data
  }//end loop stations

	//## set station sim data in EventSimData
	//fCurrentEventSimData->fStationSimDataCollection= fStationSimDataCollection;
	fEventSimDataCollection.push_back(*fCurrentEventSimData);

	//## fill sim data tree
	if(fSaveSimInfo)
		fSimInfo->Fill();

  return eSuccess;
}


VModule::ResultFlag
G4MuonCounterSimulator::RunFast(Event& event)
{
	INFO("Fast muon counter simulation");
	cout<<"No fast-mode implemented yet!"<<endl;
  return eSuccess;
}



VModule::ResultFlag
G4MuonCounterSimulator::Finish()
{

	if(fRecorder){
		delete fRecorder;
		fRecorder= 0;	
	}

	//#ifdef G4VIS_USE
	if(fVisManager) {
		delete fVisManager;
		fVisManager= 0;	
	}
	//#endif

  if (fRunManager) {
    delete fRunManager;
    fRunManager = 0;
  }


	if(fSaveSimInfo && fOutputFile) {
		fOutputFile->cd();
		fDetInfo->Write();
		fSimInfo->Write();
		fOutputFile->Close();
	}

  return eSuccess;
}



