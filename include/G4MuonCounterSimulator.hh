/**
* @mainpage G4MuonCounterSimulatorUSC tool
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
* @file G4MuonCounterSimulator.hh
* @class module
* @brief main module used to handle the entire project
*
* The simulation can run either in a batch mode or in interactive mode using visualization of detectors, tracks, etc...
* To select the run mode, just give in the command line the arguments "--batch" or "--interactive". The program will look for the "batchConfigFile.mac" or "interactiveConfigFile.mac" config files in the config_macros/ directory.
* As example, one can also handle a simulation output, for example giving as argument "--process".
* @author Dr. Simone Riggi
* @date 05/04/2010
*/



#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterSimulator_h_
#define _G4MuonCounterSimulatorUSC_G4MuonCounterSimulator_h_

#include <G4MuonCounterReconstructor.hh>
#include <G4MuonCounterStackingAction.hh>
#include <MuonDetector.hh>
#include <TStationSimData.hh>
#include <TEventSimData.hh>

#include <fwk/VModule.h>

#include <utl/TimeDistribution.h>

#include <sevt/SEvent.h>
#include <sevt/StationSimData.h>

#include <TFile.h>
#include <TTree.h>

class G4RunManager;
class G4UImanager;
class G4VisManager;

namespace sdet {
  class Station;
}

namespace det {
  class Detector;
}

/*
namespace mdet {
  class MStation;
	class MDetector;
}
*/

//using namespace G4MuonCounterReconstructorUSC;


namespace G4MuonCounterSimulatorUSC {

  class G4MuonCounterPrimaryGenerator;
  class G4MuonCounterStackingAction;
  class G4MuonCounterConstruction;
		
  // container of station id, weight
  typedef std::map<int, double> AccumulatedWeights;

  class G4MuonCounterSimulator : public fwk::VModule {

  public:
    G4MuonCounterSimulator();
    virtual ~G4MuonCounterSimulator();

    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& theEvent);
    fwk::VModule::ResultFlag Finish();

    std::string GetSVNId() const
    { return std::string("$Id$"); }

		static std::vector<TStationSimData> GetStationSimDataCollection()
    { return fStationSimDataCollection; }

		static std::vector<TEventSimData> GetEventSimDataCollection()
    { return fEventSimDataCollection; }

		static MuonDetector* GetCurrentMuonDetector()
    { return fMuonDetector; }
	
		static TEventSimData* GetCurrentEventSimData()
    { return fCurrentEventSimData; }


  private:
    void ConstructTraces(sevt::Station& station) const;

    static const sdet::Station* GetCurrentDetectorStation()
    { return fCurrentDetectorStation; }

    static sevt::SEvent::StationIterator GetCurrentEventStationIt()
    { return fCurrentEventStationIt; }

		static const evt::Event* GetCurrentEvent()
    { return fCurrentEvent; }

    static sevt::StationSimData::ParticleIterator GetCurrentParticleIt()
    { return fCurrentParticleIt; }

		static std::vector<sevt::StationSimData::ParticleIterator> GetCurrentParticleItList()
    { return fCurrentParticleItList; }

    static sevt::Station::SignalComponent GetCurrentComponent()
    { return GetComponentId(fCurrentParticleIt); }

		static TStationSimData* GetCurrentStationSimData()
    { return fCurrentStationSimData; }

		static bool OpticalOn()
    { return fgOptical; }

	  static bool ScintillationOn()
    { return fgScintillation; }

		static bool CherenkovOn()
    { return fgCherenkov; }

		static bool WLSOn()
    { return fgWLS; }

		static bool AbsorptionOn()
    { return fgAbsorption; }

		static bool RayleighOn()
    { return fgRayleigh; }

		static bool BoundaryOn()
    { return fgBoundary; }

    static bool MuCaptureOn()
    { return fgMuCapture; }


		
    /// retrieve the station signal component from particle type
    static sevt::Station::SignalComponent 
    GetComponentId(const sevt::StationSimData::ConstParticleIterator currentParticle);

    fwk::VModule::ResultFlag RunFull(evt::Event& theEvent);
    fwk::VModule::ResultFlag RunFast(evt::Event& theEvent);
		
    static const sdet::Station* fCurrentDetectorStation;
    static sevt::SEvent::StationIterator fCurrentEventStationIt;
    static sevt::StationSimData::ParticleIterator fCurrentParticleIt;
    static sevt::Station::SignalComponent fgCurrentComponent; 
		static const evt::Event* fCurrentEvent;

		static std::vector<sevt::StationSimData::ParticleIterator> fCurrentParticleItList;
		

    G4RunManager* fRunManager;
    G4UImanager* fUImanager;
    G4VisManager* fVisManager;
		G4MuonCounterRecorderBase* fRecorder;	
		
    G4MuonCounterStackingAction* fStackingAction;

    bool fGeoVisOn;
    bool fTrajVisOn;
		std::string fG4MacroFile;
		bool fUseG4StandAloneConfiguration;
		int fVisDriverType;
		int fRunVerbosity;
		int fTrackingVerbosity;
		int fEventVerbosity;		

    bool fDetectorConstructed;
    bool fFastMode;
		bool fFullMode;
	  static bool fOneParticlePerStation;

		static bool fgOptical;
    static bool fgMuCapture;
		static bool fgScintillation;
		static bool fgWLS;
		static bool fgCherenkov;
		static bool fgAbsorption;
		static bool fgRayleigh;
		static bool fgBoundary;
		static bool fgAllProcess;	
		
    static G4MuonCounterConstruction* fMuonCounterConstruction;
		

    std::string fEventId;

		bool fSaveSimInfo;
		bool fSaveFullInfo;
		TFile* fOutputFile;
		std::string fOutputFileName;
		TTree* fSimInfo;
		TTree* fDetInfo;				

	  static TStationSimData* fCurrentStationSimData;
		static std::vector<TStationSimData> fStationSimDataCollection;
		static TEventSimData* fCurrentEventSimData;
		static std::vector<TEventSimData> fEventSimDataCollection;
		static MuonDetector* fMuonDetector;

    friend class G4MuonCounterPrimaryGenerator;
    friend class G4MuonCounterConstruction;
    friend class G4MuonCounterPhysicsList;
		friend class G4MuonCounterOpticalPhysics;
		friend class G4MuonCounterRunAction;
		
    REGISTER_MODULE("G4MuonCounterSimulatorUSC", G4MuonCounterSimulator);

  };

}

#endif
