/**
* @file G4MuonCounterParticleInjector.cc
* @class G4MuonCounterParticleInjector
* @brief Inject particles in the muon station
* @author S. Riggi
* @date 08/02/2011
*/

#include <utl/config.h>

#include "G4MuonCounterParticleInjector.hh"

#include <det/Detector.h>

#include <sdet/SDetector.h>
#include <sdet/Station.h>

#include <evt/Event.h>

#include <fwk/CentralConfig.h>
#include <fwk/RandomEngineRegistry.h>

#include <sevt/SEvent.h>
#include <sevt/Station.h>
#include <sevt/StationSimData.h>
#include <sevt/Header.h>

#include <utl/AugerCoordinateSystem.h>
#include <utl/AugerUnits.h>
#include <utl/MathConstants.h>
#include <utl/Particle.h>
#include <utl/PhysicalConstants.h>
#include <utl/Point.h>
#include <utl/TimeStamp.h>
#include <utl/ErrorLogger.h>
#include <utl/RandomEngine.h>
#include <utl/Reader.h>
#include <utl/Vector.h>

#include <CLHEP/Random/Randomize.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TRandom3.h>

#include <cstddef>
#include <sstream>
#include <vector>
#include <iostream>
#include <time.h>


using namespace det;
using namespace sdet;
using namespace evt;
using namespace fwk;
using namespace sevt;
using namespace utl;
using namespace std;

using CLHEP::RandFlat;


static utl::Particle::Type spectrumParticles[7] = { utl::Particle::ePhoton,
                                                    utl::Particle::eElectron,
                                                    utl::Particle::ePositron,
                                                    utl::Particle::eAntiMuon,
                                                    utl::Particle::eMuon,
																										utl::Particle::eProton,
																										utl::Particle::eNeutron
																									 };

static utl::Particle::Type spectrumMainParticles[3] = { utl::Particle::eMuon,
                                                    		utl::Particle::eElectron,
                                                    		utl::Particle::eProton
																									 		};


G4MuonCounterParticleInjector::G4MuonCounterParticleInjector() :
  fUseSingleTank(false),
  fUseDiscreteDirection(true),
  fUseDiscreteSpectrum(true),
  fUseInjectAllParticles(false),
  fUseSinglePosition(false),
  fUseDiscreteTime(true),
	fUseEnergyThetaDistributionFromFile(false),
	fUseEnergyDistributionFromFile(false),
	fUsePositionXDistributionFromFile(false),
	fUsePositionYDistributionFromFile(false),
	fUseZenithDistributionFromFile(false),
	fUseAzimuthDistributionFromFile(false),
  fNumberOfParticles(0),
  fSingleTankID(0),
  fX(0.0), fY(0.0), fZ(0.0),
	fGridX(0.0), fGridY(0.0), fGridZ(0.0),
	fRandomAroundPosition(false),
  //fContinuousZenithString(""),
  //fContinuousAzimuthString(""),
  //fContinuousTimeString(""),
  //fMuonSpectrumString(""),
  //fElectronSpectrumString(""),
  //fPhotonSpectrumString(""),
  fContinuousZenithDistribution(NULL),
  fContinuousAzimuthDistribution(NULL),
  fTimeSpectrum(NULL),
  fMuonSpectrum(NULL),
  fElectronSpectrum(NULL),
  fPhotonSpectrum(NULL),
	fHadronSpectrum(NULL),
	fInputFile(NULL),
	fMuonEnergyThetaDistribution(NULL),
	fMuonXDistribution(NULL),
	fMuonYDistribution(NULL),
	fMuonZenithDistribution(NULL),
	fMuonAzimuthDistribution(NULL),
	fMuonEnergyDistribution(NULL),
	fAntiMuonXDistribution(NULL),
	fAntiMuonYDistribution(NULL),
	fAntiMuonZenithDistribution(NULL),
	fAntiMuonAzimuthDistribution(NULL),
	fAntiMuonEnergyDistribution(NULL),
	fElectronXDistribution(NULL),
	fElectronYDistribution(NULL),
	fElectronZenithDistribution(NULL),
	fElectronAzimuthDistribution(NULL),
	fElectronEnergyDistribution(NULL),
	fPositronXDistribution(NULL),
	fPositronYDistribution(NULL),
	fPositronZenithDistribution(NULL),
	fPositronAzimuthDistribution(NULL),
	fPositronEnergyDistribution(NULL),
	fPhotonXDistribution(NULL),
	fPhotonYDistribution(NULL),
	fPhotonZenithDistribution(NULL),
	fPhotonAzimuthDistribution(NULL),
	fPhotonEnergyDistribution(NULL),
	fProtonXDistribution(NULL),
	fProtonYDistribution(NULL),
	fProtonZenithDistribution(NULL),
	fProtonAzimuthDistribution(NULL),
	fProtonEnergyDistribution(NULL),
	fNeutronXDistribution(NULL),
	fNeutronYDistribution(NULL),
	fNeutronZenithDistribution(NULL),
	fNeutronAzimuthDistribution(NULL),
	fNeutronEnergyDistribution(NULL),
  fMuonEnergyMin(0.0), fMuonEnergyMax(0.0),
  fElectronEnergyMin(0.0), fElectronEnergyMax(0.0),
  fPhotonEnergyMin(0.0), fPhotonEnergyMax(0.0),
	fHadronEnergyMin(0.0), fHadronEnergyMax(0.0),
  fMuonSpectrumMax(0.0), fElectronSpectrumMax(0.0), fPhotonSpectrumMax(0.0), fHadronSpectrumMax(0.0),
  fParticleType(utl::Particle::eUndefined),
  fRandomEngine(0),
  fCurrentDetectorStation(NULL), fSEvent(NULL)
{
}


G4MuonCounterParticleInjector::~G4MuonCounterParticleInjector()
{
  delete fTimeSpectrum;
  delete fContinuousAzimuthDistribution;
  delete fMuonSpectrum;
  delete fElectronSpectrum;
  delete fPhotonSpectrum;
	delete fHadronSpectrum;
  delete fContinuousZenithDistribution;

	if(fMuonEnergyThetaDistribution) delete fMuonEnergyThetaDistribution;
	delete fMuonXDistribution;
	delete fMuonYDistribution;
	delete fMuonZenithDistribution;
	delete fMuonAzimuthDistribution;
	delete fMuonEnergyDistribution;
	delete fAntiMuonXDistribution;
	delete fAntiMuonYDistribution;
	delete fAntiMuonZenithDistribution;
	delete fAntiMuonAzimuthDistribution;
	delete fAntiMuonEnergyDistribution;
	delete fElectronXDistribution;
	delete fElectronYDistribution;
	delete fElectronZenithDistribution;
	delete fElectronAzimuthDistribution;
	delete fElectronEnergyDistribution;
	delete fPositronXDistribution;
	delete fPositronYDistribution;
	delete fPositronZenithDistribution;
	delete fPositronAzimuthDistribution;
	delete fPositronEnergyDistribution;
	delete fPhotonXDistribution;
	delete fPhotonYDistribution;
	delete fPhotonZenithDistribution;
	delete fPhotonAzimuthDistribution;
	delete fPhotonEnergyDistribution;
	delete fProtonXDistribution;
	delete fProtonYDistribution;
	delete fProtonZenithDistribution;
	delete fProtonAzimuthDistribution;
	delete fProtonEnergyDistribution;
	delete fNeutronXDistribution;
	delete fNeutronYDistribution;
	delete fNeutronZenithDistribution;
	delete fNeutronAzimuthDistribution;
	delete fNeutronEnergyDistribution;
}


VModule::ResultFlag
G4MuonCounterParticleInjector::Init()
{
  INFO(".");

	//## Set random seed
	unsigned long int seed1= time(0);
	unsigned long int seed2= clock();
  
	timespec timeStruct;
	int temp;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timeStruct);
	unsigned long int seed3= timeStruct.tv_nsec;

	delete gRandom;
	gRandom= new TRandom3(seed1+seed2+seed3);


//#warning dv: ParticleInjector needs xsd!
  Branch topB =
    CentralConfig::GetInstance()->GetTopBranch("G4MuonCounterParticleInjector");

  Branch configB = topB.GetChild("config");

  // Reading config branch
  configB.GetChild("UseSingleTank").GetData(fUseSingleTank);
  configB.GetChild("UseInjectAllParticles").GetData(fUseInjectAllParticles);
  configB.GetChild("UseSinglePosition").GetData(fUseSinglePosition);
  configB.GetChild("UseDiscreteDirection").GetData(fUseDiscreteDirection);
  configB.GetChild("UseDiscreteEnergy").GetData(fUseDiscreteSpectrum);
  configB.GetChild("UseDiscreteTime").GetData(fUseDiscreteTime);
	configB.GetChild("UseRandomPositionXDistributionFromFile").GetData(fUsePositionXDistributionFromFile);
	configB.GetChild("UseRandomPositionYDistributionFromFile").GetData(fUsePositionYDistributionFromFile);
	configB.GetChild("UseRandomZenithDistributionFromFile").GetData(fUseZenithDistributionFromFile);
	configB.GetChild("UseRandomAzimuthDistributionFromFile").GetData(fUseAzimuthDistributionFromFile);
	configB.GetChild("UseRandomEnergyDistributionFromFile").GetData(fUseEnergyDistributionFromFile);
	configB.GetChild("UseRandomEnergyThetaDistributionFromFile").GetData(fUseEnergyThetaDistributionFromFile);
	configB.GetChild("randomDistributionFile").GetData(fRandomDistributionFileName);
  cerr << "UseInjectAllParticles = " << fUseInjectAllParticles << '\n'
       << "UseSinglePosition = " << fUseSinglePosition << '\n'
       << "UseDiscreteDirection = " << fUseDiscreteDirection << '\n'
       << "UseDiscreteSpectrum = " << fUseDiscreteSpectrum << '\n'
       << "UseDiscreteTime = " << fUseDiscreteTime << '\n'
			 << " UseEnergyDistributionFromFile = " << fUseEnergyDistributionFromFile << "  ["<<fRandomDistributionFileName.c_str()<<"]"<<endl;

	if(fUsePositionXDistributionFromFile || fUsePositionYDistributionFromFile ||
     fUseZenithDistributionFromFile || fUseAzimuthDistributionFromFile ||
     fUseEnergyDistributionFromFile)
	{
		fInputFile= new TFile(fRandomDistributionFileName.c_str(),"READ");
	
		fMuonEnergyThetaDistribution= (TH2D*)fInputFile->Get("MuonDistrHisto");

		fMuonXDistribution= (TH1D*)fInputFile->Get("muonX");
		fMuonYDistribution= (TH1D*)fInputFile->Get("muonY");
		fMuonZenithDistribution= (TH1D*)fInputFile->Get("muonTheta");
		fMuonAzimuthDistribution= (TH1D*)fInputFile->Get("muonPhi");
		fMuonEnergyDistribution= (TH1D*)fInputFile->Get("muonEnergy");
		fAntiMuonXDistribution= (TH1D*)fInputFile->Get("antimuonX");
		fAntiMuonYDistribution= (TH1D*)fInputFile->Get("antimuonY");
		fAntiMuonZenithDistribution= (TH1D*)fInputFile->Get("antimuonTheta");
		fAntiMuonAzimuthDistribution= (TH1D*)fInputFile->Get("antimuonPhi");
		fAntiMuonEnergyDistribution= (TH1D*)fInputFile->Get("antimuonEnergy");
		fElectronXDistribution= (TH1D*)fInputFile->Get("electronX");
		fElectronYDistribution= (TH1D*)fInputFile->Get("electronY");
		fElectronZenithDistribution= (TH1D*)fInputFile->Get("electronTheta");
		fElectronAzimuthDistribution= (TH1D*)fInputFile->Get("electronPhi"); 
		fElectronEnergyDistribution= (TH1D*)fInputFile->Get("electronEnergy");
		fPositronXDistribution= (TH1D*)fInputFile->Get("positronX");
		fPositronYDistribution= (TH1D*)fInputFile->Get("positronY");
		fPositronZenithDistribution= (TH1D*)fInputFile->Get("positronTheta");
		fPositronAzimuthDistribution= (TH1D*)fInputFile->Get("positronPhi"); 
		fPositronEnergyDistribution= (TH1D*)fInputFile->Get("positronEnergy");
		fPhotonXDistribution= (TH1D*)fInputFile->Get("gammaX");
		fPhotonYDistribution= (TH1D*)fInputFile->Get("gammaY");
		fPhotonZenithDistribution= (TH1D*)fInputFile->Get("gammaTheta");
		fPhotonAzimuthDistribution= (TH1D*)fInputFile->Get("gammaPhi");
		fPhotonEnergyDistribution= (TH1D*)fInputFile->Get("gammaEnergy");
		fProtonXDistribution= (TH1D*)fInputFile->Get("protonX");
		fProtonYDistribution= (TH1D*)fInputFile->Get("protonY");
		fProtonZenithDistribution= (TH1D*)fInputFile->Get("protonTheta");
		fProtonAzimuthDistribution= (TH1D*)fInputFile->Get("protonPhi");
		fProtonEnergyDistribution= (TH1D*)fInputFile->Get("protonEnergy");
		fNeutronXDistribution= (TH1D*)fInputFile->Get("neutronX");
		fNeutronYDistribution= (TH1D*)fInputFile->Get("neutronY");
		fNeutronZenithDistribution= (TH1D*)fInputFile->Get("neutronTheta");
		fNeutronAzimuthDistribution= (TH1D*)fInputFile->Get("neutronPhi");
		fNeutronEnergyDistribution= (TH1D*)fInputFile->Get("neutronEnergy");
	}

  topB.GetChild("NumberOfParticles").GetData(fNumberOfParticles);

  // Single Position
  if (fUseSinglePosition) {
    Branch singlePositionB = topB.GetChild("SinglePosition");
    singlePositionB.GetChild("ParticleX").GetData(fX);
    singlePositionB.GetChild("ParticleY").GetData(fY);
    singlePositionB.GetChild("ParticleZ").GetData(fZ);
  }
	// Random Position
  else {
    Branch randomPositionB = topB.GetChild("RandomPosition");
    randomPositionB.GetChild("ParticleX").GetData(fX);
    randomPositionB.GetChild("ParticleY").GetData(fY);
    randomPositionB.GetChild("ParticleZ").GetData(fZ);
		randomPositionB.GetChild("GridX").GetData(fGridX);
    randomPositionB.GetChild("GridY").GetData(fGridY);
    randomPositionB.GetChild("GridZ").GetData(fGridZ);
		randomPositionB.GetChild("RandomAroundPosition").GetData(fRandomAroundPosition);
  }


  // Direction Branch
  if (!fUseDiscreteDirection) {
    //Continuous distribution
    Branch cdsB = topB.GetChild("UseContinuousDirectionSpectrum");
    //zenith angle
    cdsB.GetChild("Zenith").GetData(fContinuousZenithString);
		cdsB.GetChild("ZenithMin").GetData(fZenithMin);
		cdsB.GetChild("ZenithMax").GetData(fZenithMax);
		fZenithMin/= utl::rad;	
		fZenithMax/= utl::rad;
    //fContinuousZenithDistribution = new TF1("Zenith", fContinuousZenithString.c_str(), 0.0, utl::kPi/2.0);
		fContinuousZenithDistribution = new TF1("Zenith", fContinuousZenithString.c_str(), fZenithMin, fZenithMax);
    
		//azimuth angle
    cdsB.GetChild("Azimuth").GetData(fContinuousAzimuthString);
		cdsB.GetChild("AzimuthMin").GetData(fAzimuthMin);
		cdsB.GetChild("AzimuthMax").GetData(fAzimuthMax);
		fAzimuthMin/= utl::rad;	
		fAzimuthMax/= utl::rad;
    //fContinuousAzimuthDistribution = new TF1("Azimuth", fContinuousAzimuthString.c_str(), 0.0, 2.0*utl::kPi);
		fContinuousAzimuthDistribution = new TF1("Azimuth", fContinuousAzimuthString.c_str(), fAzimuthMin, fAzimuthMax);
  } 
	else {
    //Discrete distribution
    Branch ddsB = topB.GetChild("UseDiscreteDirectionSpectrum");
    //zenith angle
    ddsB.GetChild("Zenith").GetData(fDiscreteZenith);
    //azimuth angle
    ddsB.GetChild("Azimuth").GetData(fDiscreteAzimuth);
  }

  // Energy Branch
  if (!fUseDiscreteSpectrum) {

    //Continuous distribution
    Branch cesB = topB.GetChild("UseContinuousEnergySpectrum");

    cesB.GetChild("MuonSpectrum").GetData(fMuonSpectrumString);
    cesB.GetChild("MuonEnergyMin").GetData(fMuonEnergyMin);
    cesB.GetChild("MuonEnergyMax").GetData(fMuonEnergyMax);
    fMuonSpectrum = new TF1("MuonSpectrum", fMuonSpectrumString.c_str(),fMuonEnergyMin/utl::GeV, fMuonEnergyMax/utl::GeV);

    cesB.GetChild("ElectronSpectrum").GetData(fElectronSpectrumString);
    cesB.GetChild("ElectronEnergyMin").GetData(fElectronEnergyMin);
    cesB.GetChild("ElectronEnergyMax").GetData(fElectronEnergyMax);
    fElectronSpectrum = new TF1("ElectronSpectrum", fElectronSpectrumString.c_str(),fElectronEnergyMin/utl::GeV, fElectronEnergyMax/utl::GeV);

    cesB.GetChild("PhotonSpectrum").GetData(fPhotonSpectrumString);
    cesB.GetChild("PhotonEnergyMin").GetData(fPhotonEnergyMin);
    cesB.GetChild("PhotonEnergyMax").GetData(fPhotonEnergyMax);
    fPhotonSpectrum = new TF1("PhotonSpectrum", fPhotonSpectrumString.c_str(),fPhotonEnergyMin/utl::GeV, fPhotonEnergyMax/utl::GeV);

		cesB.GetChild("HadronSpectrum").GetData(fHadronSpectrumString);
    cesB.GetChild("HadronEnergyMin").GetData(fHadronEnergyMin);
    cesB.GetChild("HadronEnergyMax").GetData(fHadronEnergyMax);
    fHadronSpectrum = new TF1("HadronSpectrum", fHadronSpectrumString.c_str(),fHadronEnergyMin/utl::GeV, fHadronEnergyMax/utl::GeV);

  } 
	else {

    //Discrete distribution
    Branch desB = topB.GetChild("UseDiscreteEnergySpectrum");

    desB.GetChild("MuonSpectrum").GetData(fDiscreteMuonFlux);
    fMuonSpectrumMax = -10;
    for (unsigned int i = 0; i < fDiscreteMuonFlux.size(); ++i)
      if (fDiscreteMuonFlux[i] > fMuonSpectrumMax)
        fMuonSpectrumMax = fDiscreteMuonFlux[i];
    desB.GetChild("MuonEnergy").GetData(fDiscreteMuonEnergy);

    desB.GetChild("ElectronSpectrum").GetData(fDiscreteElectronFlux);
    fElectronSpectrumMax = -10;
    for (unsigned int i = 0; i < fDiscreteElectronFlux.size(); ++i)
      if (fDiscreteElectronFlux[i] > fElectronSpectrumMax)
        fElectronSpectrumMax = fDiscreteElectronFlux[i];
    desB.GetChild("ElectronEnergy").GetData(fDiscreteElectronEnergy);

    desB.GetChild("PhotonSpectrum").GetData(fDiscretePhotonFlux);
    fPhotonSpectrumMax = -10;
    for (unsigned int i = 0; i < fDiscretePhotonFlux.size(); i++)
      if (fDiscretePhotonFlux[i] > fPhotonSpectrumMax)
        fPhotonSpectrumMax = fDiscretePhotonFlux[i];
    desB.GetChild("PhotonEnergy").GetData(fDiscretePhotonEnergy);

		desB.GetChild("HadronSpectrum").GetData(fDiscreteHadronFlux);
    fHadronSpectrumMax = -10;
    for (unsigned int i = 0; i < fDiscreteHadronFlux.size(); i++)
      if (fDiscreteHadronFlux[i] > fHadronSpectrumMax)
        fHadronSpectrumMax = fDiscreteHadronFlux[i];
    desB.GetChild("HadronEnergy").GetData(fDiscreteHadronEnergy);

  }

  topB.GetChild("SingleTankID").GetData(fSingleTankID);

  int particleType;
  topB.GetChild("ParticleType").GetData(particleType);
  fParticleType = static_cast<utl::Particle::Type>(particleType);

  if (fUseDiscreteTime)
    topB.GetChild("UseDiscreteTimeSpectrum").GetChild("ParticleTime").GetData(fDiscreteParticleTime);

  fRandomEngine = &RandomEngineRegistry::GetInstance().Get(RandomEngineRegistry::eDetector);

  return eSuccess;
}


VModule::ResultFlag
G4MuonCounterParticleInjector::Run(Event& theEvent)
{
  INFO(".");

  if (!theEvent.HasSEvent()) {
    ERROR("No SEvent present. Use EventGenerator module to create "
          "an event before running G4MuonCounterParticleInjector.");
    return eFailure;
  }

  fSEvent = &theEvent.GetSEvent();

  const sdet::SDetector& theSDetector = det::Detector::GetInstance().GetSDetector();

  if (fUseSingleTank) {
    fCurrentDetectorStation = &theSDetector.GetStation(fSingleTankID);
    InjectParticles();
  } 
	else{
    for (SDetector::StationIterator sdIt = theSDetector.StationsBegin();sdIt != theSDetector.StationsEnd(); ++sdIt) {
      fCurrentDetectorStation = &theSDetector.GetStation(sdIt->GetId());
      InjectParticles();
    }
	}

  return eSuccess;
}


VModule::ResultFlag
G4MuonCounterParticleInjector::Finish()
{
  INFO(".");
  delete fContinuousZenithDistribution;
  fContinuousZenithDistribution = 0;
  delete fContinuousAzimuthDistribution;
  fContinuousAzimuthDistribution = 0;
  delete fMuonSpectrum;
  fMuonSpectrum = 0;
  delete fElectronSpectrum;
  fElectronSpectrum = 0;
  delete fPhotonSpectrum;
  fPhotonSpectrum = 0;
	delete fHadronSpectrum;
  fHadronSpectrum = 0;

	if( (fUsePositionXDistributionFromFile || fUsePositionYDistributionFromFile ||
     	 fUseZenithDistributionFromFile || fUseAzimuthDistributionFromFile ||
       fUseEnergyDistributionFromFile) && fInputFile)
	{	

		
		if(fMuonEnergyThetaDistribution) {
			delete fMuonEnergyThetaDistribution;
			fMuonEnergyThetaDistribution= 0;
		}
		
		delete fMuonXDistribution;
		fMuonXDistribution= 0;
		delete fMuonYDistribution;
		fMuonYDistribution= 0;	
		delete fMuonZenithDistribution;
		fMuonZenithDistribution= 0;
		delete fMuonAzimuthDistribution;
		fMuonAzimuthDistribution= 0;
		delete fMuonEnergyDistribution;
		fMuonEnergyDistribution= 0;

		delete fAntiMuonXDistribution;
		fAntiMuonXDistribution= 0;
		delete fAntiMuonYDistribution;
		fAntiMuonYDistribution= 0;	
		delete fAntiMuonZenithDistribution;
		fAntiMuonZenithDistribution= 0;
		delete fAntiMuonAzimuthDistribution;
		fAntiMuonAzimuthDistribution= 0;
		delete fAntiMuonEnergyDistribution;
		fAntiMuonEnergyDistribution= 0;

		delete fElectronXDistribution;
		fElectronXDistribution= 0;
		delete fElectronYDistribution;
		fElectronYDistribution= 0;
		delete fElectronZenithDistribution;
		fElectronZenithDistribution= 0;
		delete fElectronAzimuthDistribution;
		fElectronAzimuthDistribution= 0;
		delete fElectronEnergyDistribution;
		fElectronEnergyDistribution= 0;

		delete fPositronXDistribution;
		fPositronXDistribution= 0;
		delete fPositronYDistribution;
		fPositronYDistribution= 0;
		delete fPositronZenithDistribution;
		fPositronZenithDistribution= 0;
		delete fPositronAzimuthDistribution;
		fPositronAzimuthDistribution= 0;
		delete fPositronEnergyDistribution;
		fPositronEnergyDistribution= 0;

		delete fPhotonXDistribution;
		fPhotonXDistribution= 0;
		delete fPhotonYDistribution;
		fPhotonYDistribution= 0;
		delete fPhotonZenithDistribution;
		fPhotonZenithDistribution= 0;
		delete fPhotonAzimuthDistribution;
		fPhotonAzimuthDistribution= 0;
		delete fPhotonEnergyDistribution;
		fPhotonEnergyDistribution= 0;
		delete fProtonXDistribution;
		fProtonXDistribution= 0;
		delete fProtonYDistribution;
		fProtonYDistribution= 0;
		delete fProtonZenithDistribution;
		fProtonZenithDistribution= 0;
		delete fProtonAzimuthDistribution;
		fProtonAzimuthDistribution= 0;
		delete fProtonEnergyDistribution;
		fProtonEnergyDistribution= 0;
		delete fNeutronXDistribution;
		fNeutronXDistribution= 0;
		delete fNeutronYDistribution;
		fNeutronYDistribution= 0;
		delete fNeutronZenithDistribution;
		fNeutronZenithDistribution= 0;
		delete fNeutronAzimuthDistribution;
		fNeutronAzimuthDistribution= 0;
		delete fNeutronEnergyDistribution;
		fNeutronEnergyDistribution= 0;

		//fInputFile->Close();
	}

  return eSuccess;
}


void
G4MuonCounterParticleInjector::InjectParticles()
{
  if (!fSEvent->HasStation(fCurrentDetectorStation->GetId()))
    fSEvent->MakeStation(fCurrentDetectorStation->GetId());

  sevt::Station& currentEventStation =
    fSEvent->GetStation(fCurrentDetectorStation->GetId());

  if (!currentEventStation.HasSimData())
    currentEventStation.MakeSimData();

  const CoordinateSystemPtr cs = fCurrentDetectorStation->GetLocalCoordinateSystem();

  for (unsigned int num = 0; num < fNumberOfParticles; ++num) {

    double azimuth = GetAzimuth();
    double zenith = GetZenith();
		double energy = GetEnergy();
		
		if(fUseEnergyThetaDistributionFromFile){
			double randLgE, randTheta;
			fMuonEnergyThetaDistribution->GetRandom2(randTheta,randLgE);
			energy= pow(10,randLgE)*utl::GeV;
			zenith= randTheta*utl::deg;
		}	
		
		//cout<<"theta="<<zenith/utl::degree<<"  phi="<<azimuth/utl::degree<<endl;
			
		
    double x;
    double y;
    double z;
    SetPosition(x, y, z);

    const Point position(x, y, z, cs);

    const double sinZenith = sin(zenith);
    const double pX = -sinZenith*cos(azimuth);
    const double pY = -sinZenith*sin(azimuth);
    const double pZ = -cos(zenith);

    const Vector direction(pX, pY, pZ, cs);
		
		cout<<"G4MuonCounterParticleInjector::InjectParticles(): INFO: E="<<energy/utl::GeV<<endl;
		cout<<"G4MuonCounterParticleInjector::InjectParticles(): INFO: theta= "<<zenith/utl::degree<<"  phi="<<azimuth/utl::degree<<endl;
		cout<<"G4MuonCounterParticleInjector::InjectParticles(): INFO: gen particle dir: ("<<pX<<","<<pY<<","<<pZ<<")"<<endl;

    const TimeInterval time(GetTime());

    if (fUseInjectAllParticles) {
      //const int index = int(RandFlat::shoot(&fRandomEngine->GetEngine(), 0, 7));
      //fParticleType = spectrumParticles[index];
	
			const int index = int(RandFlat::shoot(&fRandomEngine->GetEngine(), 0, 2));
			
			fParticleType = spectrumMainParticles[index];
			cout<<"G4MuonCounterParticleInjector::InjectParticles(): index="<<index<<"  particle="<<fParticleType<<endl;
    }
		cout<<"G4MuonCounterParticleInjector::InjectParticles(): INFO: ParticleType="<<fParticleType<<endl;
    
		

    const Particle newParticle(fParticleType, utl::Particle::eBackground,
                               position, direction, time, 1, energy);

    currentEventStation.GetSimData().AddParticle(newParticle);

  }
}


double
G4MuonCounterParticleInjector::GetAzimuth(){

  if (fUseDiscreteDirection) {
    const unsigned int num = fDiscreteAzimuth.size();
    unsigned int index;

    do
      index = (unsigned int)(RandFlat::shoot(&fRandomEngine->GetEngine(), 0, num));
    while (index >= num);

    return fDiscreteAzimuth[index];

  } 
	else{	
		if(fUseAzimuthDistributionFromFile) {
			switch (fParticleType) {
  			case utl::Particle::eMuon:
    			return fMuonAzimuthDistribution->GetRandom()*utl::deg;
    		break;

				case utl::Particle::eAntiMuon:
    			return fAntiMuonAzimuthDistribution->GetRandom()*utl::deg;
    		break;

  			case utl::Particle::eElectron:
  				return fElectronAzimuthDistribution->GetRandom()*utl::deg;
    		break;

				case utl::Particle::ePositron:
    			return fPositronAzimuthDistribution->GetRandom()*utl::deg;
    		break;

				case utl::Particle::eProton:
					return fProtonAzimuthDistribution->GetRandom()*utl::deg;
    		break;

				case utl::Particle::eNeutron:
					return fNeutronAzimuthDistribution->GetRandom()*utl::deg;
    		break;

  			case utl::Particle::ePhoton:
					return fPhotonAzimuthDistribution->GetRandom()*utl::deg;
    		break;
  		
				default:
    			INFO("Unknown particle type used for spectrum injection in ParticleInjector.");
    			return 0;
    		break;
 	 		}//close switch
	
		}//close if UseDistributionFromFile
		else{
			return fContinuousAzimuthDistribution->GetRandom();
		}
	}//close continuous distribution

}//close GetAzimuth()


double
G4MuonCounterParticleInjector::GetZenith()
{
  if (fUseDiscreteDirection) {
    const unsigned int num = fDiscreteZenith.size();
    unsigned int index;

    do
      index = (unsigned int)(RandFlat::shoot(&fRandomEngine->GetEngine(), 0, num));
    while (index >= num);

    return fDiscreteZenith[index];

  } 
	else{	
		if(fUseZenithDistributionFromFile) {
			switch (fParticleType) {
  			case utl::Particle::eMuon:
    			return fMuonZenithDistribution->GetRandom()*utl::deg;
    		break;

				case utl::Particle::eAntiMuon:
    			return fAntiMuonZenithDistribution->GetRandom()*utl::deg;
    		break;

  			case utl::Particle::eElectron:
  				return fElectronZenithDistribution->GetRandom()*utl::deg;
    		break;

				case utl::Particle::ePositron:
    			return fPositronZenithDistribution->GetRandom()*utl::deg;
    		break;

				case utl::Particle::eProton:
  				return fProtonZenithDistribution->GetRandom()*utl::deg;
    		break;

				case utl::Particle::eNeutron:
					return fNeutronZenithDistribution->GetRandom()*utl::deg;
    		break;

  			case utl::Particle::ePhoton:
					return fPhotonZenithDistribution->GetRandom()*utl::deg;
    		break;
  		
				default:
    			INFO("Unknown particle type used for spectrum injection in ParticleInjector.");
    			return 0;
    		break;
 	 		}//close switch
	
		}//close if UseDistributionFromFile
		else{
			return fContinuousZenithDistribution->GetRandom();
		}
	}//close continuous distribution
	
}


double
G4MuonCounterParticleInjector::GetEnergy(){

  switch (fParticleType) {

		case utl::Particle::eMuon:
			return fUseDiscreteSpectrum ?
      	GetDiscreteEnergy(fDiscreteMuonFlux, fDiscreteMuonEnergy, fMuonSpectrumMax) :
					fUseEnergyDistributionFromFile ?
						pow(10,fMuonEnergyDistribution->GetRandom())*utl::GeV :
        		fMuonSpectrum->GetRandom()*GeV;
    	break;

  	case utl::Particle::eAntiMuon:
    	return fUseDiscreteSpectrum ?
      	GetDiscreteEnergy(fDiscreteMuonFlux, fDiscreteMuonEnergy, fMuonSpectrumMax) :
					fUseEnergyDistributionFromFile ?
						pow(10,fAntiMuonEnergyDistribution->GetRandom())*utl::GeV :
        		fMuonSpectrum->GetRandom()*GeV;
    	break;

  	case utl::Particle::eElectron:
			return fUseDiscreteSpectrum ?
    	  GetDiscreteEnergy(fDiscreteElectronFlux, fDiscreteElectronEnergy, fElectronSpectrumMax) :
					fUseEnergyDistributionFromFile ?
						pow(10,fElectronEnergyDistribution->GetRandom())*utl::GeV :
        		fElectronSpectrum->GetRandom()*GeV;
    	break;

  	case utl::Particle::ePositron:

    	return fUseDiscreteSpectrum ?
    	  GetDiscreteEnergy(fDiscreteElectronFlux, fDiscreteElectronEnergy, fElectronSpectrumMax) :
					fUseEnergyDistributionFromFile ?
						pow(10,fPositronEnergyDistribution->GetRandom())*utl::GeV :
        		fElectronSpectrum->GetRandom()*GeV;
    	break;

		case utl::Particle::eProton:
			return fUseDiscreteSpectrum ?
     		GetDiscreteEnergy(fDiscreteHadronFlux, fDiscreteHadronEnergy, fHadronSpectrumMax) :
					fUseEnergyDistributionFromFile ?
						pow(10,fProtonEnergyDistribution->GetRandom())*utl::GeV :	
        		fHadronSpectrum->GetRandom()*GeV;
    	break;

  	case utl::Particle::eNeutron:

    	return fUseDiscreteSpectrum ?
     		GetDiscreteEnergy(fDiscreteHadronFlux, fDiscreteHadronEnergy, fHadronSpectrumMax) :
					fUseEnergyDistributionFromFile ?
						pow(10,fNeutronEnergyDistribution->GetRandom())*utl::GeV :	
        		fHadronSpectrum->GetRandom()*GeV;
    	break;

  	case utl::Particle::ePhoton:

    	return fUseDiscreteSpectrum ?
    	  GetDiscreteEnergy(fDiscretePhotonFlux, fDiscretePhotonEnergy, fPhotonSpectrumMax) :	
					fUseEnergyDistributionFromFile ?
						pow(10,fPhotonEnergyDistribution->GetRandom())*utl::GeV :	
        		fPhotonSpectrum->GetRandom()*GeV;
    	break;

  	default:

    	INFO("Unknown particle type used for spectrum injection in ParticleInjector.");
    	return 0;
    	break;
  }//close switch


}


double
G4MuonCounterParticleInjector::GetTime()
{
  const unsigned int num = fDiscreteParticleTime.size();
  unsigned int index;

  do
    index = (unsigned int)(RandFlat::shoot(&fRandomEngine->GetEngine(), 0, num));
  while (index >= num);

  return fDiscreteParticleTime[index];
}


void
G4MuonCounterParticleInjector::SetPosition(double& x, double& y, double& z)
{
  const double height = fCurrentDetectorStation->GetHeight();
  const double radius = fCurrentDetectorStation->GetRadius();

  if (fUseSinglePosition) {
    x = fX;
    y = fY;
    z = fZ;
  } 
	else {
		x= RandFlat::shoot(&fRandomEngine->GetEngine(), -fGridX/2., fGridX/2.);
		y= RandFlat::shoot(&fRandomEngine->GetEngine(), -fGridY/2., fGridY/2.);
		z= RandFlat::shoot(&fRandomEngine->GetEngine(), -fGridZ/2., fGridZ/2.);
		
		//randomize around chosen vertex?
		if(fRandomAroundPosition){
			x+= fX;
			y+= fY;
			z+= fZ;
		}
  }

}


double
G4MuonCounterParticleInjector::GetDiscreteEnergy(std::vector<double>& flux,
                                    std::vector<double>& energy,
                                    const double max)
{
  if (flux.size() != energy.size())

    INFO("Flux and energy array sizes not equal!");

  else if (max <= 0)

    INFO("Non-positive maximum value for flux given.");

  else {

    const unsigned int num = flux.size();
    unsigned int index;

    do
      index = (unsigned int)(RandFlat::shoot(&fRandomEngine->GetEngine(), 0, num));
    while (index >= num || RandFlat::shoot(&fRandomEngine->GetEngine(), 0, max) >= flux[index]);

    return energy[index];

  }

  return 0;
}


