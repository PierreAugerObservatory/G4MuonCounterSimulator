/**
* @file G4MuonCounterParticleInjector.hh
* @class G4MuonCounterParticleInjector
* @brief Inject particles in the muon station
* @author S. Riggi
* @date 08/02/2011
*/

#ifndef _G4MuonCounterParticleInjector_hh_
#define _G4MuonCounterParticleInjector_hh_

#include <fwk/VModule.h>

#include <utl/Particle.h>

#include <string>
#include <vector>


class TF1;
class TH1D;
class TH2D;
class TFile;

namespace sdet {
  class Station; 
}
namespace evt {
  class Event;
}
namespace sevt {
  class SEvent;
}
namespace utl {
  class RandomEngine;
}


class G4MuonCounterParticleInjector : public fwk::VModule {

  public:
    G4MuonCounterParticleInjector();
    virtual ~G4MuonCounterParticleInjector();

    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& theEvent);
    fwk::VModule::ResultFlag Finish();

    std::string GetSVNId() const 
      {return std::string("$Id$");}

  private:
    // switches
    bool fUseSingleTank;
    bool fUseDiscreteDirection;
    bool fUseDiscreteSpectrum;
    bool fUseInjectAllParticles;
    bool fUseSinglePosition;
    bool fUseDiscreteTime;
		bool fUseEnergyThetaDistributionFromFile;

		bool fUseEnergyDistributionFromFile;
		bool fUseZenithDistributionFromFile;
		bool fUseAzimuthDistributionFromFile;
		bool fUsePositionXDistributionFromFile;
		bool fUsePositionYDistributionFromFile;
		std::string fRandomDistributionFileName;

		TFile* fInputFile;
		TH2D* fMuonEnergyThetaDistribution;	
		TH1D* fMuonXDistribution;
		TH1D* fMuonYDistribution;
		TH1D* fMuonZenithDistribution;
		TH1D* fMuonAzimuthDistribution;
		TH1D* fMuonEnergyDistribution;
		TH1D* fAntiMuonXDistribution;
		TH1D* fAntiMuonYDistribution;
		TH1D* fAntiMuonZenithDistribution;
		TH1D* fAntiMuonAzimuthDistribution;
		TH1D* fAntiMuonEnergyDistribution;
		TH1D* fElectronXDistribution;
		TH1D* fElectronYDistribution;
		TH1D* fElectronZenithDistribution;
		TH1D* fElectronAzimuthDistribution;
		TH1D* fElectronEnergyDistribution;
		TH1D* fPositronXDistribution;
		TH1D* fPositronYDistribution;
		TH1D* fPositronZenithDistribution;
		TH1D* fPositronAzimuthDistribution;
		TH1D* fPositronEnergyDistribution;
		TH1D* fPhotonXDistribution;
		TH1D* fPhotonYDistribution;
		TH1D* fPhotonZenithDistribution;
		TH1D* fPhotonAzimuthDistribution;
		TH1D* fPhotonEnergyDistribution;
		TH1D* fProtonXDistribution;
		TH1D* fProtonYDistribution;
		TH1D* fProtonZenithDistribution;
		TH1D* fProtonAzimuthDistribution;
		TH1D* fProtonEnergyDistribution;
		TH1D* fNeutronXDistribution;
		TH1D* fNeutronYDistribution;
		TH1D* fNeutronZenithDistribution;
		TH1D* fNeutronAzimuthDistribution;
		TH1D* fNeutronEnergyDistribution;
		

    unsigned int fNumberOfParticles;
    unsigned int fSingleTankID;

    double fX, fY, fZ;
		double fGridX, fGridY, fGridZ;
		bool fRandomAroundPosition;

    std::vector<double> fDiscreteAzimuth;
    std::vector<double> fDiscreteZenith;
    std::vector<double> fDiscreteParticleTime;
    std::vector<double> fDiscreteMuonFlux;
    std::vector<double> fDiscreteElectronFlux;
    std::vector<double> fDiscretePhotonFlux;
		std::vector<double> fDiscreteHadronFlux;
    std::vector<double> fDiscreteMuonEnergy;
    std::vector<double> fDiscreteElectronEnergy;
    std::vector<double> fDiscretePhotonEnergy;
		std::vector<double> fDiscreteHadronEnergy;

    std::string fContinuousZenithString;
    std::string fContinuousAzimuthString;
    std::string fContinuousTimeString;
    std::string fMuonSpectrumString;
    std::string fElectronSpectrumString;
    std::string fPhotonSpectrumString;
		std::string fHadronSpectrumString;

    TF1* fContinuousZenithDistribution;
    TF1* fContinuousAzimuthDistribution;
    TF1* fTimeSpectrum;
    TF1* fMuonSpectrum;
    TF1* fElectronSpectrum;
    TF1* fPhotonSpectrum;
		TF1* fHadronSpectrum;

		double fZenithMin;
    double fZenithMax;
		double fAzimuthMin;
    double fAzimuthMax;
    double fMuonEnergyMin;
    double fMuonEnergyMax;
    double fElectronEnergyMin;
    double fElectronEnergyMax;
    double fPhotonEnergyMin;
    double fPhotonEnergyMax;
		double fHadronEnergyMin;
    double fHadronEnergyMax;
    double fMuonSpectrumMax;
    double fElectronSpectrumMax;
    double fPhotonSpectrumMax;
		double fHadronSpectrumMax;

    utl::Particle::Type fParticleType;

    utl::RandomEngine* fRandomEngine;

    void InjectParticles();

    double GetAzimuth();
    double GetZenith();
    double GetEnergy();
    double GetTime();

    void SetPosition(double& x, double& y, double& z);
    double GetDiscreteEnergy(std::vector<double>& flux, 
			     std::vector<double>& energy,
			     double max);

    const sdet::Station* fCurrentDetectorStation;
    sevt::SEvent* fSEvent;

    REGISTER_MODULE("G4MuonCounterParticleInjector", G4MuonCounterParticleInjector);

};



#endif
