/**
* @file G4MuonCounterPrimaryGenerator.hh
* @class G4MuonCounterPrimaryGenerator
* @brief Generate the primary particles in the simulation
* @author S. Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterPrimaryGenerator_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterPrimaryGenerator_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterConstruction;
class PrimaryGeneratorMessenger;

class G4MuonCounterPrimaryGenerator : public G4VUserPrimaryGeneratorAction {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    G4MuonCounterPrimaryGenerator();
		/** 
		\brief Overloaded class constructor: initialize structures.
 		*/
		G4MuonCounterPrimaryGenerator(G4MuonCounterConstruction*);
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~G4MuonCounterPrimaryGenerator();

  public:
		/**
		* \brief Generate the primary particles (vertex,energy,momentum,time,...)
		*/
    void GeneratePrimaries(G4Event*);
		/**
		* \brief A debug generator
		*/
    void DebugGenerator(G4Event*);
	
		/**
		* \brief Standard generator used to inject resampled shower particles in the muon counters
		*/
    void ShowerGenerator(G4Event*);
    /**
		* \brief Standard generator used to inject resampled shower particles from a list in the muon counters
		*/
    void ShowerGeneratorList(G4Event*);
		
		/**
		* \brief Set random polarization of optical photon
		*/
    void SetOptPhotonPolar();
		/**
		* \brief Set polarization of optical photon with given angle
		*/
    void SetOptPhotonPolar(double);
		
		/**
		* \brief Set the number of primary particles to be generated in the event
		*/
		void SetNumberOfGenParticles(int value) {fNpartPerEvent= value;}
		/**
		* \brief Get the number of primary particles to be generated in the event
		*/
		int GetNumberOfGenParticles() {return fNpartPerEvent;}

    /**
		* \brief Set theta of primary particles to be generated in the event
		*/
    void SetThetaOfGenParticles(double value) {fGenTheta= value;}
    /**
		* \brief Get theta of primary particles to be generated in the event
		*/
		double GetThetaOfGenParticles() {return fGenTheta;}
		/**
		* \brief Set the minimum theta of primary particles to be generated in the event
		*/
		void SetMinThetaOfGenParticles(double value) {fMinGenTheta= value;}
		/**
		* \brief Get the minimum theta of primary particles to be generated in the event
		*/
		double GetMinThetaOfGenParticles() {return fMinGenTheta;}
		/**
		* \brief Set the maximum theta of primary particles to be generated in the event
		*/
		void SetMaxThetaOfGenParticles(double value) {fMaxGenTheta= value;}
		/**
		* \brief Get the maximum theta of primary particles to be generated in the event
		*/
		double GetMaxThetaOfGenParticles() {return fMaxGenTheta;}

    /**
		* \brief Set phi of primary particles to be generated in the event
		*/
    void SetPhiOfGenParticles(double value) {fGenPhi= value;}
    /**
		* \brief Get phi of primary particles to be generated in the event
		*/
		double GetPhiOfGenParticles() {return fGenPhi;}
		/**
		* \brief Set the minimum phi of primary particles to be generated in the event
		*/
		void SetMinPhiOfGenParticles(double value) {fMinGenPhi= value;}
		/**
		* \brief Get the minimum phi of primary particles to be generated in the event
		*/
		double GetMinPhiOfGenParticles() {return fMinGenPhi;}
		/**
		* \brief Set the maximum phi of primary particles to be generated in the event
		*/
		void SetMaxPhiOfGenParticles(double value) {fMaxGenPhi= value;}
		/**
		* \brief Get the maximum phi of primary particles to be generated in the event
		*/
		double GetMaxPhiOfGenParticles() {return fMaxGenPhi;}



		/**
		* \brief Generate the primary energy according to the spectrum of a Sr-90 beta source
		*/
		double GetEnergyFromSr90Spectrum();

		/**
		* \brief Set the primary particle to be generated in the event
		*/
		void SetGenParticle(G4ParticleDefinition* value) {fGenParticle= value;}
		/**
		* \brief Get the primary particle to be generated in the event
		*/
		G4ParticleDefinition* GetGenParticle() {return fGenParticle;}

		/**
		* \brief Set the primary particle energy
		*/
		void SetGenEnergy(double value) {fGenEnergy= value;}
		/**
		* \brief Get the primary particle energy
		*/
		double GetGenEnergy() {return fGenEnergy;}
		/**
		* \brief Set the minimum primary particle energy
		*/
		void SetMinGenEnergy(double value) {fMinGenEnergy= value;}
		/**
		* \brief Get the minimum primary particle energy
		*/
		double GetMinGenEnergy() {return fMinGenEnergy;}
		/**
		* \brief Set the maximum primary particle energy
		*/
		void SetMaxGenEnergy(double value) {fMaxGenEnergy= value;}
		/**
		* \brief Get the maximum primary particle energy
		*/
		double GetMaxGenEnergy() {return fMaxGenEnergy;}

		/**
		* \brief Set the primary particle time
		*/
		void SetGenTime(double value) {fGenTime= value;}
		/**
		* \brief Get the primary particle time
		*/
		double GetGenTime() {return fGenTime;}

		/**
		* \brief Set the primary particle vertex position
		*/
		void SetGenPosition(G4ThreeVector value) {fGenPosition= value;}
		/**
		* \brief Get the primary particle vertex position
		*/
		G4ThreeVector GetGenPosition() {return fGenPosition;}
		
		/**
		* \brief Set the primary particle direction
		*/
		void SetGenDirection(G4ThreeVector value) {fGenDirection= value;}
		/**
		* \brief Get the primary particle direction
		*/
		G4ThreeVector GetGenDirection() {return fGenDirection;}
		


		/**
		* \brief Set the list of primary particle to be generated in the event
		*/
		void SetGenParticleList(std::vector<G4ParticleDefinition*> vect) {fGenParticleList= vect;}
		/**
		* \brief Get the list of primary particle to be generated in the event
		*/
		std::vector<G4ParticleDefinition*> GetGenParticleList() {return fGenParticleList;}
		/**
		* \brief Set the list of primary particle energies
		*/
		void SetGenEnergyList(std::vector<double> vect) {fGenEnergyList= vect;}
		/**
		* \brief Get the list of primary particle energies
		*/
		std::vector<double> GetGenEnergyList() {return fGenEnergyList;}
		/**
		* \brief Set the list of primary particle times
		*/
		void SetGenTimeList(std::vector<double> vect) {fGenTimeList= vect;}
		/**
		* \brief Get the list of primary particle times
		*/
		std::vector<double> GetGenTimeList() {return fGenTimeList;}
		/**
		* \brief Set the list of primary particle vertex positions
		*/
		void SetGenPositionList(std::vector<G4ThreeVector> vect) {fGenPositionList= vect;}
		/**
		* \brief Get the list of primary particle vertex positions
		*/
		std::vector<G4ThreeVector> GetGenPositionList() {return fGenPositionList;}

		/**
		* \brief Set the list of primary particle directions
		*/
		void SetGenDirectionList(std::vector<G4ThreeVector> vect) {fGenDirectionList= vect;}
		/**
		* \brief Get the list of primary particle directions
		*/
		std::vector<G4ThreeVector> GetGenDirectionList() {return fGenDirectionList;}

		/**
		* \brief Turns on/off the reading of primary particles from resampled shower particles attached to each station
		*/
		void SetParticleInjectionFromResampledShower(bool value) {fUseShowerGenerator= value;}
		/**
		* \brief Get the flag of read particles from resampled showers 
		*/
		bool GetParticleInjectionFromResampledShower() {return fUseShowerGenerator;}

		/**
		* \brief Enable/Disable the randomize mode of primary particle positions (in a squared grid)
		*/
		void SetRandomizingPosition(bool value) {fRandomizePosition= value;}
		/**
		* \brief Get the flag of the randomize mode of primary particle positions
		*/
		bool GetRandomizingPosition() {return fRandomizePosition;}
		/**
		* \brief Set size of the random injection area 
		*/
		void SetSizeOfRandomInjectionArea(G4ThreeVector vect) {
			fSizeOfRandomAreaX= vect.x();
			fSizeOfRandomAreaY= vect.y();
			fSizeOfRandomAreaZ= vect.z();
		}
		/**
		* \brief Get size of the random injection area 
		*/
		G4ThreeVector GetSizeOfRandomInjectionArea() {return G4ThreeVector(fSizeOfRandomAreaX,fSizeOfRandomAreaY,fSizeOfRandomAreaZ);}
		

		/**
		* \brief Enable/Disable the randomize mode of primary particle positions (in a squared grid)
		*/
		void SetRandomizingPositionAroundVertex(bool value) {fRandomizePositionAroundVertex= value;}
		/**
		* \brief Get the flag of the randomize mode of primary particle positions
		*/
		bool GetRandomizingPositionAroundVertex() {return fRandomizePositionAroundVertex;}
		/**
		* \brief Enable/Disable the randomize mode of primary particle position (viewable position according to detector plane geometrical acceptance)
		*/
		void SetRandomizingPositionWithAcceptance(bool value) {fRandomizePositionWithAcceptance= value;}
		/**
		* \brief Get the flag of the randomize mode of primary particle directions (directions according to detector plane acceptance)
		*/
		bool GetRandomizingPositionWithAcceptance() {return fRandomizePositionWithAcceptance;}
		
		/**
		* \brief Enable/Disable the randomize mode of primary particle directions (uniform directions in space)
		*/
		void SetRandomizingDirectionUniform(bool value) {fRandomizeDirectionUniform= value;}
		/**
		* \brief Get the flag of the randomize mode of primary particle directions (uniform directions in space)
		*/
		bool GetRandomizingDirectionUniform() {return fRandomizeDirectionUniform;}

		/**
		* \brief Enable/Disable the randomize mode of primary particle directions
		*/
		void SetRandomizingDirection(bool value) {fRandomizeDirection= value;}
		/**
		* \brief Get the flag of the randomize mode of primary particle directions 		
		*/
		bool GetRandomizingDirection() {return fRandomizeDirection;}


		/**
		* \brief Enable/Disable the randomize mode of primary particle directions (directions according to detector plane acceptance)
		*/
		void SetRandomizingDirectionWithAcceptance(bool value) {fRandomizeDirectionWithAcceptance= value;}
		/**
		* \brief Get the flag of the randomize mode of primary particle directions (directions according to detector plane acceptance)
		*/
		bool GetRandomizingDirectionWithAcceptance() {return fRandomizeDirectionWithAcceptance;}

		/**
		* \brief Enable/Disable the randomize mode of primary particle energy (uniform in given range)
		*/
		void SetRandomizingEnergy(bool value) {fRandomizeEnergy= value;}
		/**
		* \brief Get the flag of the randomize mode of primary particle energy (uniform in given range)
		*/
		bool GetRandomizingEnergy() {return fRandomizeEnergy;}
		/**
		* \brief Enable/Disable the cosmic muon generator
		*/
		void UseCosmicMuonGenerator(bool choice) {fUseCosmicMuonGenerator= choice;}


  private:
    G4ParticleGun* particleGun;
		G4ParticleTable* particleTable;    
    PrimaryGeneratorMessenger* gunMessenger;

		int fNpartPerEvent;
		G4ParticleDefinition* fGenParticle;
		double fGenEnergy;
		double fMinGenEnergy;
		double fMaxGenEnergy;
		double fGenTime;
		G4ThreeVector fGenPosition;
		G4ThreeVector fGenDirection;
    double fGenTheta;
		double fMinGenTheta;
		double fMaxGenTheta;
    double fGenPhi;
		double fMinGenPhi;
		double fMaxGenPhi;

    G4MuonCounterConstruction* Detector;

		bool fUseShowerGenerator;	
		bool fUseCosmicMuonGenerator;
		bool fRandomizePosition;
		bool fRandomizePositionAroundVertex;
		bool fRandomizePositionWithAcceptance;
		double fSizeOfRandomAreaX;
		double fSizeOfRandomAreaY;
		double fSizeOfRandomAreaZ;
	
    bool fDirection;
		bool fRandomizeDirection;
		bool fRandomizeDirectionUniform;
		bool fRandomizeDirectionWithAcceptance;
		bool fRandomizeEnergy;

		std::vector<double> fGenEnergyList;
		std::vector<double> fGenTimeList;
		std::vector<G4ThreeVector> fGenDirectionList;
		std::vector<G4ThreeVector> fGenPositionList;
		std::vector<G4ParticleDefinition*> fGenParticleList;
		
		
};

}//close namespace

#endif /*PrimaryGeneratorAction_h*/
