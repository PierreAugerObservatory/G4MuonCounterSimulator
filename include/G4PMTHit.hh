/**
* @file G4PMTHit.hh
* @class G4PMTHit
* @brief Define the PMT hit structure
*
* A PMT hit contains information about the PMT Id, Strip Id, Plane Id, number of hit counts, arrival time and energy of * photons at the photocathode.
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4PMTHit_h
#define _G4MuonCounterSimulatorUSC_G4PMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

class G4VTouchable;

namespace G4MuonCounterSimulatorUSC {

class G4PMTHit : public G4VHit
{
public:
  
	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4PMTHit();
	/** 
	\brief Overloaded class constructor: initialize structures.
 	*/
	G4PMTHit(const G4PMTHit &right);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~G4PMTHit();

  /**
	* \brief Operator defining the equality operation between hits
	*/
  const G4PMTHit& operator=(const G4PMTHit &right);
	/**
	* \brief Tells if two hits are equal 
	*/
  int operator==(const G4PMTHit &right) const;
	
	/**
	* \brief Allocate a new hit
	*/
  inline void *operator new(size_t);
	/**
	* \brief Delete a given hit
	*/
  inline void operator delete(void *aHit);
  
	/**
	* \brief Handle the drawing and visualization attributes of PMT hits
	*/
  void Draw();
	/**
	* \brief Print information about hits
	*/
  void Print();

	/**
	* \brief Turns on/off the drawing of this hit
	*/
  inline void SetDrawit(bool b){drawit=b;}
	/**
	* \brief Tells if this hit has to be drawn
	*/
  inline bool GetDrawit(){return drawit;}

	/**
	* \brief Increment the photon counts of this hit (+1)
	*/
  inline void IncPhotonCount(){photons++;}
	/**
	* \brief Get the current photon counts
	*/
  inline int GetPhotonCount(){return photons;}

	/**
	* \brief Increment the scintillation photon counts of this hit (+1)
	*/
	inline void IncPhotonCount_scint(){photons_scint++;}
	/**
	* \brief Get the current scintillation photon counts
	*/
  inline int GetPhotonCount_scint(){return photons_scint;}

	/**
	* \brief Increment the cerenkov photon count of this hit (+1)
	*/
  inline void IncPhotonCount_cerenk(){photons_cerenk++;}
	/**
	* \brief Get the current cerenkov photon counts
	*/
  inline int GetPhotonCount_cerenk(){return photons_cerenk;}

	/**
	* \brief Increment the WLS photon count of this hit (+1)
	*/
	inline void IncPhotonCount_wls(){photons_wls++;}
	/**
	* \brief Get the current WLS photon counts
	*/
  inline int GetPhotonCount_wls(){return photons_wls;}

	/**
	* \brief Set the type of the current photon hitting the PMT
	*/
  inline void SetPhotonType(int val) { photons_type = val; }
	/**
	* \brief Get the type of the current photon hitting the PMT
	*/
  inline int GetPhotonType() { return photons_type; }

	/**
	* \brief Set the PMT id of the current hit
	*/
  inline void SetPMTNumber(int n) { pmtNumber = n; }
	/**
	* \brief Get the PMT id of the current hit
	*/
  inline int GetPMTNumber() { return pmtNumber; }
	/**
	* \brief Set the strip id of the current hit
	*/
  inline void SetStripNumber(int n) { stripNumber = n; }
	/**
	* \brief Get the strip id of the current hit
	*/
  inline int GetStripNumber() { return stripNumber; }
	/**
	* \brief Set the plane id of the current hit
	*/
	inline void SetPlaneNumber(int n) { planeNumber = n; }
	/**
	* \brief Get the plane id of the current hit
	*/
  inline int GetPlaneNumber() { return planeNumber; }

	/**
	* \brief Set the plane id of the current hit
	*/
	inline void SetSuperPlaneNumber(int n) { superplaneNumber = n; }
	/**
	* \brief Get the plane id of the current hit
	*/
  inline int GetSuperPlaneNumber() { return superplaneNumber; }

	/**
	* \brief Set the physical volume of the current hit
	*/
  inline void SetPMTPhysVol(G4VPhysicalVolume* physVol){this->physVol=physVol;}
	/**
	* \brief Get the physical volume of the current hit
	*/
  inline G4VPhysicalVolume* GetPMTPhysVol(){return physVol;}

	/**
	* \brief Set the position of the current hit
	*/
  inline void SetPMTPos(double x,double y,double z){
    pos=G4ThreeVector(x,y,z);
  }
	/**
	* \brief Get PMT position
	*/
  inline G4ThreeVector GetPMTPos(){return pos;}
	/**
	* \brief Set the position of the current hit
	*/
  inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
	/**
	* \brief Get the position of the current hit
	*/
  inline G4ThreeVector GetPos() { return pos; }
  
	/**
	* \brief Set the photon arrival times of the current hit
	*/
  inline void SetPhotonTime(std::vector<double> vect) { photonTime = vect; }
	/**
	* \brief Get the photon arrival times of the current hit
	*/
  inline std::vector<double> GetPhotonTime() { return photonTime; }
	/**
	* \brief Add a time entry to the photon arrival times of the current hit
	*/
  inline void InsertPhotonTime(double time) { photonTime.push_back(time); }
	/**
	* \brief Set the photon energies of the current hit
	*/
  inline void SetPhotonEnergy(std::vector<double> vect) { photonEnergy = vect; }
	/**
	* \brief Get the photon energies of the current hit
	*/
  inline std::vector<double> GetPhotonEnergy() { return photonEnergy; }
	/**
	* \brief Add a energy entry to the photon energies of the current hit
	*/
  inline void InsertPhotonEnergy(double energy) { photonEnergy.push_back(energy); }
  /**
	* \brief Set the photon arrival position of the current hit
	*/
  inline void SetPhotonPosition(std::vector<G4ThreeVector> vect) { photonPosition = vect; }
	/**
	* \brief Get the photon arrival position of the current hit
	*/
  inline std::vector<G4ThreeVector> GetPhotonPosition() { return photonPosition; }
	/**
	* \brief Add a position entry to the photon arrival position of the current hit
	*/
  inline void InsertPhotonPosition(G4ThreeVector position) { photonPosition.push_back(position); }	


private:
	/**
	* \brief PMT Id of the current hit
	*/
  int pmtNumber;
	/**
	* \brief Strip Id of the current hit
	*/
  int stripNumber;
	/**
	* \brief Plane Id of the current hit
	*/
	int planeNumber;
	/**
	* \brief Super Plane Id of the current hit
	*/
	int superplaneNumber; 
	/**
	* \brief Photon counts of the current hit
	*/
  int photons;
	/**
	* \brief Scintillation photon counts of the current hit
	*/
	int photons_scint;
	/**
	* \brief Cerenkov photon counts of the current hit
	*/
	int photons_cerenk;
	/**
	* \brief WLS photon counts of the current hit
	*/
	int photons_wls;
	/**
	* \brief Photon type of the current hit
	*/
  int photons_type;
	/**
	* \brief Photon hit position of the current hit
	*/
  G4ThreeVector pos;
	/**
	* \brief Physical volume of the current hit
	*/
  G4VPhysicalVolume* physVol;
	/**
	* \brief Drawing flag of the current hit
	*/
  bool drawit;
  /**
	* \brief List of photon arrival times of the current hit
	*/
  std::vector<double> photonTime;
	/**
	* \brief List of photon energies of the current hit
	*/
  std::vector<double> photonEnergy;
	/**
	* \brief List of photon positions of the current hit
	*/
  std::vector<G4ThreeVector> photonPosition;

};

typedef G4THitsCollection<G4PMTHit> PMTHitsCollection;

}//close namespace


extern G4Allocator<G4MuonCounterSimulatorUSC::G4PMTHit> PMTHitAllocator;


namespace G4MuonCounterSimulatorUSC {

	inline void* G4PMTHit::operator new(size_t){
  	void *aHit;
  	aHit = (void *) PMTHitAllocator.MallocSingle();
  	return aHit;
	}

	inline void G4PMTHit::operator delete(void *aHit){
  	PMTHitAllocator.FreeSingle((G4PMTHit*) aHit);
	}

}//close namespace

#endif


