/**
* @file G4ScintillatorHit.hh
* @class G4ScintillatorHit
* @brief Define the scintillator hit structure
*
* A scintillator hit contains information about the Strip Id, Plane Id, number of photons per process type emitted in the scintillator, energy deposit, position and arrival time of particles hitting the strip.
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4ScintillatorHit_h
#define _G4MuonCounterSimulatorUSC_G4ScintillatorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

namespace G4MuonCounterSimulatorUSC {

class G4ScintillatorHit : public G4VHit
{
public:
  
	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4ScintillatorHit();
	/** 
	\brief Overloaded class constructor: initialize structures.
 	*/
  G4ScintillatorHit(G4VPhysicalVolume* pVol);
	/** 
	\brief Overloaded class constructor: initialize structures.
 	*/
	G4ScintillatorHit(const G4ScintillatorHit &right);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~G4ScintillatorHit();
  
	/**
	* \brief Operator defining the equality operation between hits
	*/
  const G4ScintillatorHit& operator=(const G4ScintillatorHit &right);
	/**
	* \brief Tells if two hits are equal 
	*/
  int operator==(const G4ScintillatorHit &right) const;

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
	* \brief Set the energy deposit of the current hit
	*/
  inline void SetEdep(double de) { edep = de; }
	/**
	* \brief Add the current energy deposit to a previous hit
	*/
  inline void AddEdep(double de) { edep += de; }
	/**
	* \brief Get the current energy deposit of the current hit
	*/
  inline double GetEdep() { return edep; }
	/**
	* \brief Set the particle hit position of the current hit
	*/
  inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
	/**
	* \brief Get the particle hit position of the current hit
	*/
  inline G4ThreeVector GetPos() { return pos; }

	

	/**
	* \brief Set the photocathode surface position of the current hit
	*/
  inline void SetPhotocathodeSurfacePos(std::vector<G4ThreeVector> xyz) { photocathodeSurfacePos = xyz; }
	/**
	* \brief Get the photocathode surface position of the current hit
	*/
  inline std::vector<G4ThreeVector> GetPhotocathodeSurfacePos() { return photocathodeSurfacePos; }

	/**
	* \brief Set the expected time at the PMT anode of the current hit
	*/
  //inline void SetExpPMTTime(std::vector<double> vect) { expPMTTimeHit = vect; }
	/**
	* \brief Get the photocathode surface position of the current hit
	*/
  //inline std::vector<double> GetExpPMTTime() { return expPMTTimeHit; }

	/**
	* \brief Set the strip position of the current hit
	*/
  inline void SetStripPos(G4ThreeVector xyz) { strippos = xyz; }
	/**
	* \brief Get the strip position of the current hit
	*/
  inline G4ThreeVector GetStripPos() { return strippos; }
  /**
	* \brief Get the strip id of the current hit
	*/
  inline int GetSciCopyNo() { return stripNo; }
	/**
	* \brief Set the strip id of the current hit
	*/
  inline void SetSciCopyNo(int theCopyNo) { stripNo = theCopyNo; }
	/**
	* \brief Get the plane id of the current hit
	*/
	inline int GetPlaneCopyNo() { return planeNo; }
	/**
	* \brief Get the superplane id of the current hit
	*/
	inline int GetSuperPlaneCopyNo() { return superplaneNo; }
	/**
	* \brief Set the plane id of the current hit
	*/
  inline void SetPlaneCopyNo(int theCopyNo) { planeNo = theCopyNo; }
	/**
	* \brief Set the superplane id of the current hit
	*/
  inline void SetSuperPlaneCopyNo(int theCopyNo) { superplaneNo = theCopyNo; }
  /**
	* \brief Get the particle time of the current hit
	*/
  inline double GetTimeHit() { return timehit; }
	/**
	* \brief Set the particle time of the current hit
	*/
  inline void SetTimeHit(double theTimeHit) { timehit = theTimeHit; }
  /**
	* \brief Get the particle track id of the current hit
	*/
  inline int GetTrackId() { return trackid; }
	/**
	* \brief Set the particle track id of the current hit
	*/
  inline void SetTrackId(int theTrackId) { trackid = theTrackId; }
  /**
	* \brief Get the particle track type of the current hit
	*/
  inline int GetTrackPartType() { return trackparttype; }
	/**
	* \brief Set the particle track type of the current hit
	*/
  inline void SetTrackPartType(int theTrackPartType) { trackparttype = theTrackPartType; }


	/**
	* \brief Set the particle direction of the current hit
	*/
  inline void SetTrackDirection(G4ThreeVector dir) { direction = dir; }
	/**
	* \brief Get the particle direction of the current hit
	*/
  inline G4ThreeVector GetTrackDirection() { return direction; }

	/**
	* \brief Get the physical volume of the current hit
	*/
  inline const G4VPhysicalVolume * GetPhysV() { return physVol; }
	/**
	* \brief Set the physical volume of the current hit
	*/
  inline void SetPhysV(G4VPhysicalVolume * volume) { physVol = volume; }

private:

	/**
	* \brief Energy deposit of the current hit
	*/
  double edep;
	/**
	* \brief Position of the current hit
	*/
  G4ThreeVector pos;
	/**
	* \brief Direction of the current hit
	*/
  G4ThreeVector direction;
	/**
	* \brief Position of the current hit strip (center of the strip in global coordinates)
	*/
  G4ThreeVector strippos;
	/**
	* \brief Position of the PMT photocathode surface in current strip (in global coordinates)
	*/
  std::vector<G4ThreeVector> photocathodeSurfacePos;
	/**
	* \brief Expected time of the hit at the PMT anode
	*/
  //std::vector<double> expPMTTimeHit;
	/**
	* \brief Physical volume of the current hit
	*/
  const G4VPhysicalVolume* physVol;
	/**
	* \brief Strip Id of the current hit
	*/
  int stripNo;
	/**
	* \brief Plane Id of the current hit
	*/
	int planeNo;
	/**
	* \brief SuperPlane Id of the current hit
	*/
	int superplaneNo;
	/**
	* \brief Track time of the current hit
	*/
  double timehit;
	/**
	* \brief Track id of the current hit
	*/
  int trackid;
	/**
	* \brief Track type of the current hit
	*/
  int trackparttype;
  
};

typedef G4THitsCollection<G4ScintillatorHit> ScintHitsCollection;

}//close namespace


extern G4Allocator<G4MuonCounterSimulatorUSC::G4ScintillatorHit> ScintHitAllocator;


namespace G4MuonCounterSimulatorUSC {

	inline void* G4ScintillatorHit::operator new(size_t)
	{
  	void *aHit;
  	aHit = (void *) ScintHitAllocator.MallocSingle();
 		return aHit;
	}

	inline void G4ScintillatorHit::operator delete(void *aHit)
	{
 		ScintHitAllocator.FreeSingle((G4ScintillatorHit*) aHit);
	}


}//close namespace

#endif


