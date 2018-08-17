/**
* @file G4MuonCounterTrajectory.hh
* @class G4MuonCounterTrajectory
* @brief Handle all trajectory tracks produced in the simulation 
*
* Set visualization attributes, colours, set drawing flags, ...
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterTrajectory_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterTrajectory_h 1

#include "G4Trajectory.hh"
#include "G4Allocator.hh"
#include "G4ios.hh" 
#include "globals.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"
#include "G4Track.hh"
#include "G4Step.hh"

#include "G4Colour.hh"

class G4Polyline;                   // Forward declaration.


namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterTrajectory : public G4Trajectory
{
public:
	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4MuonCounterTrajectory();
	/** 
	\brief Overloaded class constructor: initialize structures.
 	*/
  G4MuonCounterTrajectory(const G4Track* aTrack);	
	/** 
	\brief Overloaded class constructor: initialize structures.
 	*/
  G4MuonCounterTrajectory(G4MuonCounterTrajectory &);
	/**
	* \brief Class destructor: free allocated memory
	*/
  virtual ~G4MuonCounterTrajectory();
  /**
	* \brief Draw a track trajectory
	*/
  virtual void DrawTrajectory(G4int i_mode=0) const;
  /**
	* \brief Allocate new trajectory
	*/
  inline void* operator new(size_t);
	/**
	* \brief Delete a given trajectory
	*/
  inline void  operator delete(void*);
	/**
	* \brief Turns on/off the drawing of a given trajectory
	*/
  void SetDrawTrajectory(G4bool b){drawit=b;}
	/**
	* \brief Set this trajectory as a WLS photon track
	*/
  void WLS(){wls=true;}
	/**
	* \brief Set this trajectory as a Scintillation photon track
	*/
	void SCINT(){scint=true;}
	/**
	* \brief Set this trajectory as a Cerenkov photon track
	*/
	void CERENK(){cerenk=true;}
	/**
	* \brief Set this trajectory as a Rayleigh scattered photon track
	*/
  void SCATTERED(){scattered=true;}
	/**
	* \brief Force the drawing of the current trajectory
	*/
  void SetForceDrawTrajectory(G4bool b){forceDraw=b;}
	/**
	* \brief Exclude the drawing of the current trajectory
	*/
  void SetForceNoDrawTrajectory(G4bool b){forceNoDraw=b;}
  /**
	* \brief Set the colour of an optical photon trajectory
	*/
  void SetOptPhotonTrackColor(G4Colour col){OptPhotonColor=col;}
  
private:
  G4bool wls;
	G4bool cerenk;
	G4bool scint;
  G4bool scattered;
  G4bool drawit;
  G4bool forceNoDraw;
  G4bool forceDraw;
  G4ParticleDefinition* particleDefinition;
  
  mutable G4bool changeColorTrack;//define as "mutable" so that it can be modified by the "const" member DrawTrajectory(G4int i_mode=0) const
  mutable G4Colour OptPhotonColor;
  
};

}//close namespace


extern G4Allocator<G4MuonCounterSimulatorUSC::G4MuonCounterTrajectory> TrajectoryAllocator;


namespace G4MuonCounterSimulatorUSC {

inline void* G4MuonCounterTrajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)TrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void G4MuonCounterTrajectory::operator delete(void* aTrajectory)
{
  TrajectoryAllocator.FreeSingle((G4MuonCounterTrajectory*)aTrajectory);
}


}//close namespace

#endif
