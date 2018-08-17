/**
* @file G4MuonCounterPhysicsList.hh
* @class G4MuonCounterPhysicsList
* @brief Build the physical processes and particles involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterPhysicsList_h_
#define _G4MuonCounterSimulatorUSC_G4MuonCounterPhysicsList_h_ 1

#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterPhysicsList: public G4VModularPhysicsList
{
public:
	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4MuonCounterPhysicsList();
	/**
	* \brief Class destructor: free allocated memory
	*/
  virtual ~G4MuonCounterPhysicsList();
  /**
	* \brief Set cuts for physics processes
	*/
  virtual void SetCuts();
	/**
	* \brief Turns on/off the optical physics processes
	*/
	virtual void SetOpticalPhysics(bool);

private:

	bool fIncludeOpticalPhysics;


};

}//close namespace

#endif



