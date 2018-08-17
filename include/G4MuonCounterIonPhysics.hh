/**
* @file G4MuonCounterIonPhysics.hh
* @class G4MuonCounterIonPhysics
* @brief Build ion particles and their physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterIonPhysics_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterIonPhysics_h 1

#include "G4VPhysicsConstructor.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterIonPhysics : public G4VPhysicsConstructor
{
  public: 
		/** 
		\brief Class constructor: initialize structures.
 		*/
    G4MuonCounterIonPhysics(const G4String& name ="ion");
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~G4MuonCounterIonPhysics();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
		/**
		* \brief Build all particles involved in the simulation
		*/
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
		/**
		* \brief Build all processes involved in the simulation
		*/
    virtual void ConstructProcess();

};

}//close namespace

#endif





