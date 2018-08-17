/**
* @file G4MuonCounterNeutronPhysics.hh
* @class G4MuonCounterNeutronPhysics
* @brief Build neutron particles and their physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterNeutronPhysics_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterNeutronPhysics_h 1

#include "G4VPhysicsConstructor.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterNeutronPhysics : public G4VPhysicsConstructor
{
  public: 
		/** 
		\brief Class constructor: initialize structures.
 		*/
    G4MuonCounterNeutronPhysics(const G4String& name ="neutron");
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~G4MuonCounterNeutronPhysics();

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





