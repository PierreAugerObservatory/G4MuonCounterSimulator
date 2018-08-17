/**
* @file G4MuonCounterDecayPhysics.hh
* @class G4MuonCounterDecayPhysics
* @brief Build decay physical processes involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterDecayPhysics_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterDecayPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4Decay.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterDecayPhysics : public G4VPhysicsConstructor
{
  public: 
		/** 
		\brief Class constructor: initialize structures.
 		*/
    G4MuonCounterDecayPhysics(const G4String& name = "decay");
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~G4MuonCounterDecayPhysics();

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

  protected:
    G4Decay fDecayProcess;
};

}//close namespace

#endif








