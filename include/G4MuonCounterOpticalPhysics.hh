/**
* @file G4MuonCounterOpticalPhysics.hh
* @class G4MuonCounterOpticalPhysics
* @brief Build the optical physical processes and optical particles involved in the simulation
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterOpticalPhysics_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterOpticalPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterOpticalPhysics : public G4VPhysicsConstructor
{
public: 
		/** 
		\brief Class constructor: initialize structures.
 		*/
    G4MuonCounterOpticalPhysics(const G4String& name ="optical");
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~G4MuonCounterOpticalPhysics();

public: 
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

		/**
		* \brief Set verbosity
		*/
    void SetVerbose(int);
		/**
		* \brief Set the max number of cerenkov photons
		*/
    void SetNbOfPhotonsCerenkov(int);
		/**
		* \brief Set the scintillation light yield factor
		*/
    void SetScintYieldFactor(double yf);

protected:
  
  G4Scintillation* theScintillationProcess;
  G4Cerenkov* theCerenkovProcess;
  G4OpAbsorption* theAbsorptionProcess;
  G4OpRayleigh* theRayleighScatteringProcess;
  G4OpBoundaryProcess* theBoundaryProcess;
  G4OpWLS* theWLSProcess;
};

}//close namespace

#endif

