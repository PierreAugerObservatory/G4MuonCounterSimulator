/**
* @file PrimaryGeneratorMessenger.hh
* @class PrimaryGeneratorMessenger
* @brief Create the UI commands to handle the primary generator parameters, i.e. from a configuration file 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_PrimaryGeneratorMessenger_h
#define _G4MuonCounterSimulatorUSC_PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithABool;


namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterPrimaryGenerator;

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:

		/** 
		\brief Class constructor: create UI commands
 		*/
    PrimaryGeneratorMessenger(G4MuonCounterPrimaryGenerator*);
    /**
	  *\brief Class destructor: delete allocated UI commands
	  */
   ~PrimaryGeneratorMessenger();
    /**
		*\brief Set (on-fly) the values parsed by the command to the detector geometry, by calling members of the detector construction class
		*/	
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    
		G4MuonCounterPrimaryGenerator* Generator;
    G4UIdirectory*               genDir; 
    G4UIcmdWithADoubleAndUnit*   polarCmd;
		G4UIcmdWith3VectorAndUnit*   genPositionCmd;
		G4UIcmdWith3Vector*          genDirectionCmd;
		G4UIcmdWithADoubleAndUnit*   genEnergyCmd;	
		G4UIcmdWithADoubleAndUnit*   genMinEnergyCmd;	
		G4UIcmdWithADoubleAndUnit*   genMaxEnergyCmd;	
		G4UIcmdWithADoubleAndUnit*   genThetaCmd;	
		G4UIcmdWithADoubleAndUnit*   genMinThetaCmd;	
		G4UIcmdWithADoubleAndUnit*   genMaxThetaCmd;	
		G4UIcmdWithADoubleAndUnit*   genPhiCmd;	
		G4UIcmdWithADoubleAndUnit*   genMinPhiCmd;	
		G4UIcmdWithADoubleAndUnit*   genMaxPhiCmd;	
		G4UIcmdWithADoubleAndUnit*   genTimeCmd;
	  G4UIcmdWithAString*          genParticleCmd;
		G4UIcmdWithAnInteger*        genNumberOfParticlesCmd;
		//G4UIcmdWithABool*            genParticlesFromListCmd;
		G4UIcmdWithABool*            genParticlesFromShowerCmd;
		G4UIcmdWithABool*            genRandomEnergyCmd;
		G4UIcmdWithABool*            genRandomDirectionCmd;
		G4UIcmdWithABool*            genRandomDirectionUniformCmd;		
		G4UIcmdWithABool*            genRandomDirectionWithAcceptanceCmd;
		G4UIcmdWithABool*            genRandomPositionCmd;
		G4UIcmdWithABool*            genRandomPositionAroundVertexCmd;
		G4UIcmdWithABool*            genRandomPositionWithAcceptanceCmd;
		G4UIcmdWith3VectorAndUnit*   genRandomGridSizeCmd;
		G4UIcmdWithABool*            genCosmicMuonCmd;
		

};

}//close namespace

#endif

