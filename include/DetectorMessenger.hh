/**
* @file DetectorMessenger.hh
* @class DetectorMessenger
* @brief Create the UI commands to handle the detector geometry parameters, i.e. from a configuration file 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_DetectorMessenger_h
#define _G4MuonCounterSimulatorUSC_DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include "vector"


class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterConstruction;

class DetectorMessenger: public G4UImessenger
{
public:

	/** 
	\brief Class constructor: create UI commands
 	*/
  DetectorMessenger(G4MuonCounterConstruction*);
	/**
	* \brief Class destructor: delete allocated UI commands
	*/
  ~DetectorMessenger();
  
	/**
	* \brief Set (on-fly) the values parsed by the command to the detector geometry, by calling members of the detector construction class
	*/
  void SetNewValue(G4UIcommand*, G4String);
	/**
	* \brief Set (after config file reading) the values parsed by the command to the detector geometry, by calling members of the detector construction class. Called by DetectorConstruction::BuildDetectorGeometry()
	*/
	void SetFinalValues();
    
private:

  G4MuonCounterConstruction*   Detector;
  G4UIdirectory*               detectorDir;
  G4UIdirectory*               volumesDir;
	G4UIcommand*                 updateCmd;
	G4UIcmdWithAnInteger*        stripDesignModeCmd;
	G4UIcmdWithAnInteger*        stripReadoutModeCmd;
	G4UIcmdWithAnInteger*        stripOptCouplingModeCmd;
  G4UIcmdWith3VectorAndUnit*   stripSizeCmd;
	G4UIcmdWith3VectorAndUnit*   grooveSizeCmd;
	G4UIcmdWithADoubleAndUnit*   fiberLengthCmd;
	G4UIcmdWithADoubleAndUnit*   fiberRadiusCmd;
	G4UIcmdWithADoubleAndUnit*   fiberAbsLengthCmd;
  G4UIcmdWithADoubleAndUnit*   housingThicknessCmd;
  G4UIcmdWithADouble*          housingReflectivityCmd;
	G4UIcmdWithAnInteger*        housingSetSurfaceTypeCmd;
  G4UIcmdWithABool*            housingUseRealReflectivityCmd;
  G4UIcmdWith3VectorAndUnit*   pmtSizeCmd;
  G4UIcmdWith3VectorAndUnit*   photocathodeSizeCmd;
	G4UIcmdWith3VectorAndUnit*   pmtOptCouplingSizeCmd;
  G4UIcmdWithABool*            wlsCmd;
	G4UIcmdWithABool*            readoutModeCmd;
	G4UIcmdWithABool*            pmtUseRealQECmd;
	G4UIcmdWithABool*            yPlaneInSuperModuleCmd;
	G4UIcmdWithAnInteger*        stripNumberCmd;
	G4UIcmdWithAnInteger*        planeNumberCmd;
	G4UIcmdWithADoubleAndUnit*   planeDistanceCmd;
	G4UIcmdWithADoubleAndUnit*   planeXYDistanceCmd;
	G4UIcmdWithADoubleAndUnit*   planeTiltAngleCmd; 
	G4UIcmdWithADoubleAndUnit*   superplaneThicknessCmd; 
	G4UIcmdWithADoubleAndUnit*   superplaneSizeZCmd;  
	G4UIcmdWithADouble*          scintillationYieldCmd;  
  G4UIcmdWithADouble*          scintillationWeightCmd;  

	G4UIcmdWithAString*          worldMaterialCmd;	
	G4UIcmdWithAString*          superplaneCasingMaterialCmd;

	std::vector<G4double> argArray;
	
};

}//close namespace

#endif

