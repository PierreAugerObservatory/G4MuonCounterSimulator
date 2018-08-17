/**
* @file G4PlaneModule.hh
* @class G4PlaneModule
* @brief Define a detector plane
*
* The detector plane module is made up of a series of strip modules, one aside to the other. This class defines the plane design in a "object-oriented" way, so that it can be allocated by the detector construction. The allocation of the strip module class in done in this class, while the allocation of the detector planes is done in the detector costruction class.
* @author S. Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4PlaneModule_h_
#define _G4MuonCounterSimulatorUSC_G4PlaneModule_h_ 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4MultiFunctionalDetector.hh"

#include "G4MuonCounterConstruction.hh"

namespace G4MuonCounterSimulatorUSC {

class G4PlaneModule : public G4PVPlacement
{
public:

	/** 
	\brief Class constructor: create the detector structure of a plane
 	*/
  G4PlaneModule(G4RotationMatrix *pRot,
		const G4ThreeVector &tlate,
		G4LogicalVolume *pMotherLogical,
		bool pMany,
		int pCopyNo,
		G4MuonCounterConstruction* c);

	/**
	* \brief ### TO BE REMOVED ###
	*/
	void SetHousingVolume(G4LogicalVolume* volume_log){ModuleHousing_log= volume_log;}

	/**
	* \brief ### TO BE REMOVED ###
	*/
	void SetTiltAngle(double value){fTiltAngle=value;}
	/**
	* \brief ### TO BE REMOVED ###
	*/
	double GetTiltAngle(){return fTiltAngle;}


private:

	/**
	* \brief Set the visualization attributes of the plane (colours, visibility, etc...)
	*/
  void VisAttributes();
	/**
	* \brief Set the surface properties of the detector plane components
	*/
  void SurfaceProperties();
	/**
	* \brief Read the geometry parameters from the detector construction and pass to this class
	*/
  void CopyValues();

  bool fUpdated;
  
  G4MuonCounterConstruction* fDetectorConstructor;

	int fStripDesignMode;
	int fStripReadoutMode;
	int fStripOptCouplingMode;

	double fPlaneSizeX;
  double fPlaneSizeY;
  double fPlaneSizeZ;

	int fNstrips;

  double fScintSizeX;
  double fScintSizeY;
  double fScintSizeZ;

	double fFiberSizeZ;

	double fScintCoatingSizeX;
  double fScintCoatingSizeY;
  double fScintCoatingSizeZ;
	double fScintCoatingThickness;
	double fScintCoatingReflectivity;

	double fPMTOptCouplingSizeX;
  double fPMTOptCouplingSizeY;
	double fPMTOptCouplingSizeZ;

	double fPMTSizeX;
  double fPMTSizeY;
  double fPMTSizeZ;

	double fTiltAngle;
	double fSurfaceSeparationGap;

  //Basic Volumes
  //
	G4Box* ModuleHousing_box;
  

  // Logical volumes
  //
  static G4LogicalVolume* ModuleHousing_log;
  
  
  // Physical volumes
  //
  G4VPhysicalVolume* Scint_phys;
	G4VPhysicalVolume* ScintCoating_phys;
	G4VPhysicalVolume* Groove_phys;	
	G4VPhysicalVolume* Fiber_phys;
	G4VPhysicalVolume* FiberClad1_phys;
  G4VPhysicalVolume* FiberClad2_phys;
  G4VPhysicalVolume* PMT_phys;
  G4VPhysicalVolume* Photocathode_phys;

  
  //Sensitive Detectors
	static G4MultiFunctionalDetector* MFDet;

};

}//close namespace

#endif
