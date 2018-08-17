/**
* @file SuperPlaneModule.hh
* @class SuperPlaneModule
* @brief Define a super plane module
*
* The superplane module contains two XY planes placed inside a casing
* @author Dr.ssa Enrica Trovato
* @date 15/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4SuperPlaneModule_h_
#define _G4MuonCounterSimulatorUSC_G4SuperPlaneModule_h_ 1

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


class G4SuperPlaneModule : public G4PVPlacement
{
public:

	/** 
	\brief Class constructor: create the detector structure of a plane
 	*/
  G4SuperPlaneModule(G4RotationMatrix *pRot,
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

	double fSuperPlaneSizeX;
  double fSuperPlaneSizeY;
  double fSuperPlaneSizeZ;

	double fSuperPlaneInsideSizeX;
  double fSuperPlaneInsideSizeY;
  double fSuperPlaneInsideSizeZ;

	double fSuperPlaneThickness;
	G4String fSuperPlaneMaterialName;
	bool fHasYPlaneInSuperModule;
	double fXYPlaneDistance;
	double fPlaneTiltAngle;	
	double fSurfaceSeparationGap;

	int fNstrips;

  double fScintSizeX;
  double fScintSizeY;
  double fScintSizeZ;

	double fFiberSizeZ;
	
	double fPMTOptCouplingSizeX;
  double fPMTOptCouplingSizeY;
	double fPMTOptCouplingSizeZ;

	double fScintCoatingSizeX;
  double fScintCoatingSizeY;
  double fScintCoatingSizeZ;
	double fScintCoatingThickness;
	double fScintCoatingReflectivity;


	double fPMTSizeX;
  double fPMTSizeY;
  double fPMTSizeZ;

	double fTiltAngle;


  //Basic Volumes
  //
	G4Box* ModuleHousing_box;
  G4Box* CasingInside_box;

  // Logical volumes
  //
  static G4LogicalVolume* ModuleHousing_log;
  static G4LogicalVolume* CasingInside_log;
  
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
	G4VPhysicalVolume* CasingInside_phys;
  

  //Sensitive Detectors
	//static G4MultiFunctionalDetector* MFDet;

};

}//close namespace

#endif
