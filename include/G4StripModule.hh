/**
* @file G4StripModule.hh
* @class G4StripModule
* @brief Define a strip module
*
* The strip module is made up of a scintillator bar, surrounded by a reflective coating. A groove is done in the bar to host a WLS fiber. The coating is done by an upper and side layer (always present) plus the coating end-caps. When both fiber and scintillator are read out no end-caps are inserted. When only the fiber is read out an end-cap with a hole (same size of the WLS fiber) is inserted.
Two PMTs are placed at the end of the scintillator bar. This class defines the strip design in a "object-oriented" way, so that it can be allocated by the detector construction. The allocation of the strip module class in done in the plane module class.
* @author S. Riggi
* @date 05/04/2010
*/
#ifndef _G4MuonCounterSimulatorUSC_G4StripModule_h_
#define _G4MuonCounterSimulatorUSC_G4StripModule_h_ 1

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

class G4StripModule : public G4PVPlacement
{
public:

  /** 
	\brief Class constructor: create the structure of a strip module
 	*/
  G4StripModule(G4RotationMatrix *pRot,
		const G4ThreeVector &tlate,
		G4LogicalVolume *pMotherLogical,
		bool pMany,
		int pCopyNo,
		G4MuonCounterConstruction* c);

  /**
	* \brief ### TO BE REMOVED ###
	*/
	void SetHousingVolume(G4LogicalVolume* volume_log){ModuleHousing_log= volume_log;}


private:

  /**
	* \brief Set the visualization attributes of the strip (colours, visibility, etc...)
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
	int fFiberOptCouplingMode;

  double fScintSizeX;
  double fScintSizeY;
  double fScintSizeZ;

	double fScintCoatingSizeX;
  double fScintCoatingSizeY;
  double fScintCoatingSizeZ;
	double fScintCoatingThickness;
	double fScintCoatingReflectivity;

	double fScintCoatingCapSizeX;
  double fScintCoatingCapSizeY;
  double fScintCoatingCapSizeZ;

	double fScintEnergyThreshold;

	bool fUseFiber;
	bool fUseDoubleReadout;

  double fGrooveSizeX;
  double fGrooveSizeY;
  double fGrooveSizeZ;

	double fFiberRmin;
  double fFiberRmax;
	double fFiberSizeZ;
	double fFiberPhiMin;
	double fFiberPhiMax;

	double fClad1Rmin;
  double fClad1Rmax;
	double fClad1SizeZ;
	double fClad1PhiMin;
	double fClad1PhiMax;
  double fClad2Rmin;
  double fClad2Rmax;
	double fClad2SizeZ;
	double fClad2PhiMin;
	double fClad2PhiMax;
	
	double fPMTSizeX;
  double fPMTSizeY;
	double fPMTSizeZ;

	double fPMTOptCouplingSizeX;
  double fPMTOptCouplingSizeY;
	double fPMTOptCouplingSizeZ;

	double fPhotocathodeSizeX;
  double fPhotocathodeSizeY;
  double fPhotocathodeSizeZ;

  bool fUseRealReflectivity;
	bool fUseRealPhotocathode;
	int fCoatingSurfaceType;

	double fSurfaceSeparationGap;

  //Basic Volumes
  //
	G4Box* ModuleHousing_box;
  G4Box* Scint_box;
  G4Box* CoatSide_box;
	G4Box* CoatUp_box;
	G4Box* CoatBottom_box;
	G4Box* Groove_box;
	
  G4Box* PMT_box;
	G4Box* Photocathode_box;	
	G4Box* PMTOptCoupling_box;
  G4Tubs* Fiber_tube;
	G4Tubs* FiberClad1_tube;
	G4Tubs* FiberClad2_tube;
  

  // Logical volumes
  //
  G4LogicalVolume* Scint_log;	
	G4LogicalVolume* CoatUp_log;
	G4LogicalVolume* CoatBottom_log;
	G4LogicalVolume* CoatSide_log;
	G4LogicalVolume* Groove_log;
	G4LogicalVolume* Fiber_log;
	G4LogicalVolume* FiberClad1_log;
	G4LogicalVolume* FiberClad2_log;	
	G4LogicalVolume* PMT_log;
  G4LogicalVolume* Photocathode_log;
	G4LogicalVolume* PMTOptCoupling_log;

  static G4LogicalVolume* ModuleHousing_log;
  
  
  // Physical volumes
  //
  G4VPhysicalVolume* Scint_phys;
	G4VPhysicalVolume* CoatUp_phys;
	G4VPhysicalVolume* CoatSide_phys;
	G4VPhysicalVolume* CoatBottom_phys;
	G4VPhysicalVolume* Groove_phys;	
	G4VPhysicalVolume* Fiber_phys;
	G4VPhysicalVolume* FiberClad1_phys;
  G4VPhysicalVolume* FiberClad2_phys;
  G4VPhysicalVolume* PMT_phys;
  G4VPhysicalVolume* Photocathode_phys;
	G4VPhysicalVolume* PMTOptCoupling_phys;

  //Sensitive Detectors
  static G4ScintillatorSD* Scint_SD;
  static G4PMTSD* PMT_SD;
	//static G4MultiFunctionalDetector* MFDet;
};

}//close namespace
#endif
