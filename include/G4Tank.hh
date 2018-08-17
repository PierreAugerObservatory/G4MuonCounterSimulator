/**
* @file G4Tank.hh
* @class G4Tank
* @brief Define a tank
*
* @author S. Riggi
* @date 15/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4Tank_h_
#define _G4MuonCounterSimulatorUSC_G4Tank_h_ 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "G4MuonCounterConstruction.hh"

namespace G4MuonCounterSimulatorUSC {


class G4Tank : public G4PVPlacement
{
public:

	/** 
	\brief Class constructor: create the detector structure of a plane
 	*/
  G4Tank(G4RotationMatrix *pRot,
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

	


  //Basic Volumes
  //
	G4Tubs* ModuleHousing_tube;
	G4Tubs* tank_solid;
	G4Tubs* top_solid;
	G4Tubs* side_solid;
  G4Ellipsoid* inner_solid;
	G4Ellipsoid* pmt_aux;	
	G4SubtractionSolid* pmt_solid;	
	G4Ellipsoid* pmt_aux1;
	G4SubtractionSolid* pmt_solid1;
	G4Ellipsoid* interface_in_aux;
	G4Ellipsoid* interface_out_aux;
	G4SubtractionSolid* interface_solid;
	G4Ellipsoid* dome_in_aux;
	G4Ellipsoid* dome_out_aux;
	G4SubtractionSolid* dome_solid;


  // Logical volumes
  //
  static G4LogicalVolume* ModuleHousing_log;
	G4LogicalVolume* tank_log;
	G4LogicalVolume* top_log;
	G4LogicalVolume* side_log;
	G4LogicalVolume* bottom_log;
	G4LogicalVolume* inner_log;
	G4LogicalVolume* pmt1_log;
	G4LogicalVolume* pmt2_log;
	G4LogicalVolume* pmt3_log;
	G4LogicalVolume* interface1_log;
	G4LogicalVolume* interface2_log;
	G4LogicalVolume* interface3_log;
	G4LogicalVolume* dome1_log;
	G4LogicalVolume* dome2_log;
	G4LogicalVolume* dome3_log;

  // Physical volumes
  //
	G4VPhysicalVolume* tank_phys;
	G4VPhysicalVolume* top_phys;
	G4VPhysicalVolume* side_phys;
	G4VPhysicalVolume* bottom_phys;
	G4VPhysicalVolume* inner1_phys;
	G4VPhysicalVolume* inner2_phys;
	G4VPhysicalVolume* inner3_phys;
	G4VPhysicalVolume* pmt1_phys;
	G4VPhysicalVolume* pmt2_phys;
	G4VPhysicalVolume* pmt3_phys;
	G4VPhysicalVolume* interface1_phys;
	G4VPhysicalVolume* interface2_phys;
	G4VPhysicalVolume* interface3_phys;
	G4VPhysicalVolume* dome1_phys;
	G4VPhysicalVolume* dome2_phys;
	G4VPhysicalVolume* dome3_phys;

  double scaleFactor;
	double fTankRadius;
	double fTankHalfHeight;
	double fTankThickness;
	double fFaceRadius;
	double fFaceRadiusz;
	double fFaceActiveRadius;
	double fInterfaceThickness;
	double fGlassThickness;
 	double fDomeThickness;
	double fPmtRmin;
	double fPmtRmax;
	double fPmtRzmin;
	double fPmtRzmax;
	double fInterfaceRmin;
	double fInterfaceRzmin;	
	double fInterfaceRmax;
	double fInterfaceRzmax;
	double fDomeRmin;
	double fDomeRzmin;
	double fDomeRmax;
	double fDomeRzmax;
	double alpha;
	double fHeightz;
	double top_TankThickness;
	double dome1_x;
	double dome1_y;
	double dome1_z;
	double dome2_x;
	double dome2_y;
	double dome2_z;
	double dome3_x;
	double dome3_y;
	double dome3_z;
	

};

}//close namespace

#endif
