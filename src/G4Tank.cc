/**
* @file G4Tank.hh
* @class G4Tank
* @brief Define a tank
*
* @author S. Riggi
* @date 15/04/2010
*/

#include "G4Tank.hh"

#include "G4Ellipsoid.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

#include <algorithm>

using namespace std;
using namespace G4MuonCounterSimulatorUSC;

G4LogicalVolume* G4Tank::ModuleHousing_log=NULL;

G4Tank::G4Tank(G4RotationMatrix *pRot,
			                   const G4ThreeVector &tlate,
			                   G4LogicalVolume *pMotherLogical,
			                   bool pMany,
			                   int pCopyNo,
			                   G4MuonCounterConstruction* c)
  //Pass info to the G4PVPlacement constructor
  :G4PVPlacement(pRot,tlate,
		            //Temp logical volume must be created here
		            new G4LogicalVolume(new G4Box("temp",1,1,1),
				        G4Material::GetMaterial("Vacuum"),
				        "temp",0,0,0),
		            "Tank_phys",pMotherLogical,pMany,pCopyNo),fDetectorConstructor(c)
{
  CopyValues();

  if(!ModuleHousing_log || fUpdated){
 	
		//#########################
		//##  MODULE HOUSING
		//#########################
		ModuleHousing_tube = new G4Tubs("ModuleHousing_tube", 0.0*CLHEP::m, fTankRadius + fTankThickness, fTankHalfHeight, 0.0*CLHEP::deg, 360.0*CLHEP::deg);
		ModuleHousing_log = new G4LogicalVolume(ModuleHousing_tube, G4Material::GetMaterial("Vacuum"), "ModuleHousing_log", 0, 0, 0);
  		

		//## Tank
		tank_solid = new G4Tubs("tank_solid", 0.0*CLHEP::m, fTankRadius, fTankHalfHeight, 0.0*CLHEP::deg, 360.0*CLHEP::deg);
		tank_log = new G4LogicalVolume(tank_solid, G4Material::GetMaterial("Water"), "tank_log", 0, 0, 0);
  	tank_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),tank_log,"tank_phys",ModuleHousing_log, false, 0);
  
		//## Tank top, side and bottom
  	top_solid = new G4Tubs("top_solid", 0.0*CLHEP::m, fTankRadius + fTankThickness, fTankThickness/2,0.0*CLHEP::deg, 360.0*CLHEP::deg);
		top_log = new G4LogicalVolume(top_solid, G4Material::GetMaterial("HDPE"), "top_log", 0, 0, 0);
		top_phys = new G4PVPlacement(0, G4ThreeVector(0,0, fTankHalfHeight + fTankThickness/2), top_log,"top_phys", ModuleHousing_log, false, 0);



  	side_solid = new G4Tubs("side_solid", fTankRadius, fTankRadius + fTankThickness, fTankHalfHeight,0.0*CLHEP::deg, 360.0*CLHEP::deg);
		side_log = new G4LogicalVolume(side_solid, G4Material::GetMaterial("HDPE"), "side_log", 0, 0, 0);
  	side_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),side_log,"side_phys", ModuleHousing_log, false, 0);


		bottom_log = new G4LogicalVolume(top_solid, G4Material::GetMaterial("HDPE"), "bottom_log", 0, 0, 0);
  	bottom_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTankHalfHeight+fTankThickness/2)),bottom_log,"bottom_phys",  ModuleHousing_log, false, 0);




		//## Inner		
		inner_solid = new G4Ellipsoid("inner_solid", fPmtRmin, fPmtRmin, fPmtRzmin, -fPmtRzmin, 0.);
		inner_log = new G4LogicalVolume(inner_solid, G4Material::GetMaterial("Vacuum"), "inner_log", 0, 0, 0);
  
		inner1_phys = new G4PVPlacement(0, G4ThreeVector(dome1_x, dome1_y, dome1_z), inner_log,
			"inner1_phys", tank_log, false, 0);

  	inner2_phys = new G4PVPlacement(0, G4ThreeVector(dome2_x, dome2_y, dome2_z), inner_log,
			"inner2_phys", tank_log, false, 0);
  
		inner3_phys = new G4PVPlacement(0, G4ThreeVector(dome3_x, dome3_y, dome3_z), inner_log,
			"inner3_phys", tank_log, false, 0);

		//## PMTs	
		pmt_aux = new G4Ellipsoid("pmt_aux", fPmtRmax, fPmtRmax, fPmtRzmax, -fPmtRzmax, -fHeightz);

		pmt_solid = new G4SubtractionSolid("pmt_solid", pmt_aux, inner_solid);
		pmt_aux1 = new G4Ellipsoid("pmt_aux1", fPmtRmax, fPmtRmax, fPmtRzmax, -fHeightz, 0.);
		pmt_solid1 = new G4SubtractionSolid("pmt_solid1", pmt_aux1, inner_solid);
		
		pmt1_log = new G4LogicalVolume(pmt_solid, G4Material::GetMaterial("Pyrex"), "pmt1_log", 0, 0, 0);
 		pmt2_log = new G4LogicalVolume(pmt_solid, G4Material::GetMaterial("Pyrex"), "pmt2_log", 0, 0, 0);
  	pmt3_log = new G4LogicalVolume(pmt_solid, G4Material::GetMaterial("Pyrex"), "pmt3_log", 0, 0, 0);
  
  	pmt1_phys = new G4PVPlacement(0, G4ThreeVector(dome1_x, dome1_y, dome1_z), pmt1_log,"pmt1_phys", tank_log, false, 0);
  	pmt2_phys = new G4PVPlacement(0, G4ThreeVector(dome2_x, dome2_y, dome2_z), pmt2_log,"pmt2_phys", tank_log, false, 0);
  	pmt3_phys = new G4PVPlacement(0, G4ThreeVector(dome3_x, dome3_y, dome3_z), pmt3_log,"pmt3_phys", tank_log, false, 0);


		//## The dome/face interface
   	interface_in_aux = new G4Ellipsoid("interface_in_aux", fInterfaceRmin, fInterfaceRmin,
				fInterfaceRzmin,-fInterfaceRzmin,0.);
  	interface_out_aux = new G4Ellipsoid("interface_out_aux", fInterfaceRmax, fInterfaceRmax,	
		    fInterfaceRzmax, -fInterfaceRzmax, 0.0*CLHEP::cm);
  	interface_solid = new G4SubtractionSolid("interface_solid", interface_out_aux, interface_in_aux);

  	interface1_log = new G4LogicalVolume(interface_solid, G4Material::GetMaterial("Interface"), "interface1_log", 0, 0, 0);
  	interface2_log = new G4LogicalVolume(interface_solid, G4Material::GetMaterial("Interface"), "interface2_log", 0, 0, 0);
  	interface3_log = new G4LogicalVolume(interface_solid, G4Material::GetMaterial("Interface"), "interface3_log", 0, 0, 0);
  
  	interface1_phys = new G4PVPlacement(0, G4ThreeVector(dome1_x, dome1_y, dome1_z), 
				interface1_log,"interface1_phys", tank_log, false, 0);

  	interface2_phys = new G4PVPlacement(0, G4ThreeVector(dome2_x, dome2_y, dome2_z), 
				interface2_log,"interface2_phys", tank_log, false, 0);

  	interface3_phys = new G4PVPlacement(0, G4ThreeVector(dome3_x, dome3_y, dome3_z), 
				interface3_log,"interface3_phys", tank_log, false, 0);

		//## PMT domes 
		dome_in_aux = new G4Ellipsoid("dome_in_aux", fDomeRmin, fDomeRmin, fDomeRzmin, -fDomeRzmin, 0.);
  	dome_out_aux = new G4Ellipsoid("dome_out_aux", fDomeRmax, fDomeRmax, fDomeRzmax, -fDomeRzmax, 0.0*CLHEP::cm);
  	dome_solid = new G4SubtractionSolid("dome_solid", dome_out_aux, dome_in_aux);

 
  	dome1_log = new G4LogicalVolume(dome_solid, G4Material::GetMaterial("Lucite"), "dome1_log", 0, 0, 0);
  	dome2_log = new G4LogicalVolume(dome_solid, G4Material::GetMaterial("Lucite"), "dome2_log", 0, 0, 0);
  	dome3_log = new G4LogicalVolume(dome_solid, G4Material::GetMaterial("Lucite"), "dome3_log", 0, 0, 0);
    
  	dome1_phys = new G4PVPlacement(0, G4ThreeVector(dome1_x, dome1_y, dome1_z), 
				dome1_log, "dome1_phys", tank_log, false, 0);
  	dome2_phys = new G4PVPlacement(0, G4ThreeVector(dome2_x, dome2_y, dome2_z), 
				dome2_log, "dome2_phys", tank_log, false, 0);
  	dome3_phys = new G4PVPlacement(0, G4ThreeVector(dome3_x, dome3_y, dome3_z), 
				dome3_log, "dome3_phys", tank_log, false, 0);



    VisAttributes();
    SurfaceProperties();
  }//close if

  SetLogicalVolume(ModuleHousing_log);

}//close function


void G4Tank::CopyValues(){

  fUpdated= fDetectorConstructor->GetUpdated();

	//## Tank parameters
	scaleFactor= 10;

	fTankRadius= 1.8*CLHEP::m/scaleFactor;
  fTankHalfHeight= 0.5*1.2*CLHEP::m/scaleFactor;
  fTankThickness= 12.7*CLHEP::mm/scaleFactor;
  top_TankThickness = fTankThickness/2;

	fFaceRadius= 114.*CLHEP::mm/scaleFactor;
  fFaceRadiusz= 84.5781 *CLHEP::mm/scaleFactor;
  fFaceActiveRadius= 104*CLHEP::mm/scaleFactor;

	fInterfaceThickness= 1.*CLHEP::mm/scaleFactor;
  fGlassThickness= 2.5*CLHEP::mm/scaleFactor;
	fDomeThickness= 5.2*CLHEP::mm/scaleFactor;
 
	fPmtRmin =        fFaceRadius-fGlassThickness;
  fPmtRmax =        fFaceRadius;
  fPmtRzmin =       fFaceRadiusz-fGlassThickness;
  fPmtRzmax =       fFaceRadiusz;
  fInterfaceRmin =  fFaceRadius;
  fInterfaceRzmin = fFaceRadiusz;
  fInterfaceRmax =  fFaceRadius + fInterfaceThickness;
  fInterfaceRzmax = fFaceRadiusz + fInterfaceThickness;
  fDomeRmin =       fFaceRadius + fInterfaceThickness;
  fDomeRzmin =      fFaceRadiusz + fInterfaceThickness;
  fDomeRmax =       fFaceRadius + fInterfaceThickness + fDomeThickness;
  fDomeRzmax =      fFaceRadiusz + fInterfaceThickness + fDomeThickness;

	
	// Taking into account that the radius of the photocathode is fFaceActiveRadius
  if ( fFaceActiveRadius/CLHEP::mm >= 114./scaleFactor ) {
		cout<<"face active radius 1"<<endl;
  	fFaceActiveRadius = 113.926112814*CLHEP::mm/scaleFactor;
    fHeightz = 0.;
  }
  else if ( fFaceActiveRadius/CLHEP::mm > 97.27/scaleFactor ) {		
		cout<<"face active radius 2"<<endl;
  	alpha = acos(1./62.*(fFaceActiveRadius/CLHEP::mm-133.*sin(47.*deg)+62.*cos(43.*deg)));
    fHeightz = 62.*sin(alpha)*CLHEP::mm/scaleFactor;
  }
  else{
		cout<<"face active radius 3"<<endl;
  	alpha = asin(fFaceActiveRadius/CLHEP::mm/133.);
    fHeightz = (133.*cos(alpha)-(133.-62.)*cos(47.*deg))*CLHEP::mm/scaleFactor;
  }
	cout<<"alpha="<<alpha<<"  fHeightz="<<fHeightz/CLHEP::mm<<endl;

	
	dome1_x = 0*CLHEP::m/scaleFactor;  
  dome1_y = 1.2*CLHEP::m/scaleFactor;
  dome1_z = 0.6*CLHEP::m/scaleFactor;
  
  dome2_x = 1.0392*CLHEP::m/scaleFactor;
  dome2_y = -0.6*CLHEP::m/scaleFactor;
  dome2_z = 0.6*CLHEP::m/scaleFactor;
  
  dome3_x = -1.0392*CLHEP::m/scaleFactor;
  dome3_y = -0.6*CLHEP::m/scaleFactor;
  dome3_z = 0.6*CLHEP::m/scaleFactor;


}//close function


void G4Tank::VisAttributes(){
		
	//## Colors
	//   Yellow : 1 1 0
	//   Cyan   : 0 1 1
	//   Green  : 0 1 0

  G4VisAttributes* ModuleHousing_va = new G4VisAttributes(G4Colour(0.,0.,0.,0.8));//black
  ModuleHousing_log->SetVisAttributes(G4VisAttributes::Invisible);
  //ModuleHousing_log->SetVisAttributes(ModuleHousing_va);

}//close function


void G4Tank::SurfaceProperties(){ 

}//close function



