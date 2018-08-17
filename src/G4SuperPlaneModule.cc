/**
* @file G4SuperPlaneModule.cc
* @class G4SuperPlaneModule
* @brief Define a super plane module
*
* The superplane module contains two XY planes placed inside a casing
* @author S. Riggi, E. Trovato
* @date 15/04/2010
*/


#include "G4SuperPlaneModule.hh"
#include "G4PlaneModule.hh"
#include "G4StripModule.hh"
#include "G4PMTSD.hh"
#include "G4ScintillatorSD.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofStep.hh"
#include "G4PSCellFlux.hh"
#include "G4PSPassageCellFlux.hh"
#include "G4PSFlatSurfaceFlux.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDChargedFilter.hh"

//#include "G4SuperPlaneModulePSFlatSurfaceCurrent.hh"

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

G4LogicalVolume* G4SuperPlaneModule::ModuleHousing_log=NULL;
G4LogicalVolume* G4SuperPlaneModule::CasingInside_log=NULL;
//G4MultiFunctionalDetector* G4SuperPlaneModule::MFDet;


G4SuperPlaneModule::G4SuperPlaneModule(G4RotationMatrix *pRot,
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
		            "SuperPlaneHousing_phys",pMotherLogical,pMany,pCopyNo),fDetectorConstructor(c)
{
  CopyValues();

  if(!ModuleHousing_log || fUpdated){
 	
		//#########################
		//##  MODULE HOUSING
		//#########################
		double fHousingSizeX= fSuperPlaneSizeX;
    double fHousingSizeY= fSuperPlaneSizeY;
    double fHousingSizeZ= fSuperPlaneSizeZ;
	
		ModuleHousing_box = new G4Box("SuperPlaneHousing_box",fHousingSizeX/2.,fHousingSizeY/2.,fHousingSizeZ/2.);
		G4Material* casingMaterial= NULL;
		if(fSuperPlaneMaterialName=="PVC")
			casingMaterial= G4Material::GetMaterial("G4_POLYVINYL_CHLORIDE");
		else if(fSuperPlaneMaterialName=="GlassFiber")
			casingMaterial= G4Material::GetMaterial("G4_SILICON_DIOXIDE");
		else if(fSuperPlaneMaterialName=="Vacuum")
			casingMaterial= G4Material::GetMaterial("Vacuum");
		else if(fSuperPlaneMaterialName=="Air")
			casingMaterial= G4Material::GetMaterial("Air");
		else{
			G4cerr<<"Invalid casing material!"<<G4endl;
			exit(1);	
		}

		ModuleHousing_log = new G4LogicalVolume(ModuleHousing_box,casingMaterial,"SuperPlaneHousing_log",0,0,0);
		
		//###########################
		//##  CASING INSIDE
		//###########################
		CasingInside_box = new G4Box("CasingInside_box",fSuperPlaneInsideSizeX/2,fSuperPlaneInsideSizeY/2,fSuperPlaneInsideSizeZ/2);
		
		G4Material* insidecasingMaterial= NULL;
		if(fSuperPlaneMaterialName=="Vacuum")
			insidecasingMaterial= G4Material::GetMaterial("Vacuum");
		else
			insidecasingMaterial= G4Material::GetMaterial("Air");	
	
		CasingInside_log= new G4LogicalVolume(CasingInside_box,insidecasingMaterial,"CasingInside_log",0,0,0);		

		CasingInside_phys= new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),CasingInside_log,"CasingInside_phys",ModuleHousing_log,false,0);


		//###########################
		//##  SUPER PLANE MODULE
		//###########################
		//PLANE ROTATION MATRIX
		G4RotationMatrix* planeRotMatrix = new G4RotationMatrix;	
		planeRotMatrix->rotateZ(fPlaneTiltAngle);

		//BUILD DETECTOR PLANES
		int planeCopyNo=0;
		double planeOffset= fPlaneSizeZ/2.;
		planeOffset+= fSurfaceSeparationGap;//strip module
		planeOffset+= fSurfaceSeparationGap;//plane module
		planeOffset+= fXYPlaneDistance/2.;//x-y distance
		
		if(fHasYPlaneInSuperModule){
			//## X+Y planes
			new G4PlaneModule(0,G4ThreeVector(0.,0.,+planeOffset),CasingInside_log,false,planeCopyNo,fDetectorConstructor);			
			planeCopyNo++;
		
			new G4PlaneModule(planeRotMatrix,G4ThreeVector(0.,0.,-planeOffset),CasingInside_log,false,planeCopyNo,fDetectorConstructor);
		}
		else{
			//## X plane
			new G4PlaneModule(0,G4ThreeVector(0.,0.,0.),CasingInside_log,false,planeCopyNo,fDetectorConstructor);			
		}

	/*
		//###############################
		//## SENSITIVE VOLUMES & SCORERS
		//###############################
		G4SDManager* SDman = G4SDManager::GetSDMpointer();

		
  	//## Define MultiFunctionalDetector with name
		if(!MFDet){
  		MFDet = new G4MultiFunctionalDetector("SuperPlaneModuleSD");
  		SDman->AddNewDetector(MFDet);                 // Register SD to SDManager
		}
  	ModuleHousing_log->SetSensitiveDetector(MFDet);  // Assign SD to the logical volume.

		//## Define sensitive detector filters
  	G4String fltName,particleName;
  	//-- proton filter
  	G4SDParticleFilter* protonFilter = new G4SDParticleFilter(fltName="protonFilter", particleName="proton");
		//-- gamma filter
  	G4SDParticleFilter* gammaFilter = new G4SDParticleFilter(fltName="gammaFilter", particleName="gamma");
  	//-- electron filter
  	G4SDParticleFilter* electronFilter = new G4SDParticleFilter(fltName="electronFilter");
  	electronFilter->add(particleName="e+");   // accept electrons
  	electronFilter->add(particleName="e-");   // accept positrons
		//-- muon filter		
		G4SDParticleFilter* muonFilter = new G4SDParticleFilter(fltName="muonFilter");
  	muonFilter->add(particleName="mu+");   // accept mu+
  	muonFilter->add(particleName="mu-");   // accept mu-
		//-- em filter
  	G4SDParticleFilter* emFilter = new G4SDParticleFilter(fltName="emFilter");
  	emFilter->add(particleName="e+");   // accept electrons
  	emFilter->add(particleName="e-");   // accept positrons
		emFilter->add(particleName="gamma");   // accept gamma
  	//-- charged particle filter
 	 	G4SDChargedFilter* chargedFilter = new G4SDChargedFilter(fltName="chargedFilter");

  	//## Define primitive scorers 
		//G4SuperPlaneModulePSFlatSurfaceCurrent* muonSurfaceCurrentScorer = new G4SuperPlaneModulePSFlatSurfaceCurrent("muonSurfaceCurrent", fCurrent_Out);
		G4PSFlatSurfaceCurrent* muonSurfaceCurrentScorer = new G4PSFlatSurfaceCurrent("muonSurfaceCurrent", fCurrent_Out);
		muonSurfaceCurrentScorer->SetFilter(muonFilter);
		muonSurfaceCurrentScorer->DivideByArea(false);

		//G4SuperPlaneModulePSFlatSurfaceCurrent* emSurfaceCurrentScorer = new G4SuperPlaneModulePSFlatSurfaceCurrent("emSurfaceCurrent", fCurrent_Out);
		G4PSFlatSurfaceCurrent* emSurfaceCurrentScorer = new G4PSFlatSurfaceCurrent("emSurfaceCurrent", fCurrent_Out);
		emSurfaceCurrentScorer->SetFilter(emFilter);
		emSurfaceCurrentScorer->DivideByArea(false);
  	
  	//## Register primitive scorers to MultiFunctionalDetector
  	MFDet->RegisterPrimitive(muonSurfaceCurrentScorer);
		MFDet->RegisterPrimitive(emSurfaceCurrentScorer);
	*/

    VisAttributes();
    SurfaceProperties();
  }//close if

  SetLogicalVolume(ModuleHousing_log);

}//close function


void G4SuperPlaneModule::CopyValues(){

  fUpdated= fDetectorConstructor->GetUpdated();

	fStripDesignMode= fDetectorConstructor->GetStripDesignMode();
	fStripReadoutMode= fDetectorConstructor->GetStripReadoutMode();
	fStripOptCouplingMode= fDetectorConstructor->GetStripOptCouplingMode();

	fNstrips= fDetectorConstructor->GetNumberOfStripsInPlane();

  fScintSizeX= (fDetectorConstructor->GetScintSize()).x();
  fScintSizeY= (fDetectorConstructor->GetScintSize()).y();
  fScintSizeZ= (fDetectorConstructor->GetScintSize()).z();

	fScintCoatingThickness= fDetectorConstructor->GetHousingThickness();
  fScintCoatingSizeX= (fDetectorConstructor->GetScintCoatingSize()).x();
	fScintCoatingSizeY= (fDetectorConstructor->GetScintCoatingSize()).y();
	fScintCoatingSizeZ= (fDetectorConstructor->GetScintCoatingSize()).z();
	fScintCoatingReflectivity= fDetectorConstructor->GetHousingReflectivity();

	fFiberSizeZ= fDetectorConstructor->GetFiberSizeZ();
	
	fPMTSizeX= (fDetectorConstructor->GetPMTSize()).x();
  fPMTSizeY= (fDetectorConstructor->GetPMTSize()).y();
  fPMTSizeZ= (fDetectorConstructor->GetPMTSize()).z();

	fPMTOptCouplingSizeX= (fDetectorConstructor->GetPMTOpticalCouplingSize()).x();
  fPMTOptCouplingSizeY= (fDetectorConstructor->GetPMTOpticalCouplingSize()).y();
  fPMTOptCouplingSizeZ= (fDetectorConstructor->GetPMTOpticalCouplingSize()).z();

	if(fStripOptCouplingMode==G4MuonCounterConstruction::kNoCoupling) 
		fPlaneSizeX= max(fScintSizeX+2.*fScintCoatingThickness, fFiberSizeZ) + 2.*fPMTSizeX;
	else
		fPlaneSizeX= max(fScintSizeX+2.*fScintCoatingThickness, fFiberSizeZ) + 2.*fPMTOptCouplingSizeX + 2.*fPMTSizeX;
  
	fPlaneSizeY= fNstrips*fScintCoatingSizeY;
  fPlaneSizeZ= fScintCoatingSizeZ;

	fXYPlaneDistance= fDetectorConstructor->GetDistanceAmongXYPlanes();
	fPlaneTiltAngle= fDetectorConstructor->GetPlaneTiltAngle();

	fSuperPlaneThickness= fDetectorConstructor->GetSuperPlaneModuleThickness();
	
	fSuperPlaneMaterialName= fDetectorConstructor->GetSuperPlaneMaterialName();
	fHasYPlaneInSuperModule= fDetectorConstructor->IsYPlaneInSuperModule();

	fSurfaceSeparationGap= fDetectorConstructor->GetSurfaceSeparationGap();

	//## set superplane size
	fSuperPlaneSizeX= fPlaneSizeX + 2.0*fSuperPlaneThickness;
  fSuperPlaneSizeY= fPlaneSizeX + 2.0*fSuperPlaneThickness;

	//check if given superplane size Z is enough to accomodate all the detector
	double minSuperPlaneSizeZ;
	minSuperPlaneSizeZ= fPlaneSizeZ;
	minSuperPlaneSizeZ+= 2.*fSurfaceSeparationGap;//strip module gap
	minSuperPlaneSizeZ+= 2.*fSurfaceSeparationGap;//plane module module gap
	
	//if(fHasYPlaneInSuperModule) minSuperPlaneSizeZ= 2.0*fPlaneSizeZ+ 2.0*fSuperPlaneThickness + fXYPlaneDistance;
	//else minSuperPlaneSizeZ= fPlaneSizeZ+ 2.0*fSuperPlaneThickness;
	/*	
	if(fHasYPlaneInSuperModule) {
		minSuperPlaneSizeZ= 2.0*fPlaneSizeZ+ 2.0*fSuperPlaneThickness + fXYPlaneDistance;
		minSuperPlaneSizeZ+= 2.*fSurfaceSeparationGap;//strip module gap
		minSuperPlaneSizeZ+= 4.*fSurfaceSeparationGap;//2 plane module module gap
	}
	else {
		minSuperPlaneSizeZ= fPlaneSizeZ+ 2.0*fSuperPlaneThickness;
		minSuperPlaneSizeZ+= 2.*fSurfaceSeparationGap;//strip module gap
		minSuperPlaneSizeZ+= 2.*fSurfaceSeparationGap;//plane module module gap
	}
	*/

	
	//accomodate Y plane?
	if(fHasYPlaneInSuperModule){
		minSuperPlaneSizeZ*= 2.;
		minSuperPlaneSizeZ+= fXYPlaneDistance;	
	}	
	//accomodate superplane thickness?
	minSuperPlaneSizeZ+= 2.*fSuperPlaneThickness;
	
	
	if(fDetectorConstructor->GetSuperPlaneModuleSizeZ()<minSuperPlaneSizeZ){
		cout<<"## WARNING ##: given superplane sizeZ too small...using minimum size!"<<endl;
  	fSuperPlaneSizeZ= minSuperPlaneSizeZ;
	}
	else{
		fSuperPlaneSizeZ= fDetectorConstructor->GetSuperPlaneModuleSizeZ();
	}
	
	//## set inside casing size
	fSuperPlaneInsideSizeX= fPlaneSizeX;
  fSuperPlaneInsideSizeY= fPlaneSizeX;
	fSuperPlaneInsideSizeZ= fSuperPlaneSizeZ - 2.*fSuperPlaneThickness;

}//close function


void G4SuperPlaneModule::VisAttributes(){
		
	//## Colors
	//   Yellow : 1 1 0
	//   Cyan   : 0 1 1
	//   Green  : 0 1 0

  G4VisAttributes* ModuleHousing_va = new G4VisAttributes(G4Colour(0.,0.,0.,0.3));//black
  //ModuleHousing_log->SetVisAttributes(G4VisAttributes::Invisible);
  ModuleHousing_log->SetVisAttributes(ModuleHousing_va);

	G4VisAttributes* CasingInside_va = new G4VisAttributes(G4Colour(0., 1., 1.));//cyan
  CasingInside_log->SetVisAttributes(G4VisAttributes::Invisible);
  //CasingInside_log->SetVisAttributes(CasingInside_va);	

}//close function


void G4SuperPlaneModule::SurfaceProperties(){ 

}//close function



