/**
* @file G4PlaneModule.cc
* @class G4PlaneModule
* @brief Define a detector plane
*
* The detector plane module is made up of a series of strip modules, one aside to the other. This class defines the plane design in a "object-oriented" way, so that it can be allocated by the detector construction. The allocation of the strip module class in done in this class, while the allocation of the detector planes is done in the detector costruction class.
* @author S. Riggi
* @date 05/04/2010
*/


#include "G4PlaneModule.hh"
#include "G4StripModule.hh"
#include "G4SuperPlaneModule.hh"
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

#include "G4PlaneModulePSFlatSurfaceCurrent.hh"

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

G4LogicalVolume* G4PlaneModule::ModuleHousing_log=NULL;
G4MultiFunctionalDetector* G4PlaneModule::MFDet;

G4PlaneModule::G4PlaneModule(G4RotationMatrix *pRot,
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
		            "PlaneHousing_phys",pMotherLogical,pMany,pCopyNo),fDetectorConstructor(c)
{
  CopyValues();

  if(!ModuleHousing_log || fUpdated){
 	
		//#########################
		//##  MODULE HOUSING
		//#########################
		double fHousingSizeX;
		double fHousingSizeY;
		double fHousingSizeZ;
		//double fHousingSizeX= fPlaneSizeX;
    //double fHousingSizeY= fPlaneSizeY;
    //double fHousingSizeZ= fPlaneSizeZ;

		//## Set module housing size
		if(fPMTSizeY>fScintSizeY+2.*fScintCoatingThickness || fPMTSizeZ>fScintSizeZ+2.*fScintCoatingThickness){
			cout<<"G4PlaneModule(): ERROR: pmt size larger than strip coating size! This would create a geometry overlap...exit!"<<endl;
			exit(1);
		}
		fHousingSizeY= fNstrips*fScintCoatingSizeY;
    fHousingSizeZ= fScintCoatingSizeZ;
		

		//accomodate fiber?
		if(fStripDesignMode == G4MuonCounterConstruction::kScintPlusFiber) 
			fHousingSizeX= max(fScintSizeX+2.*fScintCoatingThickness, fFiberSizeZ);
		else
			fHousingSizeX= fScintSizeX + 2.*fScintCoatingThickness;

		//accomodate pmt?
		fHousingSizeX+= 2.*fPMTSizeX;		

		//accomodate optical coupling?
		if(fStripOptCouplingMode != G4MuonCounterConstruction::kNoCoupling)
			fHousingSizeX+= 2.*fPMTOptCouplingSizeX;	

		//add small gap to avoid surface sharing?
		fHousingSizeZ+= 2.*fSurfaceSeparationGap;//strip module gap
		fHousingSizeZ+= 2.*fSurfaceSeparationGap;//plane module gap


		ModuleHousing_box = new G4Box("PlaneHousing_box",fHousingSizeX/2.,fHousingSizeY/2.,fHousingSizeZ/2.);
		ModuleHousing_log = new G4LogicalVolume(ModuleHousing_box,G4Material::GetMaterial("Air"),"PlaneHousing_log",0,0,0);
		
		//## issue: which material for the housing? air/aluminium?
 
		//###################
		//##  PLANE MODULE
		//###################
		double offset= -(fPlaneSizeY/2.-fScintCoatingSizeY/2.);
		int copyNo= 0;
		for(int i=0;i<fNstrips;i++){
			new G4StripModule(0,G4ThreeVector(0,offset,0),ModuleHousing_log,false,copyNo,fDetectorConstructor);
      offset+= fScintCoatingSizeY;
			copyNo++;
		}


		//###############################
		//## SENSITIVE VOLUMES & SCORERS
		//###############################
		G4SDManager* SDman = G4SDManager::GetSDMpointer();

	
  	//## Define MultiFunctionalDetector with name
		if(!MFDet){
  		MFDet = new G4MultiFunctionalDetector("PlaneModuleSD");
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
		G4PlaneModulePSFlatSurfaceCurrent* muonSurfaceInCurrentScorer = new G4PlaneModulePSFlatSurfaceCurrent("muonSurfaceInCurrent", fCurrent_In);
		muonSurfaceInCurrentScorer->SetFilter(muonFilter);
		muonSurfaceInCurrentScorer->DivideByArea(false);

		G4PlaneModulePSFlatSurfaceCurrent* muonSurfaceOutCurrentScorer = new G4PlaneModulePSFlatSurfaceCurrent("muonSurfaceOutCurrent", fCurrent_Out);
		muonSurfaceOutCurrentScorer->SetFilter(muonFilter);
		muonSurfaceOutCurrentScorer->DivideByArea(false);

		G4PlaneModulePSFlatSurfaceCurrent* emSurfaceInCurrentScorer = new G4PlaneModulePSFlatSurfaceCurrent("emSurfaceInCurrent", fCurrent_In);
		emSurfaceInCurrentScorer->SetFilter(emFilter);
		emSurfaceInCurrentScorer->DivideByArea(false);
  	
		G4PlaneModulePSFlatSurfaceCurrent* emSurfaceOutCurrentScorer = new G4PlaneModulePSFlatSurfaceCurrent("emSurfaceOutCurrent", fCurrent_Out);
		emSurfaceOutCurrentScorer->SetFilter(emFilter);
		emSurfaceOutCurrentScorer->DivideByArea(false);
  	

  	//## Register primitive scorers to MultiFunctionalDetector
  	MFDet->RegisterPrimitive(muonSurfaceInCurrentScorer);
		MFDet->RegisterPrimitive(emSurfaceInCurrentScorer);
  	MFDet->RegisterPrimitive(muonSurfaceOutCurrentScorer);
		MFDet->RegisterPrimitive(emSurfaceOutCurrentScorer);
  	
		
		//###############################
		//## VIS ATTRIBUTES & SURFACES
		//###############################
    VisAttributes();
    SurfaceProperties();
  }//close if

  SetLogicalVolume(ModuleHousing_log);

}//close function


void G4PlaneModule::CopyValues(){

  fUpdated= fDetectorConstructor->GetUpdated();

	fStripDesignMode= fDetectorConstructor->GetStripDesignMode();
	fStripReadoutMode= fDetectorConstructor->GetStripReadoutMode();
	fStripOptCouplingMode= fDetectorConstructor->GetStripOptCouplingMode();

	fNstrips= fDetectorConstructor->GetNumberOfStripsInPlane();

  fScintSizeX= (fDetectorConstructor->GetScintSize()).x();
  fScintSizeY= (fDetectorConstructor->GetScintSize()).y();
  fScintSizeZ= (fDetectorConstructor->GetScintSize()).z();

	fFiberSizeZ= fDetectorConstructor->GetFiberSizeZ();

	fScintCoatingThickness= fDetectorConstructor->GetHousingThickness();
  fScintCoatingSizeX= (fDetectorConstructor->GetScintCoatingSize()).x();
	fScintCoatingSizeY= (fDetectorConstructor->GetScintCoatingSize()).y();
	fScintCoatingSizeZ= (fDetectorConstructor->GetScintCoatingSize()).z();
	fScintCoatingReflectivity= fDetectorConstructor->GetHousingReflectivity();

	fPMTOptCouplingSizeX= (fDetectorConstructor->GetPMTOpticalCouplingSize()).x();
  fPMTOptCouplingSizeY= (fDetectorConstructor->GetPMTOpticalCouplingSize()).y();
  fPMTOptCouplingSizeZ= (fDetectorConstructor->GetPMTOpticalCouplingSize()).z();

	fPMTSizeX= (fDetectorConstructor->GetPMTSize()).x();
  fPMTSizeY= (fDetectorConstructor->GetPMTSize()).y();
  fPMTSizeZ= (fDetectorConstructor->GetPMTSize()).z();

	fSurfaceSeparationGap= fDetectorConstructor->GetSurfaceSeparationGap();

	if(fStripOptCouplingMode==G4MuonCounterConstruction::kNoCoupling) 
		fPlaneSizeX= max(fScintSizeX+2.*fScintCoatingThickness, fFiberSizeZ) + 2.*fPMTSizeX;
	else
		fPlaneSizeX= max(fScintSizeX+2.*fScintCoatingThickness, fFiberSizeZ) + 2.*fPMTOptCouplingSizeX + 2.*fPMTSizeX;

	fPlaneSizeY= fNstrips*fScintCoatingSizeY;
  fPlaneSizeZ= fScintCoatingSizeZ;

}//close function


void G4PlaneModule::VisAttributes(){
	
  G4VisAttributes* ModuleHousing_va = new G4VisAttributes(G4Colour(0.,1.,0.));//green
  ModuleHousing_log->SetVisAttributes(G4VisAttributes::Invisible);
  //ModuleHousing_log->SetVisAttributes(ModuleHousing_va);	

}//close function


void G4PlaneModule::SurfaceProperties(){ 

}//close function



