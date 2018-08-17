/**
* @file G4StripModule.cc
* @class G4StripModule
* @brief Define a strip module
*
* The strip module is made up of a scintillator bar, surrounded by a reflective coating. A groove is done in the bar to host a WLS fiber. Two PMTs are placed at the end of the scintillator bar. This class defines the strip design in a "object-oriented" way, so that it can be allocated by the detector construction. The allocation of the strip module class in done in the plane module class.
* @author S. Riggi
* @date 05/04/2010
*/

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

#include "G4StripModulePSFlatSurfaceCurrent.hh"


#include "globals.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"


#include <algorithm>

using namespace std;
using namespace G4MuonCounterSimulatorUSC;

G4ScintillatorSD* G4StripModule::Scint_SD;
G4PMTSD* G4StripModule::PMT_SD;
G4LogicalVolume* G4StripModule::ModuleHousing_log=NULL;
//G4MultiFunctionalDetector* G4StripModule::MFDet;


G4StripModule::G4StripModule(G4RotationMatrix *pRot,
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
		            "StripHousing_phys",pMotherLogical,pMany,pCopyNo),fDetectorConstructor(c)
{
  CopyValues();

  if(!ModuleHousing_log || fUpdated){
 	

		//#########################
		//##  MODULE HOUSING
		//#########################
		double fHousingSizeX;
		double fHousingSizeY;
		double fHousingSizeZ;

		/*
		//set module housing size
		if(fStripOptCouplingMode==G4MuonCounterConstruction::kNoCoupling) 
			fHousingSizeX= max(fScintSizeX+2.*fScintCoatingThickness, fFiberSizeZ) + 2.*fPMTSizeX;
		else
			fHousingSizeX= max(fScintSizeX+2.*fScintCoatingThickness, fFiberSizeZ) + 2.*fPMTOptCouplingSizeX + 2.*fPMTSizeX;

		if(fPMTSizeY<fScintSizeY+2.*fScintCoatingThickness || fPMTSizeZ<fScintSizeZ+2.*fScintCoatingThickness){
			cout<<"### WARNING ### : pmt size smaller than strip coating size...sizing strip housing to the coating"<<endl;
		}
		*/

		//## Set module housing size
		if(fPMTSizeY>fScintSizeY+2.*fScintCoatingThickness || fPMTSizeZ>fScintSizeZ+2.*fScintCoatingThickness){
			cout<<"G4StripModule(): ERROR: pmt size larger than strip coating size! This would create a geometry overlap...exit!"<<endl;
			exit(1);
		}
		//fHousingSizeY= max(fScintSizeY + 2.*fScintCoatingThickness, fPMTSizeY);
    //fHousingSizeZ= max(fScintSizeZ + 2.*fScintCoatingThickness, fPMTSizeZ);
		fHousingSizeY= fScintSizeY + 2.*fScintCoatingThickness;
    fHousingSizeZ= fScintSizeZ + 2.*fScintCoatingThickness;

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
		fHousingSizeZ+= 2.*fSurfaceSeparationGap;


		//create solid & logical volume
		ModuleHousing_box = new G4Box("StripHousing_box",fHousingSizeX/2.,fHousingSizeY/2.,fHousingSizeZ/2.);
		ModuleHousing_log = new G4LogicalVolume(ModuleHousing_box,G4Material::GetMaterial("Air"),"StripHousing_log",0,0,0);


		//#########################
		//##  SCINT BAR + GROOVE
		//#########################
		Scint_box = new G4Box("Scint_box",fScintSizeX/2,fScintSizeY/2,fScintSizeZ/2);
		G4Box* Hole_box = new G4Box("Hole_box",fGrooveSizeX,fGrooveSizeY/2,fGrooveSizeZ/2);

		//scint bar+groove solid
		double offset_hole= fScintSizeZ/2-fGrooveSizeZ/2;
		G4SubtractionSolid* ScintWithGroove = new G4SubtractionSolid("ScintWithGroove", Scint_box, Hole_box, 0, G4ThreeVector(0,0,offset_hole));

		if(fStripDesignMode==G4MuonCounterConstruction::kScint) {
			//no groove
			//Scint_log= new G4LogicalVolume(Scint_box,G4Material::GetMaterial("PolyvinylToluene"),"Scint_log",0,0,0);	
			Scint_log= new G4LogicalVolume(Scint_box,G4Material::GetMaterial("Polystyrene_doped"),"Scint_log",0,0,0);	
		}
		else{
			//with groove
			//Scint_log= new G4LogicalVolume(ScintWithGroove,G4Material::GetMaterial("PolyvinylToluene"),"Scint_log",0,0,0);
			Scint_log= new G4LogicalVolume(ScintWithGroove,G4Material::GetMaterial("Polystyrene_doped"),"Scint_log",0,0,0);
		}
		Scint_phys= new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),Scint_log,"Scint_phys",ModuleHousing_log,false,0);

	
		//#########################
		//##  GROOVE ENVELOP
		//#########################
		Groove_box = new G4Box("Groove_box",fGrooveSizeX/2,fGrooveSizeY/2,fGrooveSizeZ/2);

	
		G4Tubs* FiberHole_tube = new G4Tubs("FiberHole_tube",fClad2Rmin,fClad2Rmax,fClad2SizeZ,fClad2PhiMin,fClad2PhiMax);
		G4RotationMatrix* FiberRot = new G4RotationMatrix;
  	FiberRot->rotateY(M_PI/2.*CLHEP::rad); 
 
  	G4SubtractionSolid* GrooveWithFiberHole = new G4SubtractionSolid("GrooveWithFiberHole", Groove_box, FiberHole_tube, FiberRot, G4ThreeVector(0,0,0));	

		G4Material* GrooveEnvelopMaterial=NULL;

		//set the envelop material
		if(fFiberOptCouplingMode==G4MuonCounterConstruction::kNoCoupling) 
			GrooveEnvelopMaterial= G4Material::GetMaterial("Vacuum");
		else if(fFiberOptCouplingMode==G4MuonCounterConstruction::kAirCoupling)
			GrooveEnvelopMaterial= G4Material::GetMaterial("Air");		
		else if(fFiberOptCouplingMode==G4MuonCounterConstruction::kGreaseCoupling)
			GrooveEnvelopMaterial= G4Material::GetMaterial("Silicone");	
		else if(fFiberOptCouplingMode==G4MuonCounterConstruction::kGlueCoupling)
			GrooveEnvelopMaterial= G4Material::GetMaterial("EpoxyTek");		

		//set the logical volume
		Groove_log= new G4LogicalVolume(Groove_box,GrooveEnvelopMaterial,"Groove_log",0,0,0);
		
		if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusGroove)
			Groove_log= new G4LogicalVolume(Groove_box,GrooveEnvelopMaterial,"Groove_log",0,0,0);
		else if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusFiber)
			Groove_log= new G4LogicalVolume(GrooveWithFiberHole,GrooveEnvelopMaterial,"Groove_log",0,0,0);
		
		//place groove envelop
		G4double offset_groove= fScintSizeZ/2-fGrooveSizeZ/2;
		if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusGroove||fStripDesignMode==G4MuonCounterConstruction::kScintPlusFiber) Groove_phys= new G4PVPlacement(0,G4ThreeVector(0,0,+offset_groove),Groove_log,"Groove_phys",ModuleHousing_log,false,0);


		//###################
		//##  SCINT COATING
		//###################
		//## up coating
		CoatUp_box = new G4Box("CoatUp_box",fScintSizeX/2,fScintSizeY/2,fScintCoatingThickness/2);
		CoatSide_box = new G4Box("CoatSide_box",fScintSizeX/2,fScintCoatingThickness/2,fScintCoatingSizeZ/2);
		CoatBottom_box = new G4Box("CoatBottom_box",fScintCoatingThickness/2,fScintCoatingSizeY/2,fScintCoatingSizeZ/2);


		//make end-cap with hole for fiber use
		G4SubtractionSolid* CoatBottomWithFiberHole = new G4SubtractionSolid("CoatBottomWithFiberHole", CoatBottom_box, FiberHole_tube, FiberRot, G4ThreeVector(0,0,offset_groove));	
		
		//place up & side coating
		//CoatUp_log= new G4LogicalVolume(CoatUp_box,G4Material::GetMaterial("Al"),"CoatUp_log",0,0,0);	
		CoatUp_log= new G4LogicalVolume(CoatUp_box,G4Material::GetMaterial("Polystyrene_TiO2mixed"),"CoatUp_log",0,0,0);	
		double offset_upcoat=	fScintSizeZ/2 + fScintCoatingThickness/2;
		CoatUp_phys= new G4PVPlacement(0,G4ThreeVector(0.,0.,+offset_upcoat),CoatUp_log,"CoatUp_phys",ModuleHousing_log,false,0);
		CoatUp_phys= new G4PVPlacement(0,G4ThreeVector(0.,0.,-offset_upcoat),CoatUp_log,"CoatUp_phys",ModuleHousing_log,false,1);
		
		//CoatSide_log= new G4LogicalVolume(CoatSide_box,G4Material::GetMaterial("Al"),"CoatSide_log",0,0,0);	
		CoatSide_log= new G4LogicalVolume(CoatSide_box,G4Material::GetMaterial("Polystyrene_TiO2mixed"),"CoatSide_log",0,0,0);	
		double offset_sidecoat=	fScintSizeY/2+fScintCoatingThickness/2;
		CoatSide_phys= new G4PVPlacement(0,G4ThreeVector(0.,+offset_sidecoat,0),CoatSide_log,"CoatSide_phys",ModuleHousing_log,false,0);
		CoatSide_phys= new G4PVPlacement(0,G4ThreeVector(0.,-offset_sidecoat,0),CoatSide_log,"CoatSide_phys",ModuleHousing_log,false,1);


		//place bottom coating
		G4double offset_bottomcoat=	fScintSizeX/2+fScintCoatingThickness/2;

		if(fStripReadoutMode==G4MuonCounterConstruction::kSingleReadout) {
			//full readout at one edge and reflective coating at other edge
			if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusFiber){
				//put an end-cap with fiber hole
				//CoatBottom_log= new G4LogicalVolume(CoatBottomWithFiberHole,G4Material::GetMaterial("Al"),"CoatBottom_log",0,0,0);	
				CoatBottom_log= new G4LogicalVolume(CoatBottomWithFiberHole,G4Material::GetMaterial("Polystyrene_TiO2mixed"),"CoatBottom_log",0,0,0);
				CoatBottom_phys= new G4PVPlacement(0,G4ThreeVector(-offset_bottomcoat,0,0),CoatBottom_log,"CoatBottom_phys",ModuleHousing_log,false,0);
			}
			else{
				//put an end-cap
				//CoatBottom_log= new G4LogicalVolume(CoatBottom_box,G4Material::GetMaterial("Al"),"CoatBottom_log",0,0,0);	
				CoatBottom_log= new G4LogicalVolume(CoatBottom_box,G4Material::GetMaterial("Polystyrene_TiO2mixed"),"CoatBottom_log",0,0,0);	
				CoatBottom_phys= new G4PVPlacement(0,G4ThreeVector(-offset_bottomcoat,0,0),CoatBottom_log,"CoatBottom_phys",ModuleHousing_log,false,0);
			}
		}//close if single readout
		else if(fStripReadoutMode==G4MuonCounterConstruction::kDoubleFiberReadout){
			//put end-cap with fiber hole at both edges
			//CoatBottom_log= new G4LogicalVolume(CoatBottomWithFiberHole,G4Material::GetMaterial("Al"),"CoatBottom_log",0,0,0);	
			CoatBottom_log= new G4LogicalVolume(CoatBottomWithFiberHole,G4Material::GetMaterial("Polystyrene_TiO2mixed"),"CoatBottom_log",0,0,0);	
			CoatBottom_phys= new G4PVPlacement(0,G4ThreeVector(+offset_bottomcoat,0,0),CoatBottom_log,"CoatBottom_phys",ModuleHousing_log,false,0);	
			CoatBottom_phys= new G4PVPlacement(0,G4ThreeVector(-offset_bottomcoat,0,0),CoatBottom_log,"CoatBottom_phys",ModuleHousing_log,false,1);	
		}
		else if(fStripReadoutMode==G4MuonCounterConstruction::kSingleFiberReadout){
			//put end-cap with fiber hole at both edges
			//CoatBottom_log= new G4LogicalVolume(CoatBottomWithFiberHole,G4Material::GetMaterial("Al"),"CoatBottom_log",0,0,0);	
			CoatBottom_log= new G4LogicalVolume(CoatBottomWithFiberHole,G4Material::GetMaterial("Polystyrene_TiO2mixed"),"CoatBottom_log",0,0,0);	
			CoatBottom_phys= new G4PVPlacement(0,G4ThreeVector(+offset_bottomcoat,0,0),CoatBottom_log,"CoatBottom_phys",ModuleHousing_log,false,0);	
			CoatBottom_phys= new G4PVPlacement(0,G4ThreeVector(-offset_bottomcoat,0,0),CoatBottom_log,"CoatBottom_phys",ModuleHousing_log,false,1);	
			//## add a reflective coating at the end of the fiber - TO BE DONE		
		}
			


		//###################
		//##  WLS FIBER
		//###################
		if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusFiber){
			//create the fiber inside the groove
			//core
			Fiber_tube= new G4Tubs("Fiber_tube",fFiberRmin,fFiberRmax,fFiberSizeZ/2,fFiberPhiMin,fFiberPhiMax);
  		Fiber_log = new G4LogicalVolume(Fiber_tube,G4Material::GetMaterial("Polystyrene"),"Fiber_log",0,0,0);
  		
  		//cladding 1 (first layer)
  		FiberClad1_tube= new G4Tubs("FiberClad1_tube",fClad1Rmin,fClad1Rmax,fClad1SizeZ/2,fClad1PhiMin,fClad1PhiMax);
  		FiberClad1_log = new G4LogicalVolume(FiberClad1_tube,G4Material::GetMaterial("PMMA"),"FiberClad1_log",0,0,0);
  		
			//cladding 2 (second layer)
  		FiberClad2_tube = new G4Tubs("FiberClad2_tube",fClad2Rmin,fClad2Rmax,fClad2SizeZ/2,fClad2PhiMin,fClad2PhiMax);  
  		FiberClad2_log = new G4LogicalVolume(FiberClad2_tube,G4Material::GetMaterial("fPethylene"),"FiberClad2_log",0,0,0);
  		
  		G4RotationMatrix* FiberRot = new G4RotationMatrix;
  		FiberRot->rotateY(M_PI/2.*CLHEP::rad);
			   
  		FiberClad2_phys= new G4PVPlacement(FiberRot,G4ThreeVector(0,0,+offset_groove),FiberClad2_log,"FiberClad2_phys",ModuleHousing_log,false,0); 
  		FiberClad1_phys= new G4PVPlacement(0,G4ThreeVector(),FiberClad1_log,"FiberClad1_phys",FiberClad2_log,false,0);
    	Fiber_phys= new G4PVPlacement(0,G4ThreeVector(),Fiber_log,"Fiber_phys",FiberClad1_log,false,0);
		}



		//###################
		//## PMT COUPLING
		//###################
		PMTOptCoupling_box= new G4Box("PMTOptCoupling_box",fPMTOptCouplingSizeX/2,fPMTOptCouplingSizeY/2,fPMTOptCouplingSizeZ/2);;	
		
		//select the coupling material	
		G4Material* CouplingMaterial=NULL;
		if(fStripOptCouplingMode==G4MuonCounterConstruction::kNoCoupling) 
			CouplingMaterial= G4Material::GetMaterial("Vacuum");
		else if(fStripOptCouplingMode==G4MuonCounterConstruction::kAirCoupling)
			CouplingMaterial= G4Material::GetMaterial("Air");
		else if(fStripOptCouplingMode==G4MuonCounterConstruction::kGreaseCoupling)
			CouplingMaterial= G4Material::GetMaterial("Silicone");
		
		//Insert an optical coupling layer?
		PMTOptCoupling_log= new G4LogicalVolume(PMTOptCoupling_box,CouplingMaterial,"PMTOptCoupling_log",0,0,0);;
		PMTOptCoupling_phys= NULL;
		G4double offset_pmtOptCoupling;
		if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusFiber)
			offset_pmtOptCoupling= max(fScintSizeX/2,fFiberSizeZ/2) + fPMTOptCouplingSizeX/2;
		else 
			offset_pmtOptCoupling= fScintSizeX/2 + fPMTOptCouplingSizeX/2;


		if(fStripOptCouplingMode!=G4MuonCounterConstruction::kNoCoupling){
			PMTOptCoupling_phys= new G4PVPlacement(0,G4ThreeVector(+offset_pmtOptCoupling,0.,0.),PMTOptCoupling_log,"PMTOptCoupling_phys",ModuleHousing_log,false,0);
			if(fStripReadoutMode==G4MuonCounterConstruction::kDoubleReadout||fStripReadoutMode==G4MuonCounterConstruction::kDoubleFiberReadout) 
				PMTOptCoupling_phys= new G4PVPlacement(0,G4ThreeVector(-offset_pmtOptCoupling,0.,0.),PMTOptCoupling_log,"PMTOptCoupling_phys",ModuleHousing_log,false,1);
		}//close if
			

  	//Create PMTs and photocathodes at the strip sides.
  	//Photocathode is a part of the PMT volume (placed inside PMT_log).
  	//The sensitive detector of the PMT is set on the photocathode.
    PMT_box = new G4Box("PMT_box",fPMTSizeX/2,fPMTSizeY/2,fPMTSizeZ/2);
    PMT_log= new G4LogicalVolume(PMT_box,G4Material::GetMaterial("Borosilicate"),"PMT_log",0,0,0);

		//create 2 PMTs
  	int copyNo=0;
		double offset_pmt=0.;

		if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusFiber){
			if(fStripOptCouplingMode==G4MuonCounterConstruction::kNoCoupling)
				offset_pmt= max(fScintSizeX/2,fFiberSizeZ/2) + fPMTSizeX/2;
			else 
				offset_pmt= max(fScintSizeX/2,fFiberSizeZ/2) + fPMTOptCouplingSizeX + fPMTSizeX/2;
		}
		else{
			if(fStripOptCouplingMode==G4MuonCounterConstruction::kNoCoupling)
				offset_pmt= fScintSizeX/2 + fPMTSizeX/2;
			else 
				offset_pmt= fScintSizeX/2 + fPMTOptCouplingSizeX + fPMTSizeX/2;	
		}

		Photocathode_box= new G4Box("Photocathode_box",fPhotocathodeSizeX/2,fPhotocathodeSizeY/2,fPhotocathodeSizeZ/2);   
		Photocathode_log= new G4LogicalVolume(Photocathode_box,G4Material::GetMaterial("Al"),"Photocathode_log",0,0,0);
	
  	//create Photocathodes
  	copyNo=0;
  	double offset_photocathode= fPMTSizeX/2-fPhotocathodeSizeX/2;
	
		Photocathode_phys= new G4PVPlacement(0,G4ThreeVector(+offset_photocathode,0.,0.),Photocathode_log,"Photocathode_phys",PMT_log,false,0);//inside PMT_log
		
		//PMT RIGHT
		G4RotationMatrix* zRotPMTRight = new G4RotationMatrix;	
		zRotPMTRight->rotateZ(M_PI*rad);
		PMT_phys= new G4PVPlacement(zRotPMTRight,G4ThreeVector(+offset_pmt,0,0),PMT_log,"PMT_phys",ModuleHousing_log,false,0);
    
		//PMT LEFT
		if(fStripReadoutMode==G4MuonCounterConstruction::kDoubleReadout||fStripReadoutMode==G4MuonCounterConstruction::kDoubleFiberReadout){		
      PMT_phys= new G4PVPlacement(0,G4ThreeVector(-offset_pmt,0,0),PMT_log,"PMT_phys",ModuleHousing_log,false,1);
		}
	

		//###############################
		//## SENSITIVE VOLUMES
		//###############################
		G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
		if(!Scint_SD){//determine if it has already been created
 			Scint_SD = new G4ScintillatorSD("/Detector/ScintSD");
			Scint_SD->SetScintillatorEnergyThreshold(fScintEnergyThreshold);
    	SDman->AddNewDetector(Scint_SD);    
 		}    
		Scint_log->SetSensitiveDetector(Scint_SD);

		/*
		//## Define MultiFunctionalDetector with name
		if(!MFDet){
  		MFDet = new G4MultiFunctionalDetector("StripModuleSD");
  		SDman->AddNewDetector(MFDet);                 // Register SD to SDManager
		}
		Scint_log->SetSensitiveDetector(MFDet);  // Assign SD to the logical volume.
		*/

		if(!PMT_SD){
 			PMT_SD = new G4PMTSD("/Detector/PMTSD");
    	SDman->AddNewDetector(PMT_SD);
 		}
		if(fStripReadoutMode==G4MuonCounterConstruction::kDoubleReadout)	
 			PMT_SD->InitPMTs(2); //let pmtSD know # of pmts
		else 
			PMT_SD->InitPMTs(1); //let pmtSD know # of pmts

		
 		PMT_SD->SetPMTPos(1,+offset_pmt,0,0);
		if(fStripReadoutMode==G4MuonCounterConstruction::kDoubleReadout) PMT_SD->SetPMTPos(2,-offset_pmt,0,0);

		//sensitive detector is not actually on the photocathode.
 		//processHits gets done manually by the stepping action.
 		//It is used to detect when photons hit and get absorbed&detected at the
 		//boundary to the photocathode (which doesnt get done by attaching it to a logical volume).
 		//It does however need to be attached to something or else it doesnt get
 		//reset at the begining of events   
 		Photocathode_log->SetSensitiveDetector(PMT_SD);

		/*
		//###############################
		//## PRIMITIVE SCORERS
		//############################### 	
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
		G4StripModulePSFlatSurfaceCurrent* muonSurfaceCurrentScorer = new G4StripModulePSFlatSurfaceCurrent("muonSurfaceCurrent", fCurrent_Out);
		muonSurfaceCurrentScorer->SetFilter(muonFilter);
		muonSurfaceCurrentScorer->DivideByArea(false);

		G4StripModulePSFlatSurfaceCurrent* emSurfaceCurrentScorer = new G4StripModulePSFlatSurfaceCurrent("emSurfaceCurrent", fCurrent_Out);
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

void G4StripModule::CopyValues(){

	//## Copy parameter from detector constructor
  fUpdated= fDetectorConstructor->GetUpdated();
	
	fStripDesignMode= fDetectorConstructor->GetStripDesignMode();
	fStripReadoutMode= fDetectorConstructor->GetStripReadoutMode();
	fStripOptCouplingMode= fDetectorConstructor->GetStripOptCouplingMode();
	fFiberOptCouplingMode= fDetectorConstructor->GetFiberOptCouplingMode();

  fScintSizeX= (fDetectorConstructor->GetScintSize()).x();
  fScintSizeY= (fDetectorConstructor->GetScintSize()).y();
  fScintSizeZ= (fDetectorConstructor->GetScintSize()).z();

	fScintCoatingThickness= fDetectorConstructor->GetHousingThickness();
  fScintCoatingSizeX= (fDetectorConstructor->GetScintCoatingSize()).x();
	fScintCoatingSizeY= (fDetectorConstructor->GetScintCoatingSize()).y();
	fScintCoatingSizeZ= (fDetectorConstructor->GetScintCoatingSize()).z();
	fScintCoatingReflectivity= fDetectorConstructor->GetHousingReflectivity();

	fScintCoatingCapSizeX= fScintCoatingThickness;
	fScintCoatingCapSizeY= fScintCoatingSizeY;
	fScintCoatingCapSizeZ= fScintCoatingSizeZ;

	fScintEnergyThreshold= fDetectorConstructor->GetStripEnergyThreshold();

	fGrooveSizeX= fDetectorConstructor->GetGrooveSize().x();
  fGrooveSizeY= fDetectorConstructor->GetGrooveSize().y();
  fGrooveSizeZ= fDetectorConstructor->GetGrooveSize().z();

  fFiberRmin= 0.*CLHEP::mm;
	fFiberRmax= fDetectorConstructor->GetFiberRadius();
	fFiberSizeZ= fDetectorConstructor->GetFiberSizeZ();
	fFiberPhiMin= 0.*CLHEP::deg;
	fFiberPhiMax= 360.*CLHEP::deg;
	
	fClad1Rmin= 0.*CLHEP::mm;
	fClad1Rmax= fDetectorConstructor->GetClad1Radius();
	fClad1SizeZ= fFiberSizeZ;
	fClad1PhiMin= 0.*CLHEP::deg;
	fClad1PhiMax= 360.*CLHEP::deg;
		
	fClad2Rmin= 0.*CLHEP::mm;
	fClad2Rmax= fDetectorConstructor->GetClad2Radius();
	fClad2SizeZ= fFiberSizeZ;
	fClad2PhiMin= 0.*CLHEP::deg;
	fClad2PhiMax= 360.*CLHEP::deg;


  fPMTSizeX= (fDetectorConstructor->GetPMTSize()).x();
  fPMTSizeY= (fDetectorConstructor->GetPMTSize()).y();
  fPMTSizeZ= (fDetectorConstructor->GetPMTSize()).z();

  fPhotocathodeSizeX= (fDetectorConstructor->GetPhotocathodeSize()).x();
  fPhotocathodeSizeY= (fDetectorConstructor->GetPhotocathodeSize()).y();
  fPhotocathodeSizeZ= (fDetectorConstructor->GetPhotocathodeSize()).z();
  
	fPMTOptCouplingSizeX= (fDetectorConstructor->GetPMTOpticalCouplingSize()).x();
  fPMTOptCouplingSizeY= (fDetectorConstructor->GetPMTOpticalCouplingSize()).y();
  fPMTOptCouplingSizeZ= (fDetectorConstructor->GetPMTOpticalCouplingSize()).z();

  fUseRealReflectivity= fDetectorConstructor->IsRealReflectivity();
	fUseRealPhotocathode= fDetectorConstructor->IsRealPhotocathode();
	fUseFiber= fDetectorConstructor->IsFiber();

	fUseDoubleReadout= fDetectorConstructor->IsDoubleReadout();
	fCoatingSurfaceType= fDetectorConstructor->GetHousingSurfaceType();
	
	fSurfaceSeparationGap= fDetectorConstructor->GetSurfaceSeparationGap();

}//close function


void G4StripModule::VisAttributes(){
	
  G4VisAttributes* ModuleHousing_va = new G4VisAttributes(G4Colour(0.96,0.11,0.71));//darkpurple
  ModuleHousing_log->SetVisAttributes(G4VisAttributes::Invisible);
	//ModuleHousing_log->SetVisAttributes(ModuleHousing_va);	
	
	G4VisAttributes* Scint_va = new G4VisAttributes(G4Colour(0.75,0.75,0.75,0.3));//silver
	Scint_log->SetVisAttributes(Scint_va); 	


	G4VisAttributes* Coat_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8,0.3));//gray
	CoatUp_log->SetVisAttributes(Coat_va);
	CoatSide_log->SetVisAttributes(Coat_va);
	if(fStripReadoutMode==G4MuonCounterConstruction::kSingleReadout||fStripReadoutMode==G4MuonCounterConstruction::kSingleFiberReadout||fStripReadoutMode==G4MuonCounterConstruction::kDoubleFiberReadout)
		CoatBottom_log->SetVisAttributes(Coat_va);

 
  G4VisAttributes* Groove_va= NULL;
	G4VisAttributes* Fiber_va= NULL;	
	G4VisAttributes* FiberClad1_va= NULL;
	G4VisAttributes* FiberClad2_va= NULL;
	

	if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusFiber ||
     fStripDesignMode==G4MuonCounterConstruction::kScintPlusGroove
	){
		Groove_va = new G4VisAttributes(G4Colour(1,1,1));//white
		Groove_log->SetVisAttributes(G4VisAttributes::Invisible);
		//Groove_log->SetVisAttributes(Groove_va);
	}


	if(fStripDesignMode==G4MuonCounterConstruction::kScintPlusFiber){
		Fiber_va = new G4VisAttributes(G4Colour(0,1,1,0.3));//cyan
		Fiber_log->SetVisAttributes(Fiber_va);
		
		FiberClad1_va = new G4VisAttributes(G4Colour(0,0.8,1,0.3));//dark cyan 1 
		FiberClad1_log->SetVisAttributes(FiberClad1_va);

		FiberClad2_va = new G4VisAttributes(G4Colour(0,0.7,1,0.3));//dark cyan 2
		FiberClad2_log->SetVisAttributes(FiberClad2_va);
	}


	G4VisAttributes* PMT_va = new G4VisAttributes(G4Colour(0,0,0));//black //green 010
	PMT_log->SetVisAttributes(PMT_va);

	G4VisAttributes* Photocathode_va = new G4VisAttributes(G4Colour(0,0,1));//blue//yellow 110
	Photocathode_log->SetVisAttributes(Photocathode_va);

	G4VisAttributes* PMTOptCoupling_va= NULL;
	if(fStripOptCouplingMode!=G4MuonCounterConstruction::kNoCoupling){
		PMTOptCoupling_va = new G4VisAttributes(G4Colour(1,1,0));
		PMTOptCoupling_log->SetVisAttributes(PMTOptCoupling_va);
	}

  
}//close function

void G4StripModule::SurfaceProperties(){ 

	//#############################
	//##  Scintillator Coating  ###
	//#############################
	// ==> reflectivity=fixed or real spectrum, efficiency=0 (just reflect and absorb, not detect)
	//IDEAL
  const G4int NUMENTRIES_ScintCoat_ideal = 2;
  G4double PhotonEnergy_ScintCoat_ideal[NUMENTRIES_ScintCoat_ideal] = {2.0*CLHEP::eV, 3.5*CLHEP::eV};  
  G4double REFL_ScintCoat_ideal[NUMENTRIES_ScintCoat_ideal] = {fScintCoatingReflectivity, fScintCoatingReflectivity};
  G4double EFF_ScintCoat_ideal[NUMENTRIES_ScintCoat_ideal] = {0.0, 0.0}; 
	G4MaterialPropertiesTable* IdealScintCoatingPT = new G4MaterialPropertiesTable(); 
	IdealScintCoatingPT->AddProperty("REFLECTIVITY", PhotonEnergy_ScintCoat_ideal, REFL_ScintCoat_ideal, NUMENTRIES_ScintCoat_ideal);
 	IdealScintCoatingPT->AddProperty("EFFICIENCY", PhotonEnergy_ScintCoat_ideal, EFF_ScintCoat_ideal, NUMENTRIES_ScintCoat_ideal);

	G4OpticalSurface* OpScintCoatingSurface_ideal;

	if(fCoatingSurfaceType==1) {
		OpScintCoatingSurface_ideal= new G4OpticalSurface("ScintCoatingSurface_ideal",unified,polished,dielectric_metal);//polished
	}
	else if(fCoatingSurfaceType==2){
		OpScintCoatingSurface_ideal= new G4OpticalSurface("ScintCoatingSurface_ideal",unified,ground,dielectric_metal);//polished
	}
	else{
		cerr<<"Invalid coating surface type"<<endl;
		exit(1);
	}
	OpScintCoatingSurface_ideal->SetMaterialPropertiesTable(IdealScintCoatingPT);

	//REAL
// ==> reflectivity= realistic spectrum, efficiency=0 (just reflect and absorb, not detect)	
	const G4int NUMENTRIES_ScintCoat_real = 28;//15
  G4double PhotonEnergy_ScintCoat_real[NUMENTRIES_ScintCoat_real] ={
1.97742*CLHEP::eV,2.02093*CLHEP::eV,2.04831*CLHEP::eV,2.08552*CLHEP::eV,2.11433*CLHEP::eV,2.13398*CLHEP::eV,2.16415*CLHEP::eV,2.19558*CLHEP::eV,2.21677*CLHEP::eV,2.23838*CLHEP::eV,2.2716*CLHEP::eV,2.31746*CLHEP::eV,2.35309*CLHEP::eV,2.41543*CLHEP::eV,2.4282*CLHEP::eV,2.52154*CLHEP::eV,2.59327*CLHEP::eV,2.66862*CLHEP::eV,2.73273*CLHEP::eV,2.78303*CLHEP::eV,2.87067*CLHEP::eV,2.96329*CLHEP::eV,3.04256*CLHEP::eV,3.14681*CLHEP::eV,3.25931*CLHEP::eV,3.35546*CLHEP::eV,3.43066*CLHEP::eV,3.53735*CLHEP::eV};
  
  G4double REFL_ScintCoat_real[NUMENTRIES_ScintCoat_real] ={
0.8841,0.886,0.8841,0.8841,0.8878,0.8841,0.886,0.8822,0.8859,0.8878,0.8897,0.8915,0.8896,0.8859,0.8878,0.8877,0.8858,0.884,0.8821,0.8783,0.8708,0.8689,0.8652,0.8652,0.8652,0.8595,0.8595,0.8558};
  
	G4double EFF_ScintCoat_real[NUMENTRIES_ScintCoat_real] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  G4MaterialPropertiesTable* RealScintCoatingPT = new G4MaterialPropertiesTable(); 
	RealScintCoatingPT->AddProperty("REFLECTIVITY", PhotonEnergy_ScintCoat_real, REFL_ScintCoat_real, NUMENTRIES_ScintCoat_real);
 	RealScintCoatingPT->AddProperty("EFFICIENCY", PhotonEnergy_ScintCoat_real, EFF_ScintCoat_real, NUMENTRIES_ScintCoat_real);
	
	G4OpticalSurface* OpScintCoatingSurface_real;

	if(fCoatingSurfaceType==1) {
		OpScintCoatingSurface_real= new G4OpticalSurface("RealScintCoatingSurface",unified,polished,dielectric_metal);//polished
	}
	else if(fCoatingSurfaceType==2){
		OpScintCoatingSurface_real= new G4OpticalSurface("RealScintCoatingSurface",unified,ground,dielectric_metal);//polished
	}
	else{
		cerr<<"Invalid coating surface type"<<endl;
		exit(1);
	}
	OpScintCoatingSurface_real->SetMaterialPropertiesTable(RealScintCoatingPT);  

  
  

  //########################
	//##   PMT             ###
	//########################
  const G4int NUMENTRIES_PMT= 2;
  G4double PhotonEnergy_PMT[NUMENTRIES_PMT]= {2.0*CLHEP::eV,3.5*CLHEP::eV}; 
	G4double RIND_PMT[NUMENTRIES_PMT]={ 1.474, 1.474};
  G4MaterialPropertiesTable* PMTPT = new G4MaterialPropertiesTable();
  PMTPT->AddProperty("RINDEX",PhotonEnergy_PMT,RIND_PMT,NUMENTRIES_PMT);;


	//########################
	//##   Photocathode    ###
	//########################
	//Enables 'detection' of photons
  //ideal PMT
  const G4int NUMENTRIES_Photocathode_ideal = 2;
  G4double PhotonEnergy_Photocathode_ideal[NUMENTRIES_Photocathode_ideal]= {2.0*CLHEP::eV,3.5*CLHEP::eV}; 
  G4double Photocathode_EFF_ideal[NUMENTRIES_Photocathode_ideal]={1.0,1.0};//change to a realistic QE
  G4double Photocathode_REFL_ideal[NUMENTRIES_Photocathode_ideal]={0.0,0.0};
  G4MaterialPropertiesTable* IdealPhotocathodePT = new G4MaterialPropertiesTable();
  IdealPhotocathodePT->AddProperty("EFFICIENCY",PhotonEnergy_Photocathode_ideal,Photocathode_EFF_ideal,NUMENTRIES_Photocathode_ideal);
  IdealPhotocathodePT->AddProperty("REFLECTIVITY",PhotonEnergy_Photocathode_ideal,Photocathode_REFL_ideal,NUMENTRIES_Photocathode_ideal);

  G4OpticalSurface* OpSurfacePhotocathode_ideal= new G4OpticalSurface("PhotocathodeOpSurf_ideal",glisur,polished,dielectric_metal);
  OpSurfacePhotocathode_ideal->SetMaterialPropertiesTable(IdealPhotocathodePT);


  //Realistic QE from Thorn EMI (now Phototube) model 9826B (QB)
  const G4int NUMENTRIES_Photocathode_real = 80;
	G4double PhotonEnergy_Photocathode_real[NUMENTRIES_Photocathode_real]={1.92134*CLHEP::eV,1.93604*CLHEP::eV,1.95036*CLHEP::eV,1.99556*CLHEP::eV,2.04224*CLHEP::eV,2.09965*CLHEP::eV,2.14172*CLHEP::eV,2.17631*CLHEP::eV,2.20691*CLHEP::eV,2.22035*CLHEP::eV,2.23395*CLHEP::eV,2.24243*CLHEP::eV,2.27077*CLHEP::eV,2.27954*CLHEP::eV,2.29388*CLHEP::eV,2.32267*CLHEP::eV,2.34242*CLHEP::eV,2.3625*CLHEP::eV,2.38294*CLHEP::eV,2.39907*CLHEP::eV,2.40372*CLHEP::eV,2.41967*CLHEP::eV,2.44689*CLHEP::eV,2.46391*CLHEP::eV,2.48067*CLHEP::eV,2.49767*CLHEP::eV,2.50271*CLHEP::eV,2.51438*CLHEP::eV,2.53132*CLHEP::eV,2.54902*CLHEP::eV,2.56643*CLHEP::eV,2.58516*CLHEP::eV,2.61625*CLHEP::eV,2.64133*CLHEP::eV,2.66117*CLHEP::eV,2.68771*CLHEP::eV,2.71478*CLHEP::eV,2.73575*CLHEP::eV,2.75643*CLHEP::eV,2.78553*CLHEP::eV,2.83069*CLHEP::eV,2.84563*CLHEP::eV,2.88469*CLHEP::eV,2.91659*CLHEP::eV,2.95834*CLHEP::eV,2.99189*CLHEP::eV,3.01811*CLHEP::eV,3.0538*CLHEP::eV,3.0996*CLHEP::eV,3.14681*CLHEP::eV,3.20538*CLHEP::eV,3.2774*CLHEP::eV,3.34189*CLHEP::eV,3.37555*CLHEP::eV,3.40897*CLHEP::eV,3.46809*CLHEP::eV,3.56789*CLHEP::eV,3.66059*CLHEP::eV,3.70323*CLHEP::eV,3.73334*CLHEP::eV,3.7628*CLHEP::eV,3.77885*CLHEP::eV,3.79737*CLHEP::eV,3.8314*CLHEP::eV,3.84924*CLHEP::eV,3.86725*CLHEP::eV,3.88787*CLHEP::eV,3.90624*CLHEP::eV,3.92603*CLHEP::eV,3.94603*CLHEP::eV,3.96495*CLHEP::eV,3.96876*CLHEP::eV,3.9713*CLHEP::eV,3.99048*CLHEP::eV,4.00854*CLHEP::eV,4.04912*CLHEP::eV,4.07039*CLHEP::eV,4.07708*CLHEP::eV,4.08245*CLHEP::eV,4.08514*CLHEP::eV};

  G4double Photocathode_EFF_real[NUMENTRIES_Photocathode_real]={0.002064,0.002565,0.004543,0.009493,0.01592,0.02481,0.03419,0.04357,0.05343,0.05836,0.06231,0.06675,0.07464,0.07957,0.0845,0.09387,0.1003,0.1072,0.1121,0.1151,0.1185,0.1235,0.1294,0.1328,0.1368,0.1407,0.1447,0.1466,0.1521,0.157,0.1619,0.1654,0.1723,0.1772,0.1812,0.1851,0.1901,0.1925,0.196,0.1989,0.2029,0.2054,0.2078,0.2108,0.2133,0.2153,0.2167,0.2177,0.2202,0.2217,0.2227,0.2242,0.2242,0.2242,0.2237,0.2222,0.2193,0.2139,0.21,0.207,0.2036,0.2001,0.1937,0.1868,0.1824,0.177,0.1701,0.1652,0.1583,0.1533,0.1479,0.1435,0.1396,0.1361,0.1317,0.1208,0.1149,0.1066,0.1011,0.03662};

  G4double Photocathode_REFL_real[NUMENTRIES_Photocathode_real]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


  G4MaterialPropertiesTable* RealPhotocathodePT = new G4MaterialPropertiesTable();
  RealPhotocathodePT->AddProperty("EFFICIENCY",PhotonEnergy_Photocathode_real,Photocathode_EFF_real,NUMENTRIES_Photocathode_real);
  RealPhotocathodePT->AddProperty("REFLECTIVITY",PhotonEnergy_Photocathode_real,Photocathode_REFL_real,NUMENTRIES_Photocathode_real);

  G4OpticalSurface* OpSurfacePhotocathode_real= new G4OpticalSurface("PhotocathodeOpSurf_real",glisur,polished,dielectric_metal);
  OpSurfacePhotocathode_real->SetMaterialPropertiesTable(RealPhotocathodePT);

	//#############################
	//##  SURFACES             ####
	//#############################
	//### coating interface
	G4OpticalSurface* CoatSurfaceType=NULL;
	if(fUseRealReflectivity) {
		CoatSurfaceType= OpScintCoatingSurface_real;		
	}
	else {
		CoatSurfaceType= OpScintCoatingSurface_ideal;
	}

	new G4LogicalSkinSurface("CoatUpSurface",CoatUp_log,CoatSurfaceType);
	new G4LogicalSkinSurface("CoatSideSurface",CoatSide_log,CoatSurfaceType);
	if(fStripReadoutMode==G4MuonCounterConstruction::kSingleReadout||fStripReadoutMode==G4MuonCounterConstruction::kSingleFiberReadout||fStripReadoutMode==G4MuonCounterConstruction::kDoubleFiberReadout) 
		new G4LogicalSkinSurface("CoatBottomSurface",CoatBottom_log,CoatSurfaceType);

	//### PMT photocathode interface
	G4OpticalSurface* PhotocathodeSurfaceType=NULL;
	if(fUseRealPhotocathode)
		PhotocathodeSurfaceType= OpSurfacePhotocathode_real; 		
	else
		PhotocathodeSurfaceType= OpSurfacePhotocathode_ideal; 

	new G4LogicalSkinSurface("PhotocathodeSurface",Photocathode_log,PhotocathodeSurfaceType); 
 
}//close function



