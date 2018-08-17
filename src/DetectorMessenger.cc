/**
* @file DetectorMessenger.cc
* @class DetectorMessenger
* @brief Create the UI commands to handle the detector geometry parameters, i.e. from a configuration file  
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "DetectorMessenger.hh"
#include "G4MuonCounterConstruction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4Scintillation.hh"

#include "vector"

using namespace G4MuonCounterSimulatorUSC;

DetectorMessenger::DetectorMessenger(G4MuonCounterConstruction* det)
:Detector(det)
{
  //Setup a command directory for detector controls with guidance
  detectorDir = new G4UIdirectory("/StripSimulation/det/");
  detectorDir->SetGuidance("Detector geometry control");

  volumesDir = new G4UIdirectory("/StripSimulation/det/volumes/");
  volumesDir->SetGuidance("Enable/disable volumes");

	updateCmd = new G4UIcommand("/StripSimulation/det/update",this);
  updateCmd->SetGuidance("Update the detector geometry with changed values.");
  updateCmd->SetGuidance("Must be run before beamOn if detector has been changed.");
   
  //Various commands for modifying detector geometry
	stripDesignModeCmd = new G4UIcmdWithAnInteger("/StripSimulation/det/StripDesignMode",this);
  stripDesignModeCmd->SetGuidance("Set the type of strip design: 1-scint, 2-scint+groove, 3-scint+groove+fiber");
  stripDesignModeCmd->SetParameterName("StripDesignMode",true,true);
	stripDesignModeCmd->SetDefaultValue(3);
	stripDesignModeCmd->SetRange("StripDesignMode==1||StripDesignMode==2||StripDesignMode==3");

	stripReadoutModeCmd = new G4UIcmdWithAnInteger("/StripSimulation/det/StripReadoutMode",this);
  stripReadoutModeCmd->SetGuidance("Set the type of strip readout: 1-singleread, 2-doubleread");
  stripReadoutModeCmd->SetParameterName("StripReadoutMode",true,true);
	stripReadoutModeCmd->SetDefaultValue(2);
	stripReadoutModeCmd->SetRange("StripReadoutMode==1||StripReadoutMode==2||StripReadoutMode==3||StripReadoutMode==4");	

	stripOptCouplingModeCmd = new G4UIcmdWithAnInteger("/StripSimulation/det/StripOptCouplingMode",this);
  stripOptCouplingModeCmd->SetGuidance("Set the type of strip readout: 0-nocoupling, 1-aircoupling, 2-greasecoupling");
  stripOptCouplingModeCmd->SetParameterName("StripOptCouplingMode",true,true);
	stripOptCouplingModeCmd->SetDefaultValue(0);
	stripOptCouplingModeCmd->SetRange("StripOptCouplingMode==0||StripOptCouplingMode==1||StripOptCouplingMode==2");	

  stripSizeCmd = new G4UIcmdWith3VectorAndUnit("/StripSimulation/det/StripSize",this);
  stripSizeCmd->SetGuidance("Set the dimensions of the scintillator strip.");
  stripSizeCmd->SetParameterName("strip_x","strip_y","strip_z",true,true);
  stripSizeCmd->SetDefaultUnit("cm");
	stripSizeCmd->SetDefaultValue(G4ThreeVector(100.0*CLHEP::cm,1.0*CLHEP::cm,1.0*CLHEP::cm));
	stripSizeCmd->SetRange("strip_x>0 && strip_y>0 && strip_z>0");

	grooveSizeCmd = new G4UIcmdWith3VectorAndUnit("/StripSimulation/det/GrooveSize",this);
  grooveSizeCmd->SetGuidance("Set the dimensions of the groove in the scintillator.");
  grooveSizeCmd->SetParameterName("groove_x","groove_y","groove_z",true,true);
  grooveSizeCmd->SetDefaultUnit("cm");
	grooveSizeCmd->SetDefaultValue(G4ThreeVector(100.0*CLHEP::cm,0.12*CLHEP::cm,0.12*CLHEP::cm));
	grooveSizeCmd->SetRange("groove_x>0 && groove_y>0 && groove_z>0");

	fiberLengthCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/fiberLength",this);
  fiberLengthCmd->SetGuidance("Set the length of the fiber.");
  fiberLengthCmd->SetParameterName("FiberLength",true,true);
  fiberLengthCmd->SetDefaultUnit("cm");
	fiberLengthCmd->SetDefaultValue(100.0*cm);
	fiberLengthCmd->SetRange("FiberLength>0");

	fiberRadiusCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/fiberRadius",this);
  fiberRadiusCmd->SetGuidance("Set the radius of the fiber.");
  fiberRadiusCmd->SetParameterName("FiberRadius",true,true);
  fiberRadiusCmd->SetDefaultUnit("mm");
	fiberRadiusCmd->SetDefaultValue(0.5*CLHEP::mm);
	fiberRadiusCmd->SetRange("FiberRadius>0");

	fiberAbsLengthCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/fiberAbsLength",this);
  fiberAbsLengthCmd->SetGuidance("Set the absorption length of the fiber in the wavelength region of absorption.");
  fiberAbsLengthCmd->SetParameterName("FiberAbsLength",true,true);
  fiberAbsLengthCmd->SetDefaultUnit("mm");
	fiberAbsLengthCmd->SetDefaultValue(0.1*CLHEP::mm);
	fiberAbsLengthCmd->SetRange("FiberAbsLength>0");

  housingThicknessCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/housingThickness",this);
  housingThicknessCmd->SetGuidance("Set the thickness of the housing.");
  housingThicknessCmd->SetParameterName("HousingThickness",true,true);
  housingThicknessCmd->SetDefaultUnit("mm");
	housingThicknessCmd->SetDefaultValue(0.12);

  housingReflectivityCmd = new G4UIcmdWithADouble("/StripSimulation/det/housingReflectivity",this);
  housingReflectivityCmd->SetGuidance("Set the reflectivity of the housing.");
	housingReflectivityCmd->SetParameterName("reflectivity",true,true);
	housingReflectivityCmd->SetDefaultValue(0.9);
  housingReflectivityCmd->SetRange("reflectivity>=0.0 && reflectivity<=1.0");

  housingUseRealReflectivityCmd = new G4UIcmdWithABool("/StripSimulation/det/housingUseRealReflectivity",this);
  housingUseRealReflectivityCmd->SetGuidance("Enable/Disable the use of real coating reflectivity");
  housingUseRealReflectivityCmd->SetParameterName("UseRealReflectivity",true,true);
	housingUseRealReflectivityCmd->SetDefaultValue(true);

	housingSetSurfaceTypeCmd = new G4UIcmdWithAnInteger("/StripSimulation/det/housingSetSurfaceType",this);
  housingSetSurfaceTypeCmd->SetGuidance("Set the type of coating surface: 1-polished, 2-ground, ...");
  housingSetSurfaceTypeCmd->SetParameterName("SetSurfaceType",true,true);
	housingSetSurfaceTypeCmd->SetDefaultValue(1);
  housingSetSurfaceTypeCmd->SetRange("SetSurfaceType==1||SetSurfaceType==2");	


  pmtSizeCmd = new G4UIcmdWith3VectorAndUnit("/StripSimulation/det/pmtSize",this);
  pmtSizeCmd->SetGuidance("Set the size of the PMTs.");
  pmtSizeCmd->SetParameterName("pmt_x","pmt_y","pmt_z",true,true);
  pmtSizeCmd->SetDefaultUnit("mm");
  pmtSizeCmd->SetDefaultValue(G4ThreeVector(10.0*CLHEP::cm,1.0*CLHEP::cm,1.0*CLHEP::cm));
	pmtSizeCmd->SetRange("pmt_x>0 && pmt_y>0 && pmt_z>0");

	photocathodeSizeCmd = new G4UIcmdWith3VectorAndUnit("/StripSimulation/det/photocathodeSize",this);
  photocathodeSizeCmd->SetGuidance("Set the size of the PMT photocathode.");
  photocathodeSizeCmd->SetParameterName("photocathode_x","photocathode_y","photocathode_z",true,true);	
  photocathodeSizeCmd->SetDefaultUnit("mm");
  photocathodeSizeCmd->SetDefaultValue(G4ThreeVector(1.0*CLHEP::mm,1.0*CLHEP::cm,1.0*CLHEP::cm));
	
 	pmtOptCouplingSizeCmd = new G4UIcmdWith3VectorAndUnit("/StripSimulation/det/pmtOptCouplingSize",this);
  pmtOptCouplingSizeCmd->SetGuidance("Set the size of the PMT optical coupling region");
  pmtOptCouplingSizeCmd->SetParameterName("pmtOptCoupling_x","pmtOptCoupling_y","pmtOptCoupling_z",true,true);
  pmtOptCouplingSizeCmd->SetDefaultUnit("mm");
  pmtOptCouplingSizeCmd->SetDefaultValue(G4ThreeVector(0.5*CLHEP::mm,1.0*CLHEP::cm,1.0*CLHEP::cm));
	pmtOptCouplingSizeCmd->SetRange("pmtOptCoupling_x>0 && pmtOptCoupling_y>0 && pmtOptCoupling_z>0");


  wlsCmd = new G4UIcmdWithABool("/StripSimulation/det/wlsFiber",this);
  wlsCmd->SetGuidance("Enable/Disable the WLS fiber");
  wlsCmd->SetParameterName("UseWLSFiber",true,true);
	wlsCmd->SetDefaultValue(true);

	readoutModeCmd = new G4UIcmdWithABool("/StripSimulation/det/doubleReadout",this);
  readoutModeCmd->SetGuidance("Enable/Disable the double readout of the strip");
  readoutModeCmd->SetParameterName("UseDoubleReadout",true,true);
	readoutModeCmd->SetDefaultValue(true);

	
	pmtUseRealQECmd = new G4UIcmdWithABool("/StripSimulation/det/pmtUseRealQE",this);
  pmtUseRealQECmd->SetGuidance("Enable/Disable the use of realistic pmt QE");
  pmtUseRealQECmd->SetParameterName("pmtUseRealQE",true,true);
	pmtUseRealQECmd->SetDefaultValue(false);

	yPlaneInSuperModuleCmd = new G4UIcmdWithABool("/StripSimulation/det/enableYPlaneInSuperModule",this);
  yPlaneInSuperModuleCmd->SetGuidance("Enable/Disable the y plane in the supermodule");
  yPlaneInSuperModuleCmd->SetParameterName("UseYPlaneInSuperModule",true,true);
	yPlaneInSuperModuleCmd->SetDefaultValue(true);

  stripNumberCmd = new G4UIcmdWithAnInteger("/StripSimulation/det/stripNumber",this);
  stripNumberCmd->SetGuidance("Set the number of strips per plane");
  stripNumberCmd->SetParameterName("SetStripNumber",true,true);
	stripNumberCmd->SetDefaultValue(1);

	planeNumberCmd = new G4UIcmdWithAnInteger("/StripSimulation/det/planeNumber",this);
  planeNumberCmd->SetGuidance("Set the number of detector planes");
  planeNumberCmd->SetParameterName("SetPlaneNumber",true,true);
	planeNumberCmd->SetDefaultValue(1);

	planeDistanceCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/planeDistance",this);
  planeDistanceCmd->SetGuidance("Set the distance among detector planes.");
  planeDistanceCmd->SetParameterName("planeDistance",true,true);
  planeDistanceCmd->SetDefaultUnit("cm");
	planeDistanceCmd->SetDefaultValue(10.0);

	planeXYDistanceCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/planeXYDistance",this);
  planeXYDistanceCmd->SetGuidance("Set the distance among XY detector planes.");
  planeXYDistanceCmd->SetParameterName("planeXYDistance",true,true);
  planeXYDistanceCmd->SetDefaultUnit("cm");
	planeXYDistanceCmd->SetDefaultValue(0.0);

	planeTiltAngleCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/planeTiltAngle",this);
  planeTiltAngleCmd->SetGuidance("Set the tilt angle of detector planes.");
  planeTiltAngleCmd->SetParameterName("planeTiltAngle",true,true);
  planeTiltAngleCmd->SetDefaultUnit("deg");
	planeTiltAngleCmd->SetDefaultValue(0.0);

	superplaneThicknessCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/superplaneThickness",this);
  superplaneThicknessCmd->SetGuidance("Set the thickness of the superplane casing.");
  superplaneThicknessCmd->SetParameterName("superplaneThickness",true,true);
  superplaneThicknessCmd->SetDefaultUnit("mm");
	superplaneThicknessCmd->SetDefaultValue(3.0);

	superplaneSizeZCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/det/superplaneSizeZ",this);
  superplaneSizeZCmd->SetGuidance("Set the Z size of the superplane casing.");
  superplaneSizeZCmd->SetParameterName("superplaneSizeZ",true,true);
  superplaneSizeZCmd->SetDefaultUnit("cm");
	superplaneSizeZCmd->SetDefaultValue(6.0);

	worldMaterialCmd = new G4UIcmdWithAString("/StripSimulation/det/worldMaterial",this);
  worldMaterialCmd->SetGuidance("Set the material of the world.");
  worldMaterialCmd->SetParameterName("worldMaterial",true);
  worldMaterialCmd->SetDefaultValue("Air");
	G4String candidateList; 
 	candidateList += "Air ";
	candidateList += "Vacuum ";
	candidateList += "MalargueSoil ";	
  worldMaterialCmd->SetCandidates(candidateList);
	
  //worldMaterialCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	//worldMaterialCmd = new G4UIcommand("/StripSimulation/det/worldMaterial",this);
  //worldMaterialCmd->SetGuidance("Set the material of the world.");

	superplaneCasingMaterialCmd = new G4UIcmdWithAString("/StripSimulation/det/superplaneCasingMaterial",this);
  superplaneCasingMaterialCmd->SetGuidance("Set the material of the superplane casing.");
  superplaneCasingMaterialCmd->SetParameterName("superplaneCasingMaterial",true);
  superplaneCasingMaterialCmd->SetDefaultValue("PVC");

	G4String candidateCasingList; 
 	candidateCasingList += "Air ";
	candidateCasingList += "Vacuum ";
	candidateCasingList += "PVC ";
	candidateCasingList += "GlassFiber ";		
  superplaneCasingMaterialCmd->SetCandidates(candidateCasingList);

/*
  defaultsCmd = new G4UIcommand("/StripSimulation/det/defaults",this);
  defaultsCmd->SetGuidance("Set all detector geometry values to defaults.");
  defaultsCmd->SetGuidance("(Update still required)");
*/

  scintillationYieldCmd=new G4UIcmdWithADouble("/StripSimulation/det/ScintillationYield",this);
  scintillationYieldCmd->SetGuidance("Set scintillation yield of scintillator.");
  scintillationYieldCmd->SetGuidance("Specified in photons/MeV");
	scintillationYieldCmd->SetParameterName("ScintYield",true,true);
	scintillationYieldCmd->SetDefaultValue(11136.0);
  scintillationYieldCmd->SetRange("ScintYield>0.0");

	scintillationWeightCmd=new G4UIcmdWithADouble("/StripSimulation/det/ScintillationWeight",this);
  scintillationWeightCmd->SetGuidance("Set scintillation weight.");
  scintillationWeightCmd->SetGuidance("Used to rescale simulation results, when the yield is not the true yield");
	scintillationWeightCmd->SetParameterName("ScintWeight",true,true);
	scintillationWeightCmd->SetDefaultValue(1.0);
  scintillationWeightCmd->SetRange("ScintWeight>0.0");

/*
  WLSScintYield = new G4UIcmdWithADouble("/StripSimulation/det/WLSScintYield",this);
  WLSScintYield->SetGuidance("Set scintillation yield of WLS Slab");
  WLSScintYield->SetGuidance("Specified in photons/MeV");
  WLSScintYield->SetParameterName("WLSYield",true,true);
  WLSScintYield->SetDefaultUnit("1/MeV");
	WLSScintYield->SetDefaultValue(1.0);
  WLSScintYield->SetRange("WLSYield>0.0");
*/
}


DetectorMessenger::~DetectorMessenger()
{
	delete updateCmd;
	delete stripDesignModeCmd;
	delete stripReadoutModeCmd;
	delete stripOptCouplingModeCmd;
  delete stripSizeCmd;
	delete grooveSizeCmd;
	delete fiberLengthCmd;
	delete fiberRadiusCmd;
	delete fiberAbsLengthCmd;
  delete housingThicknessCmd;
  delete housingReflectivityCmd;
  delete housingUseRealReflectivityCmd;
  delete housingSetSurfaceTypeCmd;
  delete pmtSizeCmd;
  delete photocathodeSizeCmd;
	delete pmtOptCouplingSizeCmd;
  delete detectorDir;
  delete volumesDir;
  delete wlsCmd;
	delete readoutModeCmd;
	delete pmtUseRealQECmd;

	delete stripNumberCmd;
	delete planeNumberCmd;
	delete planeDistanceCmd;
	delete planeXYDistanceCmd;
	delete planeTiltAngleCmd;
	delete superplaneThicknessCmd;
	delete superplaneSizeZCmd;
	delete yPlaneInSuperModuleCmd;

	delete worldMaterialCmd;
	delete superplaneCasingMaterialCmd;
	delete scintillationYieldCmd;
	delete scintillationWeightCmd;
	
 
}


void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == stripSizeCmd ){ 
    Detector->SetStripSize(stripSizeCmd->GetNew3VectorValue(newValue));
  }
	else if( command == grooveSizeCmd ){ 
    Detector->SetGrooveSize(grooveSizeCmd->GetNew3VectorValue(newValue));
  }
	else if (command == stripDesignModeCmd){
    Detector->SetStripDesignMode(stripDesignModeCmd->GetNewIntValue(newValue));
  }
	else if (command == stripReadoutModeCmd){
    Detector->SetStripReadoutMode(stripReadoutModeCmd->GetNewIntValue(newValue));
  }
	else if (command == stripOptCouplingModeCmd){
    Detector->SetStripOptCouplingMode(stripOptCouplingModeCmd->GetNewIntValue(newValue));
  }
	else if (command == fiberLengthCmd){
    Detector->SetFiberSizeZ(fiberLengthCmd->GetNewDoubleValue(newValue));
  }
	else if (command == fiberRadiusCmd){
    Detector->SetFiberRadius(fiberRadiusCmd->GetNewDoubleValue(newValue));
  }
	else if (command == fiberAbsLengthCmd){
    Detector->SetFiberAbsorptionLength(fiberAbsLengthCmd->GetNewDoubleValue(newValue));
  }
  else if (command == housingThicknessCmd){
    Detector->SetHousingThickness(housingThicknessCmd->GetNewDoubleValue(newValue));
  }
	else if (command == housingReflectivityCmd){
    Detector->SetHousingReflectivity(housingReflectivityCmd->GetNewDoubleValue(newValue));
  }
  else if (command == housingUseRealReflectivityCmd){
    Detector->SetRealReflectivity(housingUseRealReflectivityCmd->GetNewBoolValue(newValue));
  }
  else if (command == housingSetSurfaceTypeCmd){
		Detector->SetHousingSurfaceType(housingSetSurfaceTypeCmd->GetNewIntValue(newValue));
  }
	
  else if (command == pmtSizeCmd){
    Detector->SetPMTSize(pmtSizeCmd->GetNew3VectorValue(newValue));
  }
	else if (command == photocathodeSizeCmd){
    Detector->SetPhotocathodeSize(photocathodeSizeCmd->GetNew3VectorValue(newValue));
  }
	else if (command == pmtOptCouplingSizeCmd){
    Detector->SetPMTOpticalCouplingSize(pmtOptCouplingSizeCmd->GetNew3VectorValue(newValue));
  }
 
  else if (command == wlsCmd){
    Detector->SetFiber(wlsCmd->GetNewBoolValue(newValue));
  }
	else if (command == readoutModeCmd){
    Detector->SetDoubleReadout(readoutModeCmd->GetNewBoolValue(newValue));
  }
	else if (command == pmtUseRealQECmd){
    Detector->SetRealPhotocathode(pmtUseRealQECmd->GetNewBoolValue(newValue));
  }
	

	else if (command == updateCmd){
    Detector->UpdateGeometry();
  }
	else if (command == stripNumberCmd){
    Detector->SetNumberOfStripsInPlane(stripNumberCmd->GetNewIntValue(newValue));
  }
	else if (command == planeNumberCmd){
    Detector->SetNumberOfPlanes(planeNumberCmd->GetNewIntValue(newValue));
  }
	else if (command == planeXYDistanceCmd){
    Detector->SetDistanceAmongXYPlanes(planeXYDistanceCmd->GetNewDoubleValue(newValue));
  }
	else if (command == planeTiltAngleCmd){
    Detector->SetPlaneTiltAngle(planeTiltAngleCmd->GetNewDoubleValue(newValue));
  }
	else if (command == superplaneThicknessCmd){
    Detector->SetSuperPlaneModuleThickness(superplaneThicknessCmd->GetNewDoubleValue(newValue));
  }
	else if (command == superplaneSizeZCmd){
    Detector->SetSuperPlaneModuleSizeZ(superplaneSizeZCmd->GetNewDoubleValue(newValue));
  }
	else if (command == scintillationYieldCmd){
    Detector->SetScintillationYield(scintillationYieldCmd->GetNewDoubleValue(newValue));
  }
	else if (command == scintillationWeightCmd){
    Detector->SetScintillationWeight(scintillationWeightCmd->GetNewDoubleValue(newValue));
  }

	else if (command == worldMaterialCmd){
		Detector->SetWorldMaterialName(newValue);    
  }
	else if (command == superplaneCasingMaterialCmd){
		Detector->SetSuperPlaneMaterialName(newValue);    
  }	
	else if (command == yPlaneInSuperModuleCmd){
		Detector->UseYPlaneInSuperModule(yPlaneInSuperModuleCmd->GetNewBoolValue(newValue));    
  }	

	else if (command == planeDistanceCmd){
		G4double distance= planeDistanceCmd->GetNewDoubleValue(newValue);
		argArray.push_back(distance);
    //Detector->SetDistanceAmongPlanes(planeDistanceCmd->GetNewDoubleValue(newValue));
  }


}//close function



void DetectorMessenger::SetFinalValues(){
	
	G4cout<<"Calling DetectorMessenger::SetFinalValues()"<<G4endl;

	//G4cout<<"argArray size="<<argArray.size()<<G4endl;
	//G4cout<<" size="<<argArray.size()<<"  "<<Detector->GetNumberOfPlanes()<<G4endl;
	//for(unsigned int i=0;i<argArray.size();i++) G4cout<<"argArray["<<i<<"]="<<argArray[i]/cm<<G4endl;

	Detector->SetDistanceAmongPlanes(argArray);
}//close function


