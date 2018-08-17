/**
* @file G4MuonCounterConstruction.cc
* @class G4MuonCounterConstruction
* @brief Build the detector geometry, material & material properties.
*
* The detector is made up of a given number of detector planes, each composed by a X- and Y- plane.
* Each plane hosts a given number of strip modules. A single strip module is a scintillator bar, with a groove hosting 
* a WLS fiber. At both strip edges, a PMT is placed for the readout.  
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4ScintillatorSD.hh"
#include "G4SDManager.hh"

#include "G4SuperPlaneModule.hh"
#include "G4PlaneModule.hh"
#include "G4StripModule.hh"
#include "G4Tank.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4NistManager.hh"

//## offline headers
#include <fwk/CentralConfig.h>
#include <fwk/RandomEngineRegistry.h>
#include <fwk/RunController.h>
#include <fwk/CentralConfig.h>

#include <utl/ErrorLogger.h>
#include <utl/CoordinateSystem.h>
#include <utl/Point.h>
#include <utl/TabulatedFunction.h>
#include <utl/Reader.h>

#include <TMath.h>

#include <vector>
#include <cstddef>
#include <iostream>
#include <sstream>

using utl::CoordinateSystem;
using utl::Point;
using utl::TabulatedFunction;
using utl::Branch;
using utl::ErrorLogger;
using namespace fwk;
using namespace std;
using namespace G4MuonCounterSimulatorUSC;

double G4MuonCounterConstruction::fStripVetoTime;
double G4MuonCounterConstruction::fStripEnergyThreshold;
double G4MuonCounterConstruction::fStripDecayTime;
double G4MuonCounterConstruction::fFiberDecayTime;
double G4MuonCounterConstruction::fFiberCriticalAngle;
double G4MuonCounterConstruction::fFiberCoreRefractiveIndex;
double G4MuonCounterConstruction::fFiberClad1RefractiveIndex;
double G4MuonCounterConstruction::fFiberClad2RefractiveIndex;
double G4MuonCounterConstruction::fPMTTransitTimeAverage;
double G4MuonCounterConstruction::fPMTTransitTimeSpread;
double G4MuonCounterConstruction::fPMTPhotoelectronYield;


G4MuonCounterConstruction::G4MuonCounterConstruction()
: H(0), B(0), C(0), N(0), O(0), Na(0), Al(0), Fe(0), Mn(0), Mg(0), Ca(0), K(0), P(0), LXe(0),Air(0),Vacuum(0),Glass(0),
  Polystyrene(0), PMMA(0),Pethylene(0),fPethylene(0),PolyvinylToluene(0),
  TiO2(0),Polystyrene_TiO2mixed(0),Polystyrene_doped(0),MalargueSoil(0),Silicone(0),SiO2(0),B2O3(0),Na2O(0),Al2O3(0),Fe2O3(0),MnO(0),MgO(0),CaO(0),K2O(0),P2O5(0),Borosilicate(0),PPO(0),POPOP(0),PVC(0),Water(0),VacuumPT(0),AirPT(0),ScintillatorPT(0),IdealScintCoatingPT(0),RealScintCoatingPT(0),WLSFiberPT(0),
  WLSClad1PT(0),WLSClad2PT(0),PMTPT(0),OpticalCouplingPT(0),IdealPhotocathodePT(0),RealPhotocathodePT(0),expHall_box(0),expHall_log(0),expHall_phys(0),expHall_va(0)
{

	//create a messenger for this class
  detMessenger = new DetectorMessenger(this);//create a messenger for this class
  
	//## Define standard parameters
	SetParameters();

	//## get parameters from XML config file
	//## overriding previous settings
  // SetXMLParameters();

}//close constructor


G4MuonCounterConstruction::~G4MuonCounterConstruction(){

	delete detMessenger;
}


void G4MuonCounterConstruction::SetParameters(){

	Nstrips= 1;
	Nplanes= 1;
	Nstrips_tot= Nplanes*Nstrips;

	//## STRIP DEFAULT OPTIONS
	fStripDesignMode= kScintPlusFiber;
	fStripReadoutMode= kDoubleReadout;
	fStripOptCouplingMode= kNoCoupling;
	fFiberOptCouplingMode= kNoCoupling;

 	//## STRIP 
  Xstrip= 1.0*CLHEP::m;
  Ystrip= 1.0*CLHEP::cm;
  Zstrip= 1.0*CLHEP::cm;
  
	//## GROOVE
	Xgroove= Xstrip;
  Ygroove= 1.2*CLHEP::mm;
  Zgroove= 1.2*CLHEP::mm;

	//## SCINT COATING
	CoatingThickness= 0.12*CLHEP::mm;
	//Xcoat= Xstrip+CoatingThickness;//SINGLE READOUT 
	Xcoat= Xstrip;//DOUBLE READOUT 
	Ycoat= Ystrip+2.0*CoatingThickness;
	Zcoat= Zstrip+2.0*CoatingThickness;

	
	//## WORLD
  expHall_x = Xstrip*1.5;
  expHall_y = Nstrips*Ycoat*1.5;
  expHall_z = 3*CLHEP::m;//typical max depth for whole BATATA detector 

	//## SUPERPLANE
	fSuperPlaneModuleThickness= 3.0*CLHEP::mm;
	fSuperPlaneModuleSizeZ= 6.0*CLHEP::cm;

	//## PLANE
	for(unsigned int k=0;k<Nplanes;k++) fPlaneDistance.push_back(10.*CLHEP::cm);
	fXYPlaneDistance= Ycoat;
	fPlaneTiltAngle= 0.*CLHEP::deg;
	
	//## WLS FIBER
  fFiberRadius= 0.5*CLHEP::mm;
	
	fiber_rmin = 0.00*CLHEP::mm;
  fiber_rmax = fFiberRadius-0.03*(2*fFiberRadius)-0.03*(2*fFiberRadius);//fiber_radius - 3%*fiber_diameter - 3%*fiber_diameter; 
  fiber_z    = Xstrip;
  fiber_sphi = 0.*CLHEP::deg;
  fiber_ephi = 360.*CLHEP::deg;	

	clad1_rmin = 0.00*CLHEP::mm;// fiber_rmax;
  clad1_rmax = fFiberRadius-0.03*(2*fFiberRadius);// fiber_radius - 3%*fiber_diameter;      
  clad1_z    = Xstrip;
  clad1_sphi = fiber_sphi;
  clad1_ephi = fiber_ephi;

  clad2_rmin = 0.0*CLHEP::mm;//clad1_rmax;
  clad2_rmax = fFiberRadius;    
  clad2_z    = Xstrip;
  clad2_sphi = fiber_sphi;
  clad2_ephi = fiber_ephi;	

	fiber_captureLength= 0.1*CLHEP::mm;
	fiber_absLength= 350*CLHEP::cm;
	fFiberCoreRefractiveIndex= 1.59;
	fFiberClad1RefractiveIndex= 1.49;
	fFiberClad2RefractiveIndex= 1.42;
	

	//## PMT
  pmt_sizeY= Ycoat;
  pmt_sizeZ= Zcoat;
  pmt_length= 8.2*CLHEP::cm;
 
  photocathode_sizeY = Ycoat;
  photocathode_sizeZ = Zcoat;
  photocathode_length = 1.0*CLHEP::mm;

	pmtOptCoupling_sizeY = Ycoat;
  pmtOptCoupling_sizeZ = Zcoat;
  pmtOptCoupling_length = 0.5*CLHEP::mm;


  CoatingReflectivity= 1.0;
	UseRealReflectivity= true;
	CoatingSurfaceType= 1;

	UseRealPhotocathode= false;
	fUseFiber= false;
	fUseDoubleReadout= true;
	HasYPlaneInSuperModule= true;

	fWorldMaterialType= 1;
	fWorldMaterialName="Air";
 	fSuperPlaneMaterialName="PVC";

	fScintillationYield= 11136./CLHEP::MeV;//true yield of BC408
	fScintillationWeight= 1.0;

	fSurfaceSeparationGap= 1.0*CLHEP::nm;

}//close G4MuonCounterConstruction::SetParameters()


void G4MuonCounterConstruction::DefineMaterials(){

	//	------------- Materials -------------
  //see all material with /material/nist/listMaterials all
  G4NistManager* man = G4NistManager::Instance();

  double a, z, density;
  double fractionmass;
  int nelements;
  
  int polyPMMA = 1;
  int nC_PMMA = 3+2*polyPMMA;
  int nH_PMMA = 6+2*polyPMMA;

  int polyeth = 1;
  int nC_eth = 2*polyeth;
  int nH_eth = 4*polyeth;

  //*** Elements	
  H = new G4Element("H", "H", z=1., a=1.01*CLHEP::g/CLHEP::mole);
	B = new G4Element("B", "B", z=5., a= 10.811*CLHEP::g/CLHEP::mole);	
  C = new G4Element("C", "C", z=6., a=12.01*CLHEP::g/CLHEP::mole);
  N = new G4Element("N", "N", z=7., a= 14.01*CLHEP::g/CLHEP::mole);
  O = new G4Element("O", "O", z=8., a= 16.00*CLHEP::g/CLHEP::mole);
	Na = new G4Element("Na","Na", z=11., a= 22.99*CLHEP::g/CLHEP::mole);	
  Al = new G4Material("Al",z=13.,a=26.98*CLHEP::g/CLHEP::mole,density=2.7*CLHEP::g/CLHEP::cm3);//Aluminum
	Fe = new G4Element("Fe","Fe", z=26., a= 55.8457*CLHEP::g/CLHEP::mole);	
	Mn = new G4Element("Mn","Mn", z=25., a= 54.938049*CLHEP::g/CLHEP::mole);	
	Mg = new G4Element("Mg","Mg", z=12., a= 24.3050*CLHEP::g/CLHEP::mole);	
	Ca = new G4Element("Ca","Ca", z=20., a= 40.078*CLHEP::g/CLHEP::mole);	
	K = new G4Element("K","K", z=19., a= 39.0983*CLHEP::g/CLHEP::mole);	
	P = new G4Element("P","P", z=15., a= 30.973761*CLHEP::g/CLHEP::mole);	
	Cl= new G4Element("Cl","Cl", z=17., a= 35.453*CLHEP::g/CLHEP::mole);	

  //***Materials
	//TiO2
  TiO2 = man->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");
	//SiO2
  SiO2 = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	//Al2O3 aluminium oxide
  Al2O3 = man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	//B2O3 diboron trioxide  
	B2O3 = new G4Material("B2O3", density= 2.550*CLHEP::g/CLHEP::cm3, 2);
  B2O3->AddElement(B, 2);
  B2O3->AddElement(O, 3);

	B2O2 = new G4Material("B2O2", density= 2.23*CLHEP::g/CLHEP::cm3, 2);
  B2O2->AddElement(B, 2);
  B2O2->AddElement(O, 2);

  //Na2O sodium monoxide
	Na2O = new G4Material("Na2O", density= 2.270*CLHEP::g/CLHEP::cm3, 2);
  Na2O->AddElement(Na, 2);
  Na2O->AddElement(O, 1);
	//Fe2O3 
	Fe2O3= man->FindOrBuildMaterial("G4_FERRIC_OXIDE");
	//MnO 
	MnO = new G4Material("MnO", density= 5.37*CLHEP::g/CLHEP::cm3, 2);
  MnO->AddElement(Mn, 1);
  MnO->AddElement(O, 1);
	//MgO 
	MgO = man->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE");
	//CaO 
	CaO = man->FindOrBuildMaterial("G4_CALCIUM_OXIDE");
	//K2O
	K2O = man->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");
	//P2O5 Phosphoric anhydride
	P2O5= new G4Material("P2O5", density= 2.39*CLHEP::g/CLHEP::cm3, 2);
	P2O5->AddElement(P, 2);
  P2O5->AddElement(O, 5);

  //Liquid Xenon
  LXe = new G4Material("LXe",z=54.,a=131.29*CLHEP::g/CLHEP::mole,density=3.020*CLHEP::g/CLHEP::cm3);
  
	//Vacuum
  Vacuum = new G4Material("Vacuum",z=1.,a=1.01*CLHEP::g/CLHEP::mole,density=universe_mean_density,kStateGas,0.1*CLHEP::kelvin,1.e-19*CLHEP::pascal); 
  
  //Air
	//Air = man->FindOrBuildMaterial("G4_Air");
 	Air = new G4Material("Air", density= 1.29*CLHEP::mg/CLHEP::cm3, 2);
  Air->AddElement(N, 70*CLHEP::perCent);
  Air->AddElement(O, 30*CLHEP::perCent);

	// Water
  Water = new G4Material("Water", density= 1.0*CLHEP::g/CLHEP::cm3, nelements=2);
  Water->AddElement(H, 2);
  Water->AddElement(O, 1);

	//Malargue soil
  MalargueSoil= new G4Material("MalargueSoil", density= 1.8*CLHEP::g/CLHEP::cm3, 10);
  MalargueSoil->AddMaterial(SiO2,fractionmass=63.33*CLHEP::perCent);
  MalargueSoil->AddMaterial(Al2O3,fractionmass=15.38*CLHEP::perCent);
	MalargueSoil->AddMaterial(Fe2O3,fractionmass=6.39*CLHEP::perCent);	
	MalargueSoil->AddMaterial(CaO,fractionmass=5.81*CLHEP::perCent);	
  MalargueSoil->AddMaterial(Na2O,fractionmass=3.66*CLHEP::perCent);	
	MalargueSoil->AddMaterial(MgO,fractionmass=2.28*CLHEP::perCent);	
	MalargueSoil->AddMaterial(K2O,fractionmass=2.06*CLHEP::perCent);	
	MalargueSoil->AddMaterial(TiO2,fractionmass=0.82*CLHEP::perCent);	
	MalargueSoil->AddMaterial(P2O5,fractionmass=0.17*CLHEP::perCent);	
	MalargueSoil->AddMaterial(MnO,fractionmass=0.10*CLHEP::perCent);	
	
  //Glass
  Glass = new G4Material("Glass", density=1.032*CLHEP::g/cm3,2);
  Glass->AddElement(C,91.533*CLHEP::perCent);
  Glass->AddElement(H,8.467*CLHEP::perCent);

	
  //Polystyrene
  Polystyrene = new G4Material("Polystyrene", density= 1.05*CLHEP::g/CLHEP::cm3, 2);
  Polystyrene->AddElement(C, 8);
  Polystyrene->AddElement(H, 8);

	//PPO (C15 H11 NO)
  PPO = new G4Material("PPO", density= 1.06*CLHEP::g/CLHEP::cm3, 4);
  PPO->AddElement(C, 15);
  PPO->AddElement(H, 11);
	PPO->AddElement(N, 1);
  PPO->AddElement(O, 1);

	//POPOP (C24 H16 N2 O2)
  POPOP = new G4Material("POPOP", density= 1.204*CLHEP::g/CLHEP::cm3, 4);
  POPOP->AddElement(C, 24);
  POPOP->AddElement(H, 16);
	POPOP->AddElement(N, 2);
  POPOP->AddElement(O, 2);
  
  //Polyvinyl Toluene & organic fluors (Material of commercial scintillator like Bicron BC-408, BC-400 (equivalent to NE-202) )
	//  chemical formula C(10)H(11) 97% , <3% organic fluors
	PolyvinylToluene = new G4Material("PolyvinylToluene", density= 1.032*CLHEP::g/CLHEP::cm3, 2);
  PolyvinylToluene->AddElement(C, 10);
  PolyvinylToluene->AddElement(H, 11);

 
  // MINOS measure a strip density (including the coating of 1.046+-0.004)
  //Polystyrene mixed with TiO2: see MINOS document http://minos-docdb.fnal.gov/cgi-bin/ShowDocument?docid=2080
  Polystyrene_TiO2mixed= new G4Material("Polystyrene_TiO2mixed", density= 1.046*CLHEP::g/CLHEP::cm3, 2);
  Polystyrene_TiO2mixed->AddMaterial(Polystyrene,fractionmass=90.0*CLHEP::perCent);
  Polystyrene_TiO2mixed->AddMaterial(TiO2,fractionmass=10.0*CLHEP::perCent);

	//Polystyrene mixed with PPO & POPOP  (leave density unchanged w.r.t PS)
  Polystyrene_doped= new G4Material("Polystyrene_doped", density= 1.05*CLHEP::g/CLHEP::cm3, 3);
  Polystyrene_doped->AddMaterial(Polystyrene,fractionmass=98.97*CLHEP::perCent);
  Polystyrene_doped->AddMaterial(PPO,fractionmass=1.0*CLHEP::perCent);
	Polystyrene_doped->AddMaterial(POPOP,fractionmass=0.03*CLHEP::perCent);	

	//Silicone optical grease (see /examples/extended/wls)
	Silicone= new G4Material("Silicone", density= 1.060*CLHEP::g/CLHEP::cm3, 2);
  Silicone->AddElement(C, 2);
  Silicone->AddElement(H, 6);

	//Borosilicate
	Borosilicate= new G4Material("Borosilicate", density= 2.230*CLHEP::g/CLHEP::cm3, 4);
	Borosilicate->AddMaterial(SiO2, fractionmass= 80.7*CLHEP::perCent);
  Borosilicate->AddMaterial(B2O3, fractionmass= 13.0*CLHEP::perCent);
	Borosilicate->AddMaterial(Na2O, fractionmass=  4.0*CLHEP::perCent);
  Borosilicate->AddMaterial(Al2O3, fractionmass= 2.3*CLHEP::perCent);
	
  //Fiber(PMMA)
  PMMA = new G4Material("PMMA", density= 1.19*CLHEP::g/CLHEP::cm3,3);
  PMMA->AddElement(H,nH_PMMA);
  PMMA->AddElement(C,nC_PMMA);
  PMMA->AddElement(O,2);
  //Cladding(polyethylene)
  Pethylene = new G4Material("Pethylene", density=1.20*CLHEP::g/CLHEP::cm3,2);
  Pethylene->AddElement(H,nH_eth);
  Pethylene->AddElement(C,nC_eth);
  //Double cladding(flourinated polyethylene)
  fPethylene = new G4Material("fPethylene", density=1.43*CLHEP::g/CLHEP::cm3,2);
  fPethylene->AddElement(H,nH_eth);
  fPethylene->AddElement(C,nC_eth);

	PVC= man->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");

	//BisphenolAEpichlorohydrin CAS Number: 25068-38-6 
	// Formula: C18H21ClO3 or (C15H16O2 - C3H5ClO)x
	// see http://www.sigmaaldrich.com/catalog/ProductDetail.do?D7=0&N5=SEARCH_CONCAT_PNO|BRAND_KEY&N4=181196|ALDRICH&N25=0&QS=ON&F=SPEC
	BisphenolAEpichlorohydrin= new G4Material("BisphenolAEpichlorohydrin", density= 1.18*CLHEP::g/CLHEP::cm3, 4);
  BisphenolAEpichlorohydrin->AddElement(C, 18);
  BisphenolAEpichlorohydrin->AddElement(H, 21);
	BisphenolAEpichlorohydrin->AddElement(Cl, 1);
	BisphenolAEpichlorohydrin->AddElement(O, 3);

	// N-ButylGlycidylEther or Oxirane,2-(butoxymethyl)-  CAS Number: 2426-O8-5
	// Formula: C7H14O2
	// see http://www.lookchem.com/cas-242/2426-08-6.html
	NButylGlycidylEther= new G4Material("NButylGlycidylEther", density= 0.912*CLHEP::g/CLHEP::cm3, 3);
  NButylGlycidylEther->AddElement(C, 7);
  NButylGlycidylEther->AddElement(H, 14);
	NButylGlycidylEther->AddElement(O, 2);

	// Shell Epon 815C resin (used in MINOS detectors)
	// Epoxy composition from DataSheet
	EpoxyShell= new G4Material("EpoxyShell", density= 1.13*CLHEP::g/CLHEP::cm3, 2);
  EpoxyShell->AddMaterial(BisphenolAEpichlorohydrin,fractionmass=86.4*CLHEP::perCent);
  EpoxyShell->AddMaterial(NButylGlycidylEther,fractionmass=13.6*CLHEP::perCent);

	// BisphenolADiglycidylEther CAS Number: 1675-54-3
	// Formula: C19H20O4
	// see http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/Epotek-301-1.html
	BisphenolADiglycidylEther= new G4Material("BisphenolADiglycidylEther", density= 1.16*CLHEP::g/CLHEP::cm3, 3);
  BisphenolADiglycidylEther->AddElement(C, 19);
  BisphenolADiglycidylEther->AddElement(H, 20);
	BisphenolADiglycidylEther->AddElement(O, 4);

	// ButanediolDiglycidylEther CAS Number: 2425-79-8 
	// Formula: C10H18O4 
	// see http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/Epotek-301-1.html
	ButanediolDiglycidylEther= new G4Material("ButanediolDiglycidylEther", density= 1.10*CLHEP::g/CLHEP::cm3, 3);
  ButanediolDiglycidylEther->AddElement(C, 10);
  ButanediolDiglycidylEther->AddElement(H, 18);
	ButanediolDiglycidylEther->AddElement(O, 4);

	// HexanediamineTrimethyl CAS Number: 3236-53-1 
	// Formula: C9H22N2
	// see http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/Epotek-301-1.html
	HexanediamineTrimethyl= new G4Material("HexanediamineTrimethyl", density= 0.865*CLHEP::g/CLHEP::cm3, 3);
  HexanediamineTrimethyl->AddElement(C, 9);
  HexanediamineTrimethyl->AddElement(H, 22);
	HexanediamineTrimethyl->AddElement(N, 2);

	// EPO-TEK® 301-1 resin (used in AUGER detectors)
	// Epoxy composition from DataSheet 
	// A=11.8991, Z=5.9883, density=1.19 g/cm^3, radLength= 34.9894 cm or 41.6374 g/cm^2 
	// see http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/Epotek-301-1.html
	EpoxyTek= new G4Material("EpoxyTek", density= 1.19*CLHEP::g/CLHEP::cm3, 3);
  EpoxyTek->AddMaterial(BisphenolADiglycidylEther,fractionmass=56.0*CLHEP::perCent);
	EpoxyTek->AddMaterial(ButanediolDiglycidylEther,fractionmass=24.0*CLHEP::perCent);
  EpoxyTek->AddMaterial(HexanediamineTrimethyl,fractionmass=20.0*CLHEP::perCent);


	HDPE = new G4Material("HDPE", 0.94*CLHEP::g/CLHEP::cm3, 2);
  HDPE->AddElement(C, 2);
  HDPE->AddElement(H, 4);


	Pyrex = new G4Material("Pyrex", 2.23*CLHEP::g/CLHEP::cm3, 3);
  Pyrex->AddMaterial(SiO2, 0.80);
  Pyrex->AddMaterial(B2O2, 0.13);
  Pyrex->AddMaterial(Na2O, 0.07);

	Interface = new G4Material("Interface", 0.97*CLHEP::g/CLHEP::cm3, 3);
  Interface->AddElement(C, 4);
  Interface->AddElement(H, 8);
  Interface->AddElement(O, 2);

	Lucite = new G4Material("Lucite", 1.20*CLHEP::g/CLHEP::cm3, 3);
  Lucite->AddElement(C, 4);
  Lucite->AddElement(H, 8);
  Lucite->AddElement(O, 2);
	

//
// ------------ Generate & Add Material Properties Table ------------
//
/*
//
// Vacuum
//
  const int NUMENTRIES_Vacuum = 2;

  double PhotonEnergy_Vacuum[NUMENTRIES_Vacuum]= { 2.0*CLHEP::eV, 3.5*CLHEP::eV};
  double RIND_Vacuum[NUMENTRIES_Vacuum]={1.0, 1.0};
  
  VacuumPT = new G4MaterialPropertiesTable();
  VacuumPT->AddProperty("RINDEX", PhotonEnergy_Vacuum, RIND_Vacuum, NUMENTRIES_Vacuum);
  
  Vacuum->SetMaterialPropertiesTable(VacuumPT);
*/
//
// Air
//
  const int NUMENTRIES_Air = 2;

  double PhotonEnergy_Air[NUMENTRIES_Air]= { 2.0*CLHEP::eV, 3.5*CLHEP::eV};
  double RIND_Air[NUMENTRIES_Air]={1.0008, 1.0008};
  
  AirPT = new G4MaterialPropertiesTable();
  AirPT->AddProperty("RINDEX", PhotonEnergy_Air, RIND_Air, NUMENTRIES_Air);
  
  Air->SetMaterialPropertiesTable(AirPT);


//
// Silicone grease
//
  const int NUMENTRIES_Silicone = 2;

  double PhotonEnergy_Silicone[NUMENTRIES_Silicone]= { 2.0*CLHEP::eV, 3.5*CLHEP::eV};
  double RIND_Silicone[NUMENTRIES_Silicone]={1.466, 1.466};
  
  OpticalCouplingPT = new G4MaterialPropertiesTable();
  OpticalCouplingPT->AddProperty("RINDEX", PhotonEnergy_Silicone, RIND_Silicone, NUMENTRIES_Silicone);
  
  Silicone->SetMaterialPropertiesTable(OpticalCouplingPT);

 
//
// PolyvinylToluene 
//  
// ***   BC-408 Characteristics ***
// ***      ScintillationYield   = 64% Anthracene (17400 photons/MeV)= 11136 photons/MeV ==> SCINTILLATIONYIELD
// ***      Refraction Index     = 1.58                                                  ==> RINDEX
// ***      Emission Spectrum    = see datasheet - peak at 425 nm                        ==>
// ***      Attenuation Length   = 210 cm 
// ***      Decay Time           = 2.1 ns
// ***      Rise Time            = 0.9 ns
//
	//## General characteristics
  const int NUMENTRIES_PolyvinylToluene= 2;
	double PhotonEnergy_PolyvinylToluene[NUMENTRIES_PolyvinylToluene]={2.0*CLHEP::eV,3.5*CLHEP::eV};
  double RIND_PolyvinylToluene[NUMENTRIES_PolyvinylToluene]={1.58,1.58};
  double ABSL_PolyvinylToluene[NUMENTRIES_PolyvinylToluene]={210.*CLHEP::cm,210.*CLHEP::cm};

	//## Emission properties
	const int NUMENTRIES_PolyvinylToluene_Emiss = 141;
	const double SlowFastRatio= 0.27;
	double FastYieldRatio= 1./(1.+SlowFastRatio);

	double PhotonEnergy_PolyvinylToluene_Emiss[NUMENTRIES_PolyvinylToluene_Emiss] ={2.38248*CLHEP::eV,2.39029*CLHEP::eV,2.40233*CLHEP::eV,2.41449*CLHEP::eV,2.43727*CLHEP::eV,2.44979*CLHEP::eV,2.4644*CLHEP::eV,2.47919*CLHEP::eV,2.49014*CLHEP::eV,2.50069*CLHEP::eV,2.51387*CLHEP::eV,2.52513*CLHEP::eV,2.53598*CLHEP::eV,2.54483*CLHEP::eV,2.55637*CLHEP::eV,2.56324*CLHEP::eV,2.5771*CLHEP::eV,2.58354*CLHEP::eV,2.59327*CLHEP::eV,2.60252*CLHEP::eV,2.6069*CLHEP::eV,2.6168*CLHEP::eV,2.62123*CLHEP::eV,2.62623*CLHEP::eV,2.63348*CLHEP::eV,2.64021*CLHEP::eV,2.64754*CLHEP::eV,2.65263*CLHEP::eV,2.65946*CLHEP::eV,2.6669*CLHEP::eV,2.67438*CLHEP::eV,2.6819*CLHEP::eV,2.6918*CLHEP::eV,2.70177*CLHEP::eV,2.70944*CLHEP::eV,2.71955*CLHEP::eV,2.72672*CLHEP::eV,2.73213*CLHEP::eV,2.7424*CLHEP::eV,2.75031*CLHEP::eV,2.75827*CLHEP::eV,2.76874*CLHEP::eV,2.77369*CLHEP::eV,2.78178*CLHEP::eV,2.78992*CLHEP::eV,2.79495*CLHEP::eV,2.8*CLHEP::eV,2.80253*CLHEP::eV,2.80571*CLHEP::eV,2.80825*CLHEP::eV,2.81079*CLHEP::eV,2.81335*CLHEP::eV,2.81846*CLHEP::eV,2.82103*CLHEP::eV,2.82617*CLHEP::eV,2.82875*CLHEP::eV,2.83198*CLHEP::eV,2.83457*CLHEP::eV,2.83976*CLHEP::eV,2.84237*CLHEP::eV,2.84825*CLHEP::eV,2.85087*CLHEP::eV,2.85349*CLHEP::eV,2.85941*CLHEP::eV,2.8647*CLHEP::eV,2.87067*CLHEP::eV,2.87599*CLHEP::eV,2.87933*CLHEP::eV,2.88805*CLHEP::eV,2.89412*CLHEP::eV,2.89954*CLHEP::eV,2.90565*CLHEP::eV,2.90906*CLHEP::eV,2.91179*CLHEP::eV,2.91796*CLHEP::eV,2.92415*CLHEP::eV,2.92692*CLHEP::eV,2.93315*CLHEP::eV,2.93662*CLHEP::eV,2.9422*CLHEP::eV,2.94569*CLHEP::eV,2.94849*CLHEP::eV,2.952*CLHEP::eV,2.95482*CLHEP::eV,2.95834*CLHEP::eV,2.96117*CLHEP::eV,2.96471*CLHEP::eV,2.96755*CLHEP::eV,2.9711*CLHEP::eV,2.97395*CLHEP::eV,2.98039*CLHEP::eV,2.9811*CLHEP::eV,2.98397*CLHEP::eV,2.99045*CLHEP::eV,2.99334*CLHEP::eV,2.99406*CLHEP::eV,2.99696*CLHEP::eV,3.00058*CLHEP::eV,3.00349*CLHEP::eV,3.00713*CLHEP::eV,3.01005*CLHEP::eV,3.01371*CLHEP::eV,3.01665*CLHEP::eV,3.01665*CLHEP::eV,3.02032*CLHEP::eV,3.02327*CLHEP::eV,3.02696*CLHEP::eV,3.02696*CLHEP::eV,3.03362*CLHEP::eV,3.03659*CLHEP::eV,3.04032*CLHEP::eV,3.0433*CLHEP::eV,3.04405*CLHEP::eV,3.05079*CLHEP::eV,3.0538*CLHEP::eV,3.05756*CLHEP::eV,3.06058*CLHEP::eV,3.06436*CLHEP::eV,3.0674*CLHEP::eV,3.07119*CLHEP::eV,3.07729*CLHEP::eV,3.08112*CLHEP::eV,3.08802*CLHEP::eV,3.09496*CLHEP::eV,3.10115*CLHEP::eV,3.10815*CLHEP::eV,3.11518*CLHEP::eV,3.12224*CLHEP::eV,3.12933*CLHEP::eV,3.13963*CLHEP::eV,3.15321*CLHEP::eV,3.16771*CLHEP::eV,3.18561*CLHEP::eV,3.19959*CLHEP::eV,3.2212*CLHEP::eV,3.25503*CLHEP::eV,3.28435*CLHEP::eV,3.31154*CLHEP::eV,3.34279*CLHEP::eV,3.37005*CLHEP::eV,3.40616*CLHEP::eV};
  double SCINT_PolyvinylToluene[NUMENTRIES_PolyvinylToluene_Emiss] ={0.03188,0.03604,0.03803,0.04428,0.05466,0.05878,0.06926,0.07761,0.086,0.09227,0.1049,0.1154,0.1259,0.1344,0.1449,0.1512,0.1702,0.1787,0.1892,0.1998,0.2103,0.223,0.2315,0.24,0.2526,0.2675,0.2823,0.2993,0.3162,0.3374,0.3565,0.3756,0.3947,0.418,0.4392,0.4582,0.4752,0.49,0.507,0.5197,0.5387,0.5578,0.5727,0.5939,0.6129,0.6363,0.6618,0.6873,0.6979,0.717,0.7319,0.7489,0.7999,0.834,0.8552,0.8679,0.8871,0.9019,0.9295,0.9402,0.9486,0.9592,0.9656,0.9805,0.9889,0.9931,1.0,0.9952,0.9909,0.9844,0.9801,0.9694,0.9588,0.9503,0.9332,0.9161,0.9012,0.8905,0.8777,0.867,0.8479,0.8393,0.8244,0.8116,0.7988,0.7882,0.7775,0.7668,0.7477,0.7285,0.705,0.6944,0.6816,0.6667,0.6496,0.6369,0.6219,0.607,0.5921,0.5687,0.5453,0.5346,0.5239,0.5112,0.4962,0.4835,0.4728,0.4622,0.4472,0.4323,0.4174,0.4025,0.3876,0.3726,0.3577,0.3407,0.3279,0.3151,0.3044,0.2895,0.2767,0.2618,0.249,0.2319,0.217,0.202,0.1871,0.1721,0.1636,0.1486,0.1336,0.1229,0.11,0.09715,0.08637,0.07127,0.06045,0.05177,0.04308,0.03866,0.0342};
	
	double SCINT_PolyvinylToluene_SLOW[NUMENTRIES_PolyvinylToluene_Emiss] ={0.03188,0.03604,0.03803,0.04428,0.05466,0.05878,0.06926,0.07761,0.086,0.09227,0.1049,0.1154,0.1259,0.1344,0.1449,0.1512,0.1702,0.1787,0.1892,0.1998,0.2103,0.223,0.2315,0.24,0.2526,0.2675,0.2823,0.2993,0.3162,0.3374,0.3565,0.3756,0.3947,0.418,0.4392,0.4582,0.4752,0.49,0.507,0.5197,0.5387,0.5578,0.5727,0.5939,0.6129,0.6363,0.6618,0.6873,0.6979,0.717,0.7319,0.7489,0.7999,0.834,0.8552,0.8679,0.8871,0.9019,0.9295,0.9402,0.9486,0.9592,0.9656,0.9805,0.9889,0.9931,1.0,0.9952,0.9909,0.9844,0.9801,0.9694,0.9588,0.9503,0.9332,0.9161,0.9012,0.8905,0.8777,0.867,0.8479,0.8393,0.8244,0.8116,0.7988,0.7882,0.7775,0.7668,0.7477,0.7285,0.705,0.6944,0.6816,0.6667,0.6496,0.6369,0.6219,0.607,0.5921,0.5687,0.5453,0.5346,0.5239,0.5112,0.4962,0.4835,0.4728,0.4622,0.4472,0.4323,0.4174,0.4025,0.3876,0.3726,0.3577,0.3407,0.3279,0.3151,0.3044,0.2895,0.2767,0.2618,0.249,0.2319,0.217,0.202,0.1871,0.1721,0.1636,0.1486,0.1336,0.1229,0.11,0.09715,0.08637,0.07127,0.06045,0.05177,0.04308,0.03866,0.0342};
  for(int i=0;i<NUMENTRIES_PolyvinylToluene_Emiss;i++) {
		SCINT_PolyvinylToluene_SLOW[i]*= SlowFastRatio;
		SCINT_PolyvinylToluene[i]/= SlowFastRatio;
	}

	//## Rayleigh scattering properties
	const int NUMENTRIES_PolyvinylToluene_RaylAbs = 151;
	double PhotonEnergy_PolyvinylToluene_RaylAbs[NUMENTRIES_PolyvinylToluene_RaylAbs]={2*CLHEP::eV,2.01*CLHEP::eV,2.02*CLHEP::eV,2.03*CLHEP::eV,2.04*CLHEP::eV,2.05*CLHEP::eV,2.06*CLHEP::eV,2.07*CLHEP::eV,2.08*CLHEP::eV,2.09*CLHEP::eV,2.1*CLHEP::eV,2.11*CLHEP::eV,2.12*CLHEP::eV,2.13*CLHEP::eV,2.14*CLHEP::eV,2.15*CLHEP::eV,2.16*CLHEP::eV,2.17*CLHEP::eV,2.18*CLHEP::eV,2.19*CLHEP::eV,2.2*CLHEP::eV,2.21*CLHEP::eV,2.22*CLHEP::eV,2.23*CLHEP::eV,2.24*CLHEP::eV,2.25*CLHEP::eV,2.26*CLHEP::eV,2.27*CLHEP::eV,2.28*CLHEP::eV,2.29*CLHEP::eV,2.3*CLHEP::eV,2.31*CLHEP::eV,2.32*CLHEP::eV,2.33*CLHEP::eV,2.34*CLHEP::eV,2.35*CLHEP::eV,2.36*CLHEP::eV,2.37*CLHEP::eV,2.38*CLHEP::eV,2.39*CLHEP::eV,2.4*CLHEP::eV,2.41*CLHEP::eV,2.42*CLHEP::eV,2.43*CLHEP::eV,2.44*CLHEP::eV,2.45*CLHEP::eV,2.46*CLHEP::eV,2.47*CLHEP::eV,2.48*CLHEP::eV,2.49*CLHEP::eV,2.5*CLHEP::eV,2.51*CLHEP::eV,2.52*CLHEP::eV,2.53*CLHEP::eV,2.54*CLHEP::eV,2.55*CLHEP::eV,2.56*CLHEP::eV,2.57*CLHEP::eV,2.58*CLHEP::eV,2.59*CLHEP::eV,2.6*CLHEP::eV,2.61*CLHEP::eV,2.62*CLHEP::eV,2.63*CLHEP::eV,2.64*CLHEP::eV,2.65*CLHEP::eV,2.66*CLHEP::eV,2.67*CLHEP::eV,2.68*CLHEP::eV,2.69*CLHEP::eV,2.7*CLHEP::eV,2.71*CLHEP::eV,2.72*CLHEP::eV,2.73*CLHEP::eV,2.74*CLHEP::eV,2.75*CLHEP::eV,2.76*CLHEP::eV,2.77*CLHEP::eV,2.78*CLHEP::eV,2.79*CLHEP::eV,2.8*CLHEP::eV,2.81*CLHEP::eV,2.82*CLHEP::eV,2.83*CLHEP::eV,2.84*CLHEP::eV,2.85*CLHEP::eV,2.86*CLHEP::eV,2.87*CLHEP::eV,2.88*CLHEP::eV,2.89*CLHEP::eV,2.9*CLHEP::eV,2.91*CLHEP::eV,2.92*CLHEP::eV,2.93*CLHEP::eV,2.94*CLHEP::eV,2.95*CLHEP::eV,2.96*CLHEP::eV,2.97*CLHEP::eV,2.98*CLHEP::eV,2.99*CLHEP::eV,3*CLHEP::eV,3.01*CLHEP::eV,3.02*CLHEP::eV,3.03*CLHEP::eV,3.04*CLHEP::eV,3.05*CLHEP::eV,3.06*CLHEP::eV,3.07*CLHEP::eV,3.08*CLHEP::eV,3.09*CLHEP::eV,3.1*CLHEP::eV,3.11*CLHEP::eV,3.12*CLHEP::eV,3.13*CLHEP::eV,3.14*CLHEP::eV,3.15*CLHEP::eV,3.16*CLHEP::eV,3.17*CLHEP::eV,3.18*CLHEP::eV,3.19*CLHEP::eV,3.2*CLHEP::eV,3.21*CLHEP::eV,3.22*CLHEP::eV,3.23*CLHEP::eV,3.24*CLHEP::eV,3.25*CLHEP::eV,3.26*CLHEP::eV,3.27*CLHEP::eV,3.28*CLHEP::eV,3.29*CLHEP::eV,3.3*CLHEP::eV,3.31*CLHEP::eV,3.32*CLHEP::eV,3.33*CLHEP::eV,3.34*CLHEP::eV,3.35*CLHEP::eV,3.36*CLHEP::eV,3.37*CLHEP::eV,3.38*CLHEP::eV,3.39*CLHEP::eV,3.4*CLHEP::eV,3.41*CLHEP::eV,3.42*CLHEP::eV,3.43*CLHEP::eV,3.44*CLHEP::eV,3.45*CLHEP::eV,3.46*CLHEP::eV,3.47*CLHEP::eV,3.48*CLHEP::eV,3.49*CLHEP::eV,3.5*CLHEP::eV};
  double RAYLABSL_PolyvinylToluene[NUMENTRIES_PolyvinylToluene_RaylAbs]={535.817*CLHEP::m,525.234*CLHEP::m,514.91*CLHEP::m,504.839*CLHEP::m,495.012*CLHEP::m,485.424*CLHEP::m,476.067*CLHEP::m,466.934*CLHEP::m,458.019*CLHEP::m,449.316*CLHEP::m,440.818*CLHEP::m,432.521*CLHEP::m,424.417*CLHEP::m,416.503*CLHEP::m,408.772*CLHEP::m,401.22*CLHEP::m,393.842*CLHEP::m,386.632*CLHEP::m,379.586*CLHEP::m,372.701*CLHEP::m,365.97*CLHEP::m,359.391*CLHEP::m,352.959*CLHEP::m,346.671*CLHEP::m,340.522*CLHEP::m,334.508*CLHEP::m,328.627*CLHEP::m,322.874*CLHEP::m,317.247*CLHEP::m,311.742*CLHEP::m,306.355*CLHEP::m,301.085*CLHEP::m,295.927*CLHEP::m,290.879*CLHEP::m,285.939*CLHEP::m,281.103*CLHEP::m,276.369*CLHEP::m,271.734*CLHEP::m,267.195*CLHEP::m,262.751*CLHEP::m,258.4*CLHEP::m,254.137*CLHEP::m,249.963*CLHEP::m,245.873*CLHEP::m,241.867*CLHEP::m,237.943*CLHEP::m,234.097*CLHEP::m,230.329*CLHEP::m,226.637*CLHEP::m,223.018*CLHEP::m,219.471*CLHEP::m,215.994*CLHEP::m,212.586*CLHEP::m,209.245*CLHEP::m,205.969*CLHEP::m,202.757*CLHEP::m,199.608*CLHEP::m,196.519*CLHEP::m,193.49*CLHEP::m,190.519*CLHEP::m,187.605*CLHEP::m,184.746*CLHEP::m,181.941*CLHEP::m,179.19*CLHEP::m,176.49*CLHEP::m,173.841*CLHEP::m,171.242*CLHEP::m,168.691*CLHEP::m,166.187*CLHEP::m,163.73*CLHEP::m,161.318*CLHEP::m,158.95*CLHEP::m,156.625*CLHEP::m,154.343*CLHEP::m,152.102*CLHEP::m,149.901*CLHEP::m,147.741*CLHEP::m,145.619*CLHEP::m,143.535*CLHEP::m,141.488*CLHEP::m,139.478*CLHEP::m,137.503*CLHEP::m,135.563*CLHEP::m,133.657*CLHEP::m,131.784*CLHEP::m,129.944*CLHEP::m,128.136*CLHEP::m,126.36*CLHEP::m,124.614*CLHEP::m,122.898*CLHEP::m,121.212*CLHEP::m,119.554*CLHEP::m,117.925*CLHEP::m,116.323*CLHEP::m,114.749*CLHEP::m,113.201*CLHEP::m,111.679*CLHEP::m,110.182*CLHEP::m,108.711*CLHEP::m,107.264*CLHEP::m,105.84*CLHEP::m,104.441*CLHEP::m,103.064*CLHEP::m,101.711*CLHEP::m,100.379*CLHEP::m,99.0689*CLHEP::m,97.7802*CLHEP::m,96.5124*CLHEP::m,95.2651*CLHEP::m,94.0379*CLHEP::m,92.8303*CLHEP::m,91.6421*CLHEP::m,90.4729*CLHEP::m,89.3222*CLHEP::m,88.1898*CLHEP::m,87.0752*CLHEP::m,85.9782*CLHEP::m,84.8984*CLHEP::m,83.8356*CLHEP::m,82.7893*CLHEP::m,81.7592*CLHEP::m,80.7452*CLHEP::m,79.7468*CLHEP::m,78.7638*CLHEP::m,77.7959*CLHEP::m,76.8428*CLHEP::m,75.9043*CLHEP::m,74.98*CLHEP::m,74.0698*CLHEP::m,73.1734*CLHEP::m,72.2905*CLHEP::m,71.4208*CLHEP::m,70.5642*CLHEP::m,69.7204*CLHEP::m,68.8892*CLHEP::m,68.0703*CLHEP::m,67.2635*CLHEP::m,66.4687*CLHEP::m,65.6856*CLHEP::m,64.9139*CLHEP::m,64.1536*CLHEP::m,63.4044*CLHEP::m,62.666*CLHEP::m,61.9384*CLHEP::m,61.2214*CLHEP::m,60.5146*CLHEP::m,59.8181*CLHEP::m,59.1315*CLHEP::m,58.4547*CLHEP::m,57.7876*CLHEP::m,57.13*CLHEP::m};
 
	//## Create property table for scintillator and add properties
  ScintillatorPT = new G4MaterialPropertiesTable();
  ScintillatorPT->AddProperty("RINDEX",PhotonEnergy_PolyvinylToluene,RIND_PolyvinylToluene,NUMENTRIES_PolyvinylToluene);
  ScintillatorPT->AddProperty("ABSLENGTH",PhotonEnergy_PolyvinylToluene,ABSL_PolyvinylToluene,NUMENTRIES_PolyvinylToluene);
  ScintillatorPT->AddProperty("FASTCOMPONENT",PhotonEnergy_PolyvinylToluene_Emiss, SCINT_PolyvinylToluene,NUMENTRIES_PolyvinylToluene_Emiss);
	ScintillatorPT->AddProperty("SLOWCOMPONENT",PhotonEnergy_PolyvinylToluene_Emiss, SCINT_PolyvinylToluene_SLOW,NUMENTRIES_PolyvinylToluene_Emiss);

  ScintillatorPT->AddConstProperty("SCINTILLATIONYIELD",11136.0/CLHEP::MeV);//fScintillationYield/MeV eventually reduce yield by a given factor, weight final results (accuracy loss but simulation speed up)
  ScintillatorPT->AddConstProperty("RESOLUTIONSCALE",1.0);//?
  ScintillatorPT->AddConstProperty("FASTTIMECONSTANT", 2.1*CLHEP::ns);//2.1 ns
	ScintillatorPT->AddConstProperty("SLOWTIMECONSTANT", 14.2*CLHEP::ns);//14.2 ns
 	//ScintillatorPT->AddConstProperty("YIELDRATIO",1.0);//only one component
	ScintillatorPT->AddConstProperty("YIELDRATIO",FastYieldRatio);// YIELDRATIO=relative strength of the fast component as a fraction of total scintillation yield (slow/fast=0.27, fast/tot=0.7874)

  ScintillatorPT->AddProperty("RAYLEIGH",PhotonEnergy_PolyvinylToluene_RaylAbs,RAYLABSL_PolyvinylToluene,NUMENTRIES_PolyvinylToluene_RaylAbs);
  
  PolyvinylToluene->SetMaterialPropertiesTable(ScintillatorPT);
	

  //##################
	//### WLS FIBER
	//##################
	//
	//## Core
	//	
	const int NUMENTRIES_WLSFiber = 2;
	double PhotonEnergy_WLSFiber[NUMENTRIES_WLSFiber] ={2.0*CLHEP::eV,3.5*CLHEP::eV};
	//double RIND_WLSFiber[NUMENTRIES_WLSFiber]={1.59,1.59}; 
	double RIND_WLSFiber[NUMENTRIES_WLSFiber]={fFiberCoreRefractiveIndex,fFiberCoreRefractiveIndex}; 


	const int NUMENTRIES_WLSFiber_Emiss = 69;
  double PhotonEnergy_WLSFiber_Emiss[NUMENTRIES_WLSFiber_Emiss] ={
2.0*CLHEP::eV,2.00297*CLHEP::eV,2.01567*CLHEP::eV,2.04561*CLHEP::eV,2.07609*CLHEP::eV,2.11145*CLHEP::eV,2.14617*CLHEP::eV,2.16945*CLHEP::eV,2.1886*CLHEP::eV,2.20181*CLHEP::eV,2.21915*CLHEP::eV,2.24406*CLHEP::eV,2.25098*CLHEP::eV,2.2808*CLHEP::eV,2.314*CLHEP::eV,2.33624*CLHEP::eV,2.35577*CLHEP::eV,2.37109*CLHEP::eV,2.38615*CLHEP::eV,2.40186*CLHEP::eV,2.41496*CLHEP::eV,2.43058*CLHEP::eV,2.43632*CLHEP::eV,2.44689*CLHEP::eV,2.46293*CLHEP::eV,2.46882*CLHEP::eV,2.47671*CLHEP::eV,2.48266*CLHEP::eV,2.48814*CLHEP::eV,2.49666*CLHEP::eV,2.50473*CLHEP::eV,2.51031*CLHEP::eV,2.52462*CLHEP::eV,2.53339*CLHEP::eV,2.54483*CLHEP::eV,2.55954*CLHEP::eV,2.5659*CLHEP::eV,2.56855*CLHEP::eV,2.57763*CLHEP::eV,2.58354*CLHEP::eV,2.58948*CLHEP::eV,2.59273*CLHEP::eV,2.60198*CLHEP::eV,2.60471*CLHEP::eV,2.608*CLHEP::eV,2.61074*CLHEP::eV,2.61404*CLHEP::eV,2.62012*CLHEP::eV,2.62623*CLHEP::eV,2.63292*CLHEP::eV,2.63908*CLHEP::eV,2.64528*CLHEP::eV,2.65491*CLHEP::eV,2.66403*CLHEP::eV,2.67092*CLHEP::eV,2.68016*CLHEP::eV,2.68655*CLHEP::eV,2.69355*CLHEP::eV,2.70295*CLHEP::eV,2.71004*CLHEP::eV,2.72973*CLHEP::eV,2.74666*CLHEP::eV,2.78116*CLHEP::eV,2.8019*CLHEP::eV,2.8236*CLHEP::eV,2.8489*CLHEP::eV,2.85*CLHEP::eV,3.40616*CLHEP::eV,3.5*CLHEP::eV};
  double WLSFiber_EMISS[NUMENTRIES_WLSFiber_Emiss] ={
0.,0.009586,0.01278,0.02237,0.02876,0.04153,0.05431,0.07029,0.09265,0.115,0.1278,0.1725,0.1885,0.2428,0.2907,0.3419,0.4026,0.4856,0.5687,0.6454,0.7125,0.7668,0.7764,0.7859,0.7668,0.7412,0.722,0.7093,0.6965,0.6869,0.6805,0.6837,0.7093,0.7508,0.7987,0.869,0.9073,0.9393,0.9681,0.9776,0.9872,0.9936,0.9936,0.984,0.9681,0.9393,0.9042,0.853,0.8115,0.7604,0.6965,0.6102,0.4984,0.3866,0.3099,0.2396,0.1821,0.1278,0.09585,0.07348,0.03834,0.02237,0.009586,0.006391,0.006391,0.006391,0.,0.,0.};

	const int NUMENTRIES_WLSFiber_Abs = 4;
  double PhotonEnergy_WLSFiber_Abs[NUMENTRIES_WLSFiber_Abs]={2.0*CLHEP::eV,2.65*CLHEP::eV,2.7*CLHEP::eV,3.5*CLHEP::eV};
	double WLSFiber_ABSL[NUMENTRIES_WLSFiber_Abs]={3.5*CLHEP::m,3.5*CLHEP::m,4.8*CLHEP::mm,4.8*CLHEP::mm};//fiber_absLength

	
	const int NUMENTRIES_WLSFiber_Rayleigh= 151;
	double PhotonEnergy_WLSFiber_Rayleigh[NUMENTRIES_WLSFiber_Rayleigh]={2.0*CLHEP::eV,2.01*CLHEP::eV,2.02*CLHEP::eV,2.03*CLHEP::eV,2.04*CLHEP::eV,2.05*CLHEP::eV,2.06*CLHEP::eV,2.07*CLHEP::eV,2.08*CLHEP::eV,2.09*CLHEP::eV,2.1*CLHEP::eV,2.11*CLHEP::eV,2.12*CLHEP::eV,2.13*CLHEP::eV,2.14*CLHEP::eV,2.15*CLHEP::eV,2.16*CLHEP::eV,2.17*CLHEP::eV,2.18*CLHEP::eV,2.19*CLHEP::eV,2.2*CLHEP::eV,2.21*CLHEP::eV,2.22*CLHEP::eV,2.23*CLHEP::eV,2.24*CLHEP::eV,2.25*CLHEP::eV,2.26*CLHEP::eV,2.27*CLHEP::eV,2.28*CLHEP::eV,2.29*CLHEP::eV,2.3*CLHEP::eV,2.31*CLHEP::eV,2.32*CLHEP::eV,2.33*CLHEP::eV,2.34*CLHEP::eV,2.35*CLHEP::eV,2.36*CLHEP::eV,2.37*CLHEP::eV,2.38*CLHEP::eV,2.39*CLHEP::eV,2.4*CLHEP::eV,2.41*CLHEP::eV,2.42*CLHEP::eV,2.43*CLHEP::eV,2.44*CLHEP::eV,2.45*CLHEP::eV,2.46*CLHEP::eV,2.47*CLHEP::eV,2.48*CLHEP::eV,2.49*CLHEP::eV,2.5*CLHEP::eV,2.51*CLHEP::eV,2.52*CLHEP::eV,2.53*CLHEP::eV,2.54*CLHEP::eV,2.55*CLHEP::eV,2.56*CLHEP::eV,2.57*CLHEP::eV,2.58*CLHEP::eV,2.59*CLHEP::eV,2.6*CLHEP::eV,2.61*CLHEP::eV,2.62*CLHEP::eV,2.63*CLHEP::eV,2.64*CLHEP::eV,2.65*CLHEP::eV,2.66*CLHEP::eV,2.67*CLHEP::eV,2.68*CLHEP::eV,2.69*CLHEP::eV,2.7*CLHEP::eV,2.71*CLHEP::eV,2.72*CLHEP::eV,2.73*CLHEP::eV,2.74*CLHEP::eV,2.75*CLHEP::eV,2.76*CLHEP::eV,2.77*CLHEP::eV,2.78*CLHEP::eV,2.79*CLHEP::eV,2.8*CLHEP::eV,2.81*CLHEP::eV,2.82*CLHEP::eV,2.83*CLHEP::eV,2.84*CLHEP::eV,2.85*CLHEP::eV,2.86*CLHEP::eV,2.87*CLHEP::eV,2.88*CLHEP::eV,2.89*CLHEP::eV,2.9*CLHEP::eV,2.91*CLHEP::eV,2.92*CLHEP::eV,2.93*CLHEP::eV,2.94*CLHEP::eV,2.95*CLHEP::eV,2.96*CLHEP::eV,2.97*CLHEP::eV,2.98*CLHEP::eV,2.99*CLHEP::eV,3*CLHEP::eV,3.01*CLHEP::eV,3.02*CLHEP::eV,3.03*CLHEP::eV,3.04*CLHEP::eV,3.05*CLHEP::eV,3.06*CLHEP::eV,3.07*CLHEP::eV,3.08*CLHEP::eV,3.09*CLHEP::eV,3.1*CLHEP::eV,3.11*CLHEP::eV,3.12*CLHEP::eV,3.13*CLHEP::eV,3.14*CLHEP::eV,3.15*CLHEP::eV,3.16*CLHEP::eV,3.17*CLHEP::eV,3.18*CLHEP::eV,3.19*CLHEP::eV,3.2*CLHEP::eV,3.21*CLHEP::eV,3.22*CLHEP::eV,3.23*CLHEP::eV,3.24*CLHEP::eV,3.25*CLHEP::eV,3.26*CLHEP::eV,3.27*CLHEP::eV,3.28*CLHEP::eV,3.29*CLHEP::eV,3.3*CLHEP::eV,3.31*CLHEP::eV,3.32*CLHEP::eV,3.33*CLHEP::eV,3.34*CLHEP::eV,3.35*CLHEP::eV,3.36*CLHEP::eV,3.37*CLHEP::eV,3.38*CLHEP::eV,3.39*CLHEP::eV,3.4*CLHEP::eV,3.41*CLHEP::eV,3.42*CLHEP::eV,3.43*CLHEP::eV,3.44*CLHEP::eV,3.45*CLHEP::eV,3.46*CLHEP::eV,3.47*CLHEP::eV,3.48*CLHEP::eV,3.49*CLHEP::eV,3.5*CLHEP::eV};
	double WLSFiber_ABSL_Rayleigh[NUMENTRIES_WLSFiber_Rayleigh]={42.7808*CLHEP::m,41.9358*CLHEP::m,41.1115*CLHEP::m,40.3074*CLHEP::m,39.5229*CLHEP::m,38.7573*CLHEP::m,38.0102*CLHEP::m,37.281*CLHEP::m,36.5692*CLHEP::m,35.8744*CLHEP::m,35.1959*CLHEP::m,34.5334*CLHEP::m,33.8864*CLHEP::m,33.2545*CLHEP::m,32.6373*CLHEP::m,32.0343*CLHEP::m,31.4452*CLHEP::m,30.8696*CLHEP::m,30.307*CLHEP::m,29.7573*CLHEP::m,29.2199*CLHEP::m,28.6946*CLHEP::m,28.1811*CLHEP::m,27.679*CLHEP::m,27.188*CLHEP::m,26.7079*CLHEP::m,26.2383*CLHEP::m,25.779*CLHEP::m,25.3297*CLHEP::m,24.8901*CLHEP::m,24.4601*CLHEP::m,24.0393*CLHEP::m,23.6275*CLHEP::m,23.2245*CLHEP::m,22.83*CLHEP::m,22.4439*CLHEP::m,22.0659*CLHEP::m,21.6958*CLHEP::m,21.3335*CLHEP::m,20.9787*CLHEP::m,20.6312*CLHEP::m,20.2909*CLHEP::m,19.9576*CLHEP::m,19.6311*CLHEP::m,19.3112*CLHEP::m,18.9979*CLHEP::m,18.6908*CLHEP::m,18.39*CLHEP::m,18.0952*CLHEP::m,17.8062*CLHEP::m,17.523*CLHEP::m,17.2454*CLHEP::m,16.9733*CLHEP::m,16.7066*CLHEP::m,16.445*CLHEP::m,16.1886*CLHEP::m,15.9371*CLHEP::m,15.6905*CLHEP::m,15.4486*CLHEP::m,15.2114*CLHEP::m,14.9788*CLHEP::m,14.7505*CLHEP::m,14.5266*CLHEP::m,14.3069*CLHEP::m,14.0914*CLHEP::m,13.8799*CLHEP::m,13.6723*CLHEP::m,13.4687*CLHEP::m,13.2688*CLHEP::m,13.0725*CLHEP::m,12.88*CLHEP::m,12.6909*CLHEP::m,12.5053*CLHEP::m,12.3231*CLHEP::m,12.1441*CLHEP::m,11.9685*CLHEP::m,11.796*CLHEP::m,11.6265*CLHEP::m,11.4601*CLHEP::m,11.2967*CLHEP::m,11.1362*CLHEP::m,10.9785*CLHEP::m,10.8236*CLHEP::m,10.6715*CLHEP::m,10.5219*CLHEP::m,10.375*CLHEP::m,10.2307*CLHEP::m,10.0889*CLHEP::m,9.94946*CLHEP::m,9.81246*CLHEP::m,9.67781*CLHEP::m,9.54547*CLHEP::m,9.41538*CLHEP::m,9.2875*CLHEP::m,9.16178*CLHEP::m,9.03818*CLHEP::m,8.91666*CLHEP::m,8.79718*CLHEP::m,8.67969*CLHEP::m,8.56416*CLHEP::m,8.45054*CLHEP::m,8.3388*CLHEP::m,8.22889*CLHEP::m,8.1208*CLHEP::m,8.01447*CLHEP::m,7.90988*CLHEP::m,7.80699*CLHEP::m,7.70577*CLHEP::m,7.60618*CLHEP::m,7.50819*CLHEP::m,7.41178*CLHEP::m,7.31691*CLHEP::m,7.22355*CLHEP::m,7.13168*CLHEP::m,7.04127*CLHEP::m,6.95228*CLHEP::m,6.86469*CLHEP::m,6.77848*CLHEP::m,6.69362*CLHEP::m,6.61008*CLHEP::m,6.52784*CLHEP::m,6.44687*CLHEP::m,6.36716*CLHEP::m,6.28868*CLHEP::m,6.2114*CLHEP::m,6.1353*CLHEP::m,6.06037*CLHEP::m,5.98657*CLHEP::m,5.9139*CLHEP::m,5.84232*CLHEP::m,5.77183*CLHEP::m,5.7024*CLHEP::m,5.634*CLHEP::m,5.56663*CLHEP::m,5.50026*CLHEP::m,5.43488*CLHEP::m,5.37047*CLHEP::m,5.30701*CLHEP::m,5.24448*CLHEP::m,5.18287*CLHEP::m,5.12217*CLHEP::m,5.06235*CLHEP::m,5.0034*CLHEP::m,4.9453*CLHEP::m,4.88805*CLHEP::m,4.83162*CLHEP::m,4.77601*CLHEP::m,4.72119*CLHEP::m,4.66716*CLHEP::m,4.61389*CLHEP::m,4.56139*CLHEP::m};

  WLSFiberPT = new G4MaterialPropertiesTable();
  WLSFiberPT->AddProperty("RINDEX",PhotonEnergy_WLSFiber,RIND_WLSFiber,NUMENTRIES_WLSFiber);
	WLSFiberPT->AddProperty("WLSABSLENGTH",PhotonEnergy_WLSFiber_Abs,WLSFiber_ABSL,NUMENTRIES_WLSFiber_Abs);
  WLSFiberPT->AddProperty("WLSCOMPONENT",PhotonEnergy_WLSFiber_Emiss,WLSFiber_EMISS,NUMENTRIES_WLSFiber_Emiss);
  WLSFiberPT->AddConstProperty("WLSTIMECONSTANT", 7.0*CLHEP::ns);
	WLSFiberPT->AddProperty("RAYLEIGH",PhotonEnergy_WLSFiber_Rayleigh,WLSFiber_ABSL_Rayleigh,NUMENTRIES_WLSFiber_Rayleigh);
  Polystyrene->SetMaterialPropertiesTable(WLSFiberPT);


	//
	//## Cladding 1
	//
  //double RIND_WLSClad1[NUMENTRIES_WLSFiber]={1.49, 1.49};
	double RIND_WLSClad1[NUMENTRIES_WLSFiber]={fFiberClad1RefractiveIndex, fFiberClad1RefractiveIndex};
	double WLSClad1_ABSL_Rayleigh[NUMENTRIES_WLSFiber_Rayleigh]={293.4*CLHEP::m,287.605*CLHEP::m,281.952*CLHEP::m,276.437*CLHEP::m,271.056*CLHEP::m,265.806*CLHEP::m,260.682*CLHEP::m,255.681*CLHEP::m,250.799*CLHEP::m,246.034*CLHEP::m,241.381*CLHEP::m,236.837*CLHEP::m,232.4*CLHEP::m,228.067*CLHEP::m,223.833*CLHEP::m,219.698*CLHEP::m,215.658*CLHEP::m,211.71*CLHEP::m,207.852*CLHEP::m,204.081*CLHEP::m,200.396*CLHEP::m,196.794*CLHEP::m,193.272*CLHEP::m,189.828*CLHEP::m,186.461*CLHEP::m,183.168*CLHEP::m,179.948*CLHEP::m,176.798*CLHEP::m,173.716*CLHEP::m,170.702*CLHEP::m,167.752*CLHEP::m,164.866*CLHEP::m,162.042*CLHEP::m,159.278*CLHEP::m,156.573*CLHEP::m,153.925*CLHEP::m,151.332*CLHEP::m,148.794*CLHEP::m,146.309*CLHEP::m,143.876*CLHEP::m,141.493*CLHEP::m,139.159*CLHEP::m,136.873*CLHEP::m,134.634*CLHEP::m,132.44*CLHEP::m,130.291*CLHEP::m,128.186*CLHEP::m,126.122*CLHEP::m,124.1*CLHEP::m,122.119*CLHEP::m,120.177*CLHEP::m,118.273*CLHEP::m,116.407*CLHEP::m,114.577*CLHEP::m,112.783*CLHEP::m,111.025*CLHEP::m,109.3*CLHEP::m,107.609*CLHEP::m,105.95*CLHEP::m,104.323*CLHEP::m,102.727*CLHEP::m,101.162*CLHEP::m,99.6265*CLHEP::m,98.1199*CLHEP::m,96.6416*CLHEP::m,95.1911*CLHEP::m,93.7677*CLHEP::m,92.3708*CLHEP::m,90.9999*CLHEP::m,89.6542*CLHEP::m,88.3334*CLHEP::m,87.0368*CLHEP::m,85.7639*CLHEP::m,84.5141*CLHEP::m,83.2871*CLHEP::m,82.0822*CLHEP::m,80.8991*CLHEP::m,79.7372*CLHEP::m,78.5961*CLHEP::m,77.4753*CLHEP::m,76.3744*CLHEP::m,75.293*CLHEP::m,74.2307*CLHEP::m,73.187*CLHEP::m,72.1617*CLHEP::m,71.1542*CLHEP::m,70.1642*CLHEP::m,69.1914*CLHEP::m,68.2354*CLHEP::m,67.2959*CLHEP::m,66.3725*CLHEP::m,65.4648*CLHEP::m,64.5726*CLHEP::m,63.6956*CLHEP::m,62.8334*CLHEP::m,61.9857*CLHEP::m,61.1523*CLHEP::m,60.3329*CLHEP::m,59.5271*CLHEP::m,58.7348*CLHEP::m,57.9555*CLHEP::m,57.1892*CLHEP::m,56.4355*CLHEP::m,55.6941*CLHEP::m,54.9649*CLHEP::m,54.2476*CLHEP::m,53.542*CLHEP::m,52.8477*CLHEP::m,52.1647*CLHEP::m,51.4927*CLHEP::m,50.8315*CLHEP::m,50.1809*CLHEP::m,49.5406*CLHEP::m,48.9106*CLHEP::m,48.2905*CLHEP::m,47.6802*CLHEP::m,47.0795*CLHEP::m,46.4882*CLHEP::m,45.9062*CLHEP::m,45.3333*CLHEP::m,44.7693*CLHEP::m,44.214*CLHEP::m,43.6673*CLHEP::m,43.129*CLHEP::m,42.599*CLHEP::m,42.0772*CLHEP::m,41.5633*CLHEP::m,41.0572*CLHEP::m,40.5587*CLHEP::m,40.0679*CLHEP::m,39.5844*CLHEP::m,39.1082*CLHEP::m,38.6392*CLHEP::m,38.1771*CLHEP::m,37.7219*CLHEP::m,37.2735*CLHEP::m,36.8318*CLHEP::m,36.3966*CLHEP::m,35.9677*CLHEP::m,35.5452*CLHEP::m,35.1289*CLHEP::m,34.7186*CLHEP::m,34.3143*CLHEP::m,33.9159*CLHEP::m,33.5233*CLHEP::m,33.1363*CLHEP::m,32.7548*CLHEP::m,32.3789*CLHEP::m,32.0083*CLHEP::m,31.643*CLHEP::m,31.283*CLHEP::m};

  WLSClad1PT = new G4MaterialPropertiesTable();
  WLSClad1PT->AddProperty("RINDEX",PhotonEnergy_WLSFiber,RIND_WLSClad1,NUMENTRIES_WLSFiber);
  WLSClad1PT->AddProperty("ABSLENGTH",PhotonEnergy_WLSFiber_Abs,WLSFiber_ABSL,NUMENTRIES_WLSFiber_Abs);
	WLSClad1PT->AddProperty("RAYLEIGH",PhotonEnergy_WLSFiber_Rayleigh,WLSClad1_ABSL_Rayleigh,NUMENTRIES_WLSFiber_Rayleigh);
  PMMA->SetMaterialPropertiesTable(WLSClad1PT);


	//
	//## Cladding 2
	//
  //double RIND_WLSClad2[NUMENTRIES_WLSFiber]={ 1.42, 1.42};
	double RIND_WLSClad2[NUMENTRIES_WLSFiber]={fFiberClad2RefractiveIndex, fFiberClad2RefractiveIndex};
  WLSClad2PT = new G4MaterialPropertiesTable();
  WLSClad2PT->AddProperty("RINDEX",PhotonEnergy_WLSFiber,RIND_WLSClad2,NUMENTRIES_WLSFiber);
  WLSClad2PT->AddProperty("ABSLENGTH",PhotonEnergy_WLSFiber_Abs,WLSFiber_ABSL,NUMENTRIES_WLSFiber_Abs);
  fPethylene->SetMaterialPropertiesTable(WLSClad2PT);


}//close function


G4VPhysicalVolume* G4MuonCounterConstruction::Construct(){
  DefineMaterials();
  return BuildDetectorGeometry();
}//close function


G4VPhysicalVolume* G4MuonCounterConstruction::BuildDetectorGeometry(){

	cout<<"### G4MuonCounterConstruction::BuildDetectorGeometry() ### "<<endl;
		
	//re-init detector parameters
	if(fUseDetectorParametersFromXML) {
		//## get parameters from XML config file
		//## overriding previous settings
  	SetXMLParameters();
	}
	
	Xcoat= Xstrip;
	Ycoat= Ystrip+2.0*CoatingThickness;
	Zcoat= Zstrip+2.0*CoatingThickness;
	
	expHall_x = Xstrip*1.5;
  expHall_y = Nstrips*Ycoat*1.5;

	Nstrips_tot= Nplanes*Nstrips;
	
	//### WLS Fiber
  fiber_rmax = fFiberRadius-0.03*(2*fFiberRadius)-0.03*(2*fFiberRadius);//fiber_radius - 3%*fiber_diameter - 3%*fiber_diameter; 
  //fiber_z    = Xstrip;
 
  clad1_rmax = fFiberRadius-0.03*(2*fFiberRadius);// fiber_radius - 3%*fiber_diameter;      
  clad1_z    = fiber_z;
  
  clad2_rmax = fFiberRadius;    
  clad2_z    = fiber_z;
 
	fFiberCriticalAngle= TMath::Pi()/2. - asin(fFiberClad1RefractiveIndex/fFiberCoreRefractiveIndex);

  //### PMT
  if(!fUseDetectorParametersFromXML) {
		pmt_sizeY= Ycoat; 
  	pmt_sizeZ= Zcoat; 
  	photocathode_sizeY = Ycoat;
  	photocathode_sizeZ = Zcoat;
	} 
 
	//Calling messenger
	if(!fUseDetectorParametersFromXML) detMessenger->SetFinalValues();

	double TotalPlaneSizeZ= 0.*CLHEP::cm;
	for(int i=0;i<Nplanes;i++) TotalPlaneSizeZ+= fPlaneDistance[i];
	expHall_z= 2.0* TotalPlaneSizeZ;
	

//
//	------------- Volumes --------------

// The experimental Hall
//
	
	//## Draw ground layer REMOVE ME!
	//expHall_x= 6.*m;
	//expHall_y= 6.*m;
	//###


  expHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
	G4Material* worldMaterial;
	if(fWorldMaterialName=="Vacuum") {
		//cout<<"World material="<<fWorldMaterialName<<endl;	
		worldMaterial= Vacuum;
	}
	else if(fWorldMaterialName=="Air") {
		//cout<<"World material="<<fWorldMaterialName<<endl;
		worldMaterial= Air;
	}
	else if(fWorldMaterialName=="MalargueSoil") {
		//cout<<"World material="<<fWorldMaterialName<<endl;
		worldMaterial= MalargueSoil;
	}
	else{
		cerr<<"The specified material "<<fWorldMaterialName<<" has not been defined!"<<endl;
		exit(1);
	}

	expHall_log= new G4LogicalVolume(expHall_box,worldMaterial,"expHall_log",0,0,0);

  expHall_va = new G4VisAttributes(G4Colour(1,1,1));
  expHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  expHall_phys= new G4PVPlacement(0,G4ThreeVector(),expHall_log,"expHall_phys",0,false,0);


	
	//BUILD SUPER PLANES
	int planeCopyNo=0;
	double planeOffset=0.*CLHEP::cm;
	
	for(int i=0;i<Nplanes;i++){
		planeOffset+= fPlaneDistance[i];
		new G4SuperPlaneModule(0,G4ThreeVector(0.,0.,-planeOffset),expHall_log,false,planeCopyNo,this);		
		

		planeCopyNo++;
	}

	//BuildTank();


//always return the physical World
  return expHall_phys;
}


void G4MuonCounterConstruction::BuildTank(){

	double scaleFactor= 10;

	double fTankPos_x= 0.*CLHEP::m;
	double fTankPos_y= 0.*CLHEP::m;
	double fTankPos_z= 0.*CLHEP::m;
	double fTankGridSize= 4.5*CLHEP::m;
	double fTankInfillGridSize= 1.*CLHEP::m;

	double fTankHalfHeight= 0.5*1.2*CLHEP::m/scaleFactor;
  double fTankThickness= 12.7*CLHEP::mm/scaleFactor;
  
	
	//## Make tank 1
	int tankCopyNo= 0;

	double Tank1PosX= fTankPos_x;
	double Tank1PosY= fTankPos_y;
	double Tank1PosZ= fTankPos_z + fTankHalfHeight + fTankThickness;

	cout<<"G4MuonCounterConstruction::BuildTank(): INFO: Making tank no. 1..."<<endl;
	new G4Tank(0,G4ThreeVector(Tank1PosX, Tank1PosY+2./3.*fTankInfillGridSize*sin(60.*TMath::Pi()/180.),Tank1PosZ), expHall_log,false,tankCopyNo,this);		
	tankCopyNo++;

	//## Make tank 2
	double Tank2PosX= Tank1PosX+fTankGridSize/2;
	double Tank2PosY= Tank1PosY+fTankGridSize*sin(60.*TMath::Pi()/180.);
	double Tank2PosZ= Tank1PosZ;	
	cout<<"G4MuonCounterConstruction::BuildTank(): INFO: Making tank no. 2..."<<endl;
	new G4Tank(0,G4ThreeVector(Tank2PosX, Tank2PosY,Tank2PosZ), expHall_log,false,tankCopyNo,this);		
	tankCopyNo++;

	//## Make tank 3
	double Tank3PosX= Tank1PosX-fTankGridSize/2;
	double Tank3PosY= Tank1PosY+fTankGridSize*sin(60.*TMath::Pi()/180.);
	double Tank3PosZ= Tank1PosZ;
	cout<<"G4MuonCounterConstruction::BuildTank(): INFO: Making tank no. 3..."<<endl;
	
	new G4Tank(0,G4ThreeVector(Tank3PosX, Tank3PosY,Tank3PosZ), expHall_log,false,tankCopyNo,this);		
	tankCopyNo++;



	//## Make tank 4
	double Tank4PosX= Tank1PosX+fTankInfillGridSize/2;
	double Tank4PosY= Tank1PosY-1./3.*fTankInfillGridSize*sin(60.*TMath::Pi()/180.);
	double Tank4PosZ= Tank1PosZ;
	cout<<"G4MuonCounterConstruction::BuildTank(): INFO: Making tank no. 4..."<<endl;
	
	new G4Tank(0,G4ThreeVector(Tank4PosX, Tank4PosY,Tank4PosZ), expHall_log,false,tankCopyNo,this);		
	tankCopyNo++;


	//## Make tank 5
	double Tank5PosX= Tank1PosX-fTankInfillGridSize/2;
	double Tank5PosY= Tank1PosY-1./3.*fTankInfillGridSize*sin(60.*TMath::Pi()/180.);
	double Tank5PosZ= Tank1PosZ;
	cout<<"G4MuonCounterConstruction::BuildTank(): INFO: Making tank no. 5..."<<endl;
	
	new G4Tank(0,G4ThreeVector(Tank5PosX, Tank5PosY,Tank5PosZ), expHall_log,false,tankCopyNo,this);		
	tankCopyNo++;


	//## Create ground
	double fGroundThickness = 1*cm;
	double fGroundSizeX= expHall_x;
	double fGroundSizeY= expHall_y;	
	G4Box* ground_solid = new G4Box("ground_solid", fGroundSizeX/2, fGroundSizeY/2, fGroundThickness/2);

	G4LogicalVolume* ground_log = new G4LogicalVolume(ground_solid, MalargueSoil, "ground_solid"); 
  G4VisAttributes* ground_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  ground_log->SetVisAttributes(ground_va);
  G4VPhysicalVolume* ground_phys = new G4PVPlacement(0, G4ThreeVector(0, +fGroundSizeY/4,-fGroundThickness), "ground", ground_log, expHall_phys, false, 0); 	
	

	


}//close G4MuonCounterConstruction::BuildTank()


void G4MuonCounterConstruction::SetStripDesignMode(int value) { 
  fStripDesignMode= value; 
	updated= true;
}

void G4MuonCounterConstruction::SetStripReadoutMode(int value) { 
  fStripReadoutMode= value; 
	updated= true;
}


void G4MuonCounterConstruction::SetStripOptCouplingMode(int value) { 
  fStripOptCouplingMode= value; 
	updated= true;
}

void G4MuonCounterConstruction::SetFiberOptCouplingMode(int value) { 
  fFiberOptCouplingMode= value; 
	updated= true;
}

void G4MuonCounterConstruction::SetFiber(bool value) { 
  fUseFiber= value; 
	updated= true;
}


void G4MuonCounterConstruction::SetDoubleReadout(bool value) { 
  fUseDoubleReadout= value; 
	updated= true;
}


void G4MuonCounterConstruction::SetFiberRadius(double value) { 
	fFiberRadius= value; 
	updated= true;
}


void G4MuonCounterConstruction::SetFiberSizeZ(double value) {  
	fiber_z= value; 
	updated= true;
}


void G4MuonCounterConstruction::SetClad1Radius(double value) { 
	clad1_rmax= value; 
	updated= true;
}


void G4MuonCounterConstruction::SetClad2Radius(double value) { 
	clad2_rmax= value; 
	updated= true;
}
  

void G4MuonCounterConstruction::SetFiberAbsorptionLength(double value) { 
	fiber_absLength= value; 
	updated= true;
}
 
void G4MuonCounterConstruction::SetFiberDecayTime(double value) { 
	fFiberDecayTime= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetFiberCriticalAngle(double value) { 
	fFiberCriticalAngle= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetFiberCoreRefractiveIndex(double value) { 
	fFiberCoreRefractiveIndex= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetFiberClad1RefractiveIndex(double value) { 
	fFiberClad1RefractiveIndex= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetFiberClad2RefractiveIndex(double value) { 
	fFiberClad2RefractiveIndex= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetStripSize(G4ThreeVector vect) { 
	Xstrip= vect.x(); 
	Ystrip= vect.y();
	Zstrip= vect.z();		
  updated=true;
}


void G4MuonCounterConstruction::SetHousingThickness(double value) { 
	CoatingThickness= value; 
	updated= true;
}


void G4MuonCounterConstruction::SetHousingReflectivity(double value) { 
	CoatingReflectivity= value; 
	updated=true;
}


void G4MuonCounterConstruction::SetRealReflectivity(bool value) { 
	UseRealReflectivity= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetHousingSurfaceType(int value) { 
	CoatingSurfaceType= value; 
	updated= true;
}

void G4MuonCounterConstruction::SetGrooveSize(G4ThreeVector vect) { 
	Xgroove= vect.x(); 
  Ygroove= vect.y();
	Zgroove= vect.z();
	updated= true;
}
	

void G4MuonCounterConstruction::SetPMTSize(G4ThreeVector vect) { 
	pmt_length= vect.x(); 
  pmt_sizeY= vect.y();
	pmt_sizeZ= vect.z();
	updated=true;
}


void G4MuonCounterConstruction::SetPhotocathodeSize(G4ThreeVector vect) { 
	photocathode_length= vect.x(); 
  photocathode_sizeY= vect.y();
	photocathode_sizeZ= vect.z();
	updated= true;
}

void G4MuonCounterConstruction::SetPMTOpticalCouplingSize(G4ThreeVector vect) { 
	pmtOptCoupling_length= vect.x(); 
  pmtOptCoupling_sizeY= vect.y();
	pmtOptCoupling_sizeZ= vect.z();
	updated= true;
}


void G4MuonCounterConstruction::SetRealPhotocathode(bool value) { 
	UseRealPhotocathode= value; 
	updated=true;
}

void G4MuonCounterConstruction::UseYPlaneInSuperModule(bool value) { 
	HasYPlaneInSuperModule= value; 
	updated=true;
}


void G4MuonCounterConstruction::SetNumberOfStripsInPlane(int value) { 
	Nstrips= value; 
	updated=true;
}


void G4MuonCounterConstruction::SetTotalNumberOfStrips(int value) { 
	Nstrips_tot= value; 
	updated=true;
}


void G4MuonCounterConstruction::SetNumberOfPlanes(int value) { 
	Nplanes= value; 
	updated=true;
}


void G4MuonCounterConstruction::SetDistanceAmongPlanes(std::vector<double> vect) { 
	fPlaneDistance= vect; 
	updated=true;
}



void G4MuonCounterConstruction::SetDistanceAmongXYPlanes(double value) { 
	fXYPlaneDistance= value; 
	updated=true;
}


void G4MuonCounterConstruction::SetPlaneTiltAngle(double value) { 
	fPlaneTiltAngle= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetSuperPlaneModuleThickness(double value) { 
	fSuperPlaneModuleThickness= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetSuperPlaneModuleSizeZ(double value) { 
	fSuperPlaneModuleSizeZ= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetWorldMaterialName(G4String value) { 
	fWorldMaterialName= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetSuperPlaneMaterialName(G4String value) { 
	fSuperPlaneMaterialName= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetScintillationYield(double value) { 
	fScintillationYield= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetScintillationWeight(double value) { 
	fScintillationWeight= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetStripVetoTime(double value) { 
	fStripVetoTime= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetStripDecayTime(double value) { 
	fStripDecayTime= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetStripEnergyThreshold(double value) { 
	fStripEnergyThreshold= value; 
	updated=true;
}



void G4MuonCounterConstruction::SetPMTPhotoelectronYield(double value) { 
	fPMTPhotoelectronYield= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetPMTTransitTimeAverage(double value) { 
	fPMTTransitTimeAverage= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetPMTTransitTimeSpread(double value) { 
	fPMTTransitTimeSpread= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetSurfaceSeparationGap(double value) { 
	fSurfaceSeparationGap= value; 
	updated=true;
}

void G4MuonCounterConstruction::SetXMLParameters() {

	//##################################
	//## get some XML config
	//##################################
	cout<<"Getting parameters from XML file..."<<endl;

  Branch TopBranch = CentralConfig::GetInstance()->GetTopBranch("G4MuonCounterSimulatorUSC");
	
  TopBranch.GetChild("worldMaterial").GetData(fWorldMaterialType);
	if(fWorldMaterialType==1) fWorldMaterialName= "Vacuum";
	else if(fWorldMaterialType==2) fWorldMaterialName= "Air";
	else if(fWorldMaterialType==3) fWorldMaterialName= "MalargueSoil";
	else{
		cerr<<"Invalid world material type given in XML config...exit"<<endl;
	}

	fPlaneDistance.clear();
	fPlaneDistance.resize(0);

	//## muon counter design
  Branch MuonCounterModuleB = TopBranch.GetChild("MuonCounterModule");
  MuonCounterModuleB.GetChild("nStrip").GetData(Nstrips);
	MuonCounterModuleB.GetChild("nPlane").GetData(Nplanes);
	MuonCounterModuleB.GetChild("superplaneThickness").GetData(fSuperPlaneModuleThickness);
	MuonCounterModuleB.GetChild("superplaneHeight").GetData(fSuperPlaneModuleSizeZ);
	MuonCounterModuleB.GetChild("superplaneCasingMaterial").GetData(fSuperPlaneMaterialType);
	MuonCounterModuleB.GetChild("superplaneDistance").GetData(fPlaneDistance);
  MuonCounterModuleB.GetChild("enableYplane").GetData(HasYPlaneInSuperModule);
  MuonCounterModuleB.GetChild("planeXYDistance").GetData(fXYPlaneDistance);
  MuonCounterModuleB.GetChild("planeXYTiltAngle").GetData(fPlaneTiltAngle);
  MuonCounterModuleB.GetChild("surfaceSeparationGap").GetData(fSurfaceSeparationGap);
  
	//convert to G4 units
	fSuperPlaneModuleThickness*= CLHEP::mm / utl::mm;
	fSuperPlaneModuleSizeZ*= CLHEP::cm / utl::cm;
	fXYPlaneDistance*= CLHEP::cm / utl::cm;
	fPlaneTiltAngle*= CLHEP::degree / utl::degree;
	for(unsigned int i=0;i<fPlaneDistance.size();i++) fPlaneDistance[i]*= CLHEP::m / utl::m;
	fSurfaceSeparationGap*= CLHEP::nm/ utl::nanometer;

	if(fSuperPlaneMaterialType==1) fSuperPlaneMaterialName= "Vacuum";
	else if(fSuperPlaneMaterialType==2) fSuperPlaneMaterialName= "Air";
	else if(fSuperPlaneMaterialType==3) fSuperPlaneMaterialName= "PVC";
	else if(fSuperPlaneMaterialType==4) fSuperPlaneMaterialName= "GlassFiber";
	else{
		cerr<<"Invalid superplane material type given in XML config...exit"<<endl;
	}
	
	//print parsed values
	cout<<"## PARSED STRIP MODULE XML CONFIG ##"<<endl;
	cout<<"nStrip "<<Nstrips<<endl;
	cout<<"nPlane "<<Nplanes<<endl;
	cout<<"SuperPlaneMaterial "<<fSuperPlaneMaterialName<<endl;
	cout<<"SuperPlaneModuleThickness [mm]="<<fSuperPlaneModuleThickness/CLHEP::mm<<endl;
	cout<<"SuperPlaneModuleSizeZ [cm]="<<fSuperPlaneModuleSizeZ/CLHEP::cm<<endl;
	cout<<"PlaneDistance size"<<fPlaneDistance.size()<<endl;
	for(unsigned int i=0;i<fPlaneDistance.size();i++) {
		cout<< "PlaneDistance["<<i<<"] [m]="<<fPlaneDistance[i]/CLHEP::m <<endl;
	}
	cout<<"HasYPlaneInSuperModule? "<<HasYPlaneInSuperModule<<endl;
	cout<<"fXYPlaneDistance [cm]="<<fXYPlaneDistance/CLHEP::cm<<endl;
	cout<<"fPlaneTiltAngle [degree]="<<fPlaneTiltAngle/CLHEP::degree<<endl;
	cout<<"fSurfaceSeparationGap [ns]="<<fSurfaceSeparationGap/CLHEP::nm<<endl;

	//## scintillator strip design
	Branch ScintillatorStripB = TopBranch.GetChild("ScintillatorStrip");
	ScintillatorStripB.GetChild("stripSizeX").GetData(Xstrip);
	ScintillatorStripB.GetChild("stripSizeY").GetData(Ystrip);
	ScintillatorStripB.GetChild("stripSizeZ").GetData(Zstrip);
	ScintillatorStripB.GetChild("stripGrooveX").GetData(Xgroove);
	ScintillatorStripB.GetChild("stripGrooveY").GetData(Ygroove);
	ScintillatorStripB.GetChild("stripGrooveZ").GetData(Zgroove);
	ScintillatorStripB.GetChild("stripHousingThickness").GetData(CoatingThickness);
	ScintillatorStripB.GetChild("stripHousingUseSpectrumReflectivity").GetData(UseRealReflectivity);
	ScintillatorStripB.GetChild("stripHousingReflectivity").GetData(CoatingReflectivity);
	ScintillatorStripB.GetChild("stripHousingSurfaceType").GetData(CoatingSurfaceType);
	ScintillatorStripB.GetChild("scintillationYield").GetData(fScintillationYield);
	ScintillatorStripB.GetChild("scintillationWeight").GetData(fScintillationWeight);		
	ScintillatorStripB.GetChild("stripReadoutMode").GetData(fStripReadoutMode);
	ScintillatorStripB.GetChild("stripDesignMode").GetData(fStripDesignMode);
	ScintillatorStripB.GetChild("stripOpticalCouplingMode").GetData(fStripOptCouplingMode);	
	ScintillatorStripB.GetChild("stripEnergyThreshold").GetData(fStripEnergyThreshold);
	ScintillatorStripB.GetChild("stripVetoTime").GetData(fStripVetoTime);
	ScintillatorStripB.GetChild("stripDecayTime").GetData(fStripDecayTime);
	ScintillatorStripB.GetChild("stripEnergyThreshold").GetData(fStripEnergyThreshold);
	

	//convert to G4 units
	Xstrip*= CLHEP::cm / utl::cm;
	Ystrip*= CLHEP::cm / utl::cm;
	Zstrip*= CLHEP::cm / utl::cm;
	Xgroove*= CLHEP::cm / utl::cm;
	Ygroove*= CLHEP::cm / utl::cm;
	Zgroove*= CLHEP::cm / utl::cm;
	CoatingThickness*= CLHEP::mm / utl::mm;
	fStripVetoTime*= CLHEP::ns/ utl::ns;
	fStripDecayTime*= CLHEP::ns/ utl::ns;
	fStripEnergyThreshold*= CLHEP::MeV/ utl::MeV;

	//print parsed values
	cout<<"## PARSED STRIP XML CONFIG ##"<<endl;
	cout<<"StripSize XYZ [cm] =("<<Xstrip/CLHEP::cm<<","<<Ystrip/CLHEP::cm<<","<<Zstrip/CLHEP::cm<<")"<<endl;
	cout<<"GrooveSize XYZ [cm] =("<<Xgroove/CLHEP::cm<<","<<Ygroove/CLHEP::cm<<","<<Zgroove/CLHEP::cm<<")"<<endl;
	cout<<"CoatingThickness [mm] ="<<CoatingThickness/CLHEP::mm<<endl;
	cout<<"RealReflectivity? "<< UseRealReflectivity<<endl;
	cout<<"FixedReflectivity = "<< CoatingReflectivity<<endl;
	cout<<"CoatingSurfaceType = "<<CoatingSurfaceType<<endl;
	cout<<"StripReadoutMode = "<<fStripReadoutMode<<endl;
	cout<<"StripDesignMode = "<<fStripDesignMode<<endl;
	cout<<"StripOptCouplingMode = "<<fStripOptCouplingMode<<endl;
	cout<<"StripDecayTime [ns] = "<<fStripDecayTime/CLHEP::ns<<endl;
	cout<<"StripVetoTime [ns] = "<<fStripVetoTime/CLHEP::ns<<endl;
	cout<<"StripEnergyThreshold [MeV] = "<<fStripEnergyThreshold/CLHEP::MeV<<endl;
	

	//## WLS fiber design
	Branch WLSFiberB = TopBranch.GetChild("WLSFiber");
	WLSFiberB.GetChild("fiberLength").GetData(fiber_z);
	WLSFiberB.GetChild("fiberRadius").GetData(fFiberRadius);
	WLSFiberB.GetChild("fiberCaptureLength").GetData(fiber_captureLength);
	WLSFiberB.GetChild("fiberAttLength").GetData(fiber_absLength);
	WLSFiberB.GetChild("fiberDecayTime").GetData(fFiberDecayTime);
	WLSFiberB.GetChild("fiberCoreRefractiveIndex").GetData(fFiberCoreRefractiveIndex);
	WLSFiberB.GetChild("fiberClad1RefractiveIndex").GetData(fFiberClad1RefractiveIndex);
	WLSFiberB.GetChild("fiberClad2RefractiveIndex").GetData(fFiberClad2RefractiveIndex);
	WLSFiberB.GetChild("fiberOpticalCouplingMode").GetData(fFiberOptCouplingMode);

	//convert to G4 units
	fiber_z*= CLHEP::cm / utl::cm;
	fFiberRadius*= CLHEP::mm / utl::mm;
	fiber_captureLength*= CLHEP::mm / utl::mm;
	fiber_absLength*= CLHEP::m / utl::m;
	fFiberDecayTime*= CLHEP::ns / utl::ns;

	//print parsed values
	cout<<"## PARSED FIBER XML CONFIG ##"<<endl;
	cout<<"FiberLenght [cm] = "<<fiber_z/CLHEP::cm<<endl;
	cout<<"FiberRadius [mm] = "<<fFiberRadius/CLHEP::mm<<endl;
	cout<<"FiberCaptureLength [mm] = "<<fiber_captureLength/CLHEP::mm<<endl;
	cout<<"FiberAbsLength [m] = "<<fiber_absLength/CLHEP::m<<endl;
	cout<<"FiberDecayTime [ns] = "<<fFiberDecayTime/CLHEP::ns<<endl;
	cout<<"FiberCoreRefractiveIndex = "<<fFiberCoreRefractiveIndex<<endl;
	cout<<"FiberClad1RefractiveIndex = "<<fFiberClad1RefractiveIndex<<endl;
	cout<<"FiberClad2RefractiveIndex = "<<fFiberClad2RefractiveIndex<<endl;
	cout<<"FiberOptCouplingMode = "<<fFiberOptCouplingMode<<endl;


	//## WLS fiber design
	Branch PMTB = TopBranch.GetChild("PMT");
	PMTB.GetChild("pmtSizeX").GetData(pmt_length);
	PMTB.GetChild("pmtSizeY").GetData(pmt_sizeY);
	PMTB.GetChild("pmtSizeZ").GetData(pmt_sizeZ);
	PMTB.GetChild("photocathodeSizeX").GetData(photocathode_length);
	PMTB.GetChild("photocathodeSizeY").GetData(photocathode_sizeY);
	PMTB.GetChild("photocathodeSizeZ").GetData(photocathode_sizeZ);
	PMTB.GetChild("optCouplingSizeX").GetData(pmtOptCoupling_length);
	PMTB.GetChild("optCouplingSizeY").GetData(pmtOptCoupling_sizeY);
	PMTB.GetChild("optCouplingSizeZ").GetData(pmtOptCoupling_sizeZ);
	PMTB.GetChild("useSpectrumQE").GetData(UseRealPhotocathode);
	PMTB.GetChild("peYield").GetData(fPMTPhotoelectronYield);
	PMTB.GetChild("averageTransitTime").GetData(fPMTTransitTimeAverage);
	PMTB.GetChild("spreadTransitTime").GetData(fPMTTransitTimeSpread);


	//convert to G4 units
	pmt_length*= CLHEP::cm / utl::cm;
	pmt_sizeY*= CLHEP::cm / utl::cm;
	pmt_sizeZ*= CLHEP::cm / utl::cm;
	photocathode_length*= CLHEP::mm / utl::mm;
	photocathode_sizeY*= CLHEP::mm / utl::mm;
	photocathode_sizeZ*= CLHEP::mm / utl::mm;
	pmtOptCoupling_length*= CLHEP::mm / utl::mm;
	pmtOptCoupling_sizeY*= CLHEP::mm / utl::mm;
	pmtOptCoupling_sizeZ*= CLHEP::mm / utl::mm;
	fPMTTransitTimeAverage*= CLHEP::ns / utl::ns;
	fPMTTransitTimeSpread*= CLHEP::ns / utl::ns;
	

	//print parsed values
	cout<<"## PARSED PMT XML CONFIG ##"<<endl;
	cout<<"PMTSize XYZ [cm] =("<<pmt_length/CLHEP::cm<<","<<pmt_sizeY/CLHEP::cm<<","<<pmt_sizeY/CLHEP::cm<<")"<<endl;
	cout<<"PMTPhotocathodeSize XYZ [mm] =("<<photocathode_length/CLHEP::mm<<","<<photocathode_sizeY/CLHEP::mm<<","<<photocathode_sizeY/CLHEP::mm<<")"<<endl;
	cout<<"PMTOptCouplingSize XYZ [mm] =("<<pmtOptCoupling_length/CLHEP::mm<<","<<pmtOptCoupling_sizeY/CLHEP::mm<<","<<pmtOptCoupling_sizeY/CLHEP::mm<<")"<<endl;
	cout<<"UseSpectrumQE? "<<UseRealPhotocathode<<endl;
	cout<<"PMTYield [pe/MeV] = "<<fPMTPhotoelectronYield<<endl;
	cout<<"PMTTransitTimeAverage [ns] = "<<fPMTTransitTimeAverage/CLHEP::ns<<endl;
	cout<<"PMTTransitTimeSpread [ns] = "<<fPMTTransitTimeSpread/CLHEP::ns<<endl;
	

}//close SetXML

void G4MuonCounterConstruction::UpdateGeometry(){

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();


  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();


  updated=false;
}


