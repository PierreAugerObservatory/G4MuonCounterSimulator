/**
* @file G4MuonCounterConstruction.hh
* @class G4MuonCounterConstruction
* @brief Build the detector geometry, material & material properties.
*
* The detector is made up of a given number of detector planes, each composed by a X- and Y- plane.
* Each plane hosts a given number of strip modules. A single strip module is a scintillator bar, with a groove hosting 
* a WLS fiber. At both strip edges, a PMT is placed for the readout.  
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterConstruction_h_
#define _G4MuonCounterSimulatorUSC_G4MuonCounterConstruction_h_ 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4ScintillatorSD.hh"
#include "G4PMTSD.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"

namespace G4MuonCounterSimulatorUSC {

class DetectorMessenger;

class G4MuonCounterConstruction : public G4VUserDetectorConstruction
{
  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    G4MuonCounterConstruction();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~G4MuonCounterConstruction();

		/**
		* \brief Strip Readout Mode: FullReadout ==> readout with both scintillator+fiber, DoubleReadout ==> readout at both strip edges
		*/
		enum stripReadoutMode { kDoubleFiberReadout= 4, kSingleFiberReadout= 3, kDoubleReadout= 2, kSingleReadout= 1};
		/**
		* \brief Strip Design Mode
		*/
		enum stripDesignMode { kScintPlusFiber= 3, kScintPlusGroove= 2, kScint= 1};
		/**
		* \brief Strip Optical Coupling Mode
		*/
		enum stripOptCouplingMode { kGlueCoupling= 3, kGreaseCoupling= 2, kAirCoupling= 1, kNoCoupling= 0};	
		


  public:

		/**
		* \brief Build the whole geometry: first define materials via DefineMaterials(), then call the detector building
		*/
    G4VPhysicalVolume* Construct();
    
		
		/**
		* \brief Update the geometry if this has been modified
		*/
		void UpdateGeometry();
		/**
		* \brief Tells if the geometry has been modified
		*/
  	bool GetUpdated(){return updated;}


		/**
		* \brief Build a surface station tank
		*/
		void BuildTank();

		/**
		* \brief Get the logical world
		*/
		G4LogicalVolume* GetWorldVolume(){return expHall_log;}

		//### WLS Fiber
		/**
 		* \brief Turns on/off the WLS fiber in the geometry
 		*/
		void SetFiber(bool value);
		
		/**
 		* \brief Set the radius of the WLS fiber in the geometry
 		*/
		void SetFiberRadius(double value);

		/**
 		* \brief Get the radius of the WLS fiber in the geometry
 		*/
		double GetFiberRadius(){return fiber_rmax;}

		/**
 		* \brief Set the length of the WLS fiber in the geometry
 		*/
		void SetFiberSizeZ(double value);

		/**
 		* \brief Get the length of the WLS fiber in the geometry
 		*/
		double GetFiberSizeZ(){return fiber_z;}

		/**
 		* \brief Set the radius of the WLS fiber clad 1 in the geometry
 		*/
		void SetClad1Radius(double value);

		/**
 		* \brief Get the radius of the WLS fiber clad 1 in the geometry
 		*/
		double GetClad1Radius(){return clad1_rmax;}

		/**
 		* \brief Set the radius of the WLS fiber clad 2 in the geometry
 		*/
		void SetClad2Radius(double value);

		/**
 		* \brief Get the radius of the WLS fiber clad 2 in the geometry
 		*/
		double GetClad2Radius(){return clad2_rmax;}
		/**
 		* \brief Set the absorption length in the absorption region of the WLS fiber in the geometry
 		*/
		void SetFiberAbsorptionLength(double value);

		/**
 		* \brief Get the absorption length in the absorption region of the WLS fiber in the geometry
 		*/
		double GetFiberAbsorptionLength(){return fiber_absLength;}
    /**
 		* \brief Set the fiber decay time in ns
 		*/
		void SetFiberDecayTime(double value);
		/**
 		* \brief Get the fiber decay time in ns
 		*/
    static double GetFiberDecayTime(){return fFiberDecayTime;}
		/**
 		* \brief Set the fiber critical angle in rad
 		*/
		void SetFiberCriticalAngle(double value);
		/**
 		* \brief Get the fiber critical angle in rad
 		*/
    static double GetFiberCriticalAngle(){return fFiberCriticalAngle;}
		/**
 		* \brief Set the fiber core refractive index
 		*/
		void SetFiberCoreRefractiveIndex(double value);
		/**
 		* \brief Get the fiber core refractive index
 		*/
    static double GetFiberCoreRefractiveIndex(){return fFiberCoreRefractiveIndex;}
		/**
 		* \brief Set the fiber clad 1 refractive index
 		*/
		void SetFiberClad1RefractiveIndex(double value);
		/**
 		* \brief Get the fiber clad 1 refractive index
 		*/
    static double GetFiberClad1RefractiveIndex(){return fFiberClad1RefractiveIndex;}
		/**
 		* \brief Set the fiber clad 2 refractive index
 		*/
		void SetFiberClad2RefractiveIndex(double value);
		/**
 		* \brief Get the fiber clad 2 refractive index
 		*/
    static double GetFiberClad2RefractiveIndex(){return fFiberClad2RefractiveIndex;}
	
		//### SCINTILLATOR
		/**
 		* \brief Set the size of the scintillator bar in the geometry
 		*/
		void SetStripSize(G4ThreeVector vect);
		/**
 		* \brief Get the size of the scintillator bar in the geometry
 		*/
		G4ThreeVector GetScintSize(){return G4ThreeVector(Xstrip,Ystrip,Zstrip);}	
		/**
 		* \brief Set the size of the groove in the scintillator bar in the geometry
 		*/
		void SetGrooveSize(G4ThreeVector vect);
		/**
 		* \brief Get the size of the groove in the scintillator bar in the geometry
 		*/
		G4ThreeVector GetGrooveSize(){return G4ThreeVector(Xgroove,Ygroove,Zgroove);}
		/**
 		* \brief Get the size of the scintillator coating in the geometry
 		*/
		G4ThreeVector GetScintCoatingSize(){return G4ThreeVector(Xcoat,Ycoat,Zcoat);}
		/**
 		* \brief Set the thickness of the scintillator coating in the geometry
 		*/	
		void SetHousingThickness(double value);
		/**
 		* \brief Get the thickness of the scintillator coating in the geometry
 		*/	
		double GetHousingThickness(){return CoatingThickness;}
		/**
 		* \brief Set the reflectivity (from 0 to 1) of the scintillator coating in the geometry
 		*/	
	  void SetHousingReflectivity(double value);
		/**
 		* \brief Get the reflectivity of the scintillator coating in the geometry
 		*/
		double GetHousingReflectivity(){return CoatingReflectivity;}
		/**
 		* \brief Set the scintillator absolute scintillation yield (in 1/MeV)
 		*/	
	  void SetScintillationYield(double value);
		/**
 		* \brief Get the scintillator absolute scintillation yield (in 1/MeV)
 		*/
		double GetScintillationYield(){return fScintillationYield;}
		/**
 		* \brief Set the scintillation weight (if yield is equal to the true yield then weight=1)
 		*/	
	  void SetScintillationWeight(double value);
		/**
 		* \brief Get the scintillation weight (if yield is equal to the true yield then weight=1)
 		*/
		double GetScintillationWeight(){return fScintillationWeight;}
	

		/**
 		* \brief Turns on/off the real reflectivity mode (realistic spectrum reflectivity) of the scintillator coating in the geometry
 		*/
		void SetRealReflectivity(bool value);
		/**
 		* \brief Set the reflectivity surface type (polished, ground, etc...) of the scintillator coating in the geometry
 		*/
    void SetHousingSurfaceType(int value);	
		/**
 		* \brief Get the reflectivity surface type (polished, ground, etc...) of the scintillator coating in the geometry
 		*/
		int GetHousingSurfaceType(){return CoatingSurfaceType;}
		
		/**
 		* \brief Set the strip readout mode in the geometry
 		*/	
	  void SetStripReadoutMode(int value);
		/**
 		* \brief Get the strip readout mode in the geometry
 		*/
		int GetStripReadoutMode(){return fStripReadoutMode;}
		/**
 		* \brief Set the strip design mode in the geometry
 		*/	
	  void SetStripDesignMode(int value);
		/**
 		* \brief Get the strip design mode in the geometry
 		*/
		int GetStripDesignMode(){return fStripDesignMode;}
		/**
 		* \brief Set the strip optical coupling mode in the geometry
 		*/	
	  void SetStripOptCouplingMode(int value);
		/**
 		* \brief Get the strip optical coupling mode in the geometry
 		*/
		int GetStripOptCouplingMode(){return fStripOptCouplingMode;}
		/**
 		* \brief Set the fiber optical coupling mode in the geometry
 		*/	
	  void SetFiberOptCouplingMode(int value);
		/**
 		* \brief Get the fiber optical coupling mode in the geometry
 		*/
		int GetFiberOptCouplingMode(){return fFiberOptCouplingMode;}
		

		//### PMT
		/**
 		* \brief Set the size of the PMT in the geometry
 		*/	
		void SetPMTSize(G4ThreeVector vect);
		/**
 		* \brief Get the size of the PMT in the geometry
 		*/	
		G4ThreeVector GetPMTSize(){return G4ThreeVector(pmt_length,pmt_sizeY,pmt_sizeZ);}	
		/**
 		* \brief Set the size of the PMT photocathode in the geometry
 		*/
		void SetPhotocathodeSize(G4ThreeVector vect);
		/**
 		* \brief Get the size of the PMT photocathode in the geometry
 		*/
		G4ThreeVector GetPhotocathodeSize(){return G4ThreeVector(photocathode_length,photocathode_sizeY,photocathode_sizeZ);}

		/**
 		* \brief Set the size of the detector-PMT optical coupling area in the geometry
 		*/	
		void SetPMTOpticalCouplingSize(G4ThreeVector vect);
		/**
 		* \brief Get the size of the detector-PMT optical coupling area in the geometry
 		*/	
		G4ThreeVector GetPMTOpticalCouplingSize(){return G4ThreeVector(pmtOptCoupling_length,pmtOptCoupling_sizeY,pmtOptCoupling_sizeZ);}	

		/**
 		* \brief Turns on/off the real QE mode (realistic QE spectrum) of the PMT in the geometry
 		*/
		void SetRealPhotocathode(bool value);
		/**
 		* \brief Turns on/off the double readout mode in the geometry (true: 2PMTs, false: 1PMT+reflective edge)
 		*/
    void SetDoubleReadout(bool value);//true ==> 2 PMTs

		/**
 		* \brief Set the world material name in the geometry
 		*/
		void SetWorldMaterialName(G4String value);
		/**
 		* \brief Get the world material name in the geometry
 		*/
		G4String GetWorldMaterialName(){return fWorldMaterialName;}
		
		

		//### Module
		/**
 		* \brief Set the number of X-Y detector planes in the geometry
 		*/
		void SetNumberOfPlanes(int value);
		/**
 		* \brief Get the number of X-Y detector planes in the geometry
 		*/
    int GetNumberOfPlanes(){return Nplanes;}
		/**
 		* \brief Set the number of scintillator strips placed in a plane in the geometry
 		*/
		void SetNumberOfStripsInPlane(int value);
		/**
 		* \brief Get the number of scintillator strips placed in a plane in the geometry
 		*/
		int GetNumberOfStripsInPlane(){return Nstrips;}
		/**
 		* \brief Set the total number of scintillator strips in the geometry
 		*/
		void SetTotalNumberOfStrips(int value);
		/**
 		* \brief Get the total number of scintillator strips in the geometry
 		*/
		int GetTotalNumberOfStrips(){return Nstrips_tot;}

		/**
 		* \brief Set the distance among the detector planes in the geometry
 		*/
		void SetDistanceAmongPlanes(std::vector<double> vect);
		/**
 		* \brief Get the distance among the detector planes in the geometry
 		*/
		std::vector<double> GetDistanceAmongPlanes(){return fPlaneDistance;}
		/**
 		* \brief Insert the distance among the detector planes in the geometry
 		*/
		void InsertDistanceAmongPlanes(double value){fPlaneDistance.push_back(value);}

		/**
 		* \brief Set the distance among the X- and Y- detector plane modules in the geometry
 		*/
		void SetDistanceAmongXYPlanes(double value);
		/**
 		* \brief Get the distance among the X- and Y- detector plane modules in the geometry
 		*/
		double GetDistanceAmongXYPlanes(){return fXYPlaneDistance;}
		/**
 		* \brief Set the tilt angle among the X- and Y- detector plane modules in the geometry
 		*/
		void SetPlaneTiltAngle(double value);
		/**
 		* \brief Get the tilt angle among the X- and Y- detector plane modules in the geometry
 		*/
		double GetPlaneTiltAngle(){return fPlaneTiltAngle;}
		
		//### Super Module
		/**
 		* \brief Set the thickness of the superplane casing
 		*/
		void SetSuperPlaneModuleThickness(double value);
		/**
 		* \brief Get the thickness of the superplane casing
 		*/
		double GetSuperPlaneModuleThickness(){return fSuperPlaneModuleThickness;}
		/**
 		* \brief Set the Z size of the superplane casing
 		*/
		void SetSuperPlaneModuleSizeZ(double value);
		/**
 		* \brief Get the Z size of the superplane casing
 		*/
		double GetSuperPlaneModuleSizeZ(){return fSuperPlaneModuleSizeZ;}
		/**
 		* \brief Set the superplane casing material name in the geometry
 		*/
		void SetSuperPlaneMaterialName(G4String value);
		/**
 		* \brief Get the superplane casing material name in the geometry
 		*/
		G4String GetSuperPlaneMaterialName(){return fSuperPlaneMaterialName;}
		/**
 		* \brief Add the y plane in the supermodule in the geometry
 		*/
		void UseYPlaneInSuperModule(bool value);
		/**
 		* \brief Has the supermodule an Y plane?
 		*/
		bool IsYPlaneInSuperModule(){return HasYPlaneInSuperModule;}


		/**
 		* \brief Tells if the WLS fiber is present in the geometry
 		*/
		bool IsFiber(){return fUseFiber;}
		/**
 		* \brief Tells if the coating reflectivity is set to the realistic spectrum mode in the geometry
 		*/
		bool IsRealReflectivity(){return UseRealReflectivity;}
		/**
 		* \brief Tells if the PMT QE is set to the realistic spectrum mode in the geometry
 		*/
		bool IsRealPhotocathode(){return UseRealPhotocathode;}	
		/**
 		* \brief Tells if the strip module has a double readout in the geometry
 		*/
		bool IsDoubleReadout(){return fUseDoubleReadout;}
		
		/**
 		* \brief Set the strip veto time in ns
 		*/
		void SetStripVetoTime(double value);
		/**
 		* \brief Get the strip veto time in ns
 		*/
    static double GetStripVetoTime(){return fStripVetoTime;}
		/**
 		* \brief Set the strip energy threshold in MeV
 		*/
		void SetStripEnergyThreshold(double value);
		/**
 		* \brief Get the strip energy threshold in MeV
 		*/
    static double GetStripEnergyThreshold(){return fStripEnergyThreshold;}			
		/**
 		* \brief Set the strip decay time in ns
 		*/
		void SetStripDecayTime(double value);
		/**
 		* \brief Get the strip decay time in ns
 		*/
    static double GetStripDecayTime(){return fStripDecayTime;}
		
		
		/**
 		* \brief Set the PMT photoelectron yield in pe/MeV
 		*/
		void SetPMTPhotoelectronYield(double value);
		/**
 		* \brief Get the PMT photoelectron yield in pe/MeV
 		*/
    static double GetPMTPhotoelectronYield(){return fPMTPhotoelectronYield;}

		

		/**
 		* \brief Set the PMT transit time average in ns
 		*/
		void SetPMTTransitTimeAverage(double value);
		/**
 		* \brief Get the PMT transit time average in ns
 		*/
    static double GetPMTTransitTimeAverage(){return fPMTTransitTimeAverage;}
		/**
 		* \brief Set the PMT transit time spread in ns
 		*/
		void SetPMTTransitTimeSpread(double value);
		/**
 		* \brief Get the PMT transit time spread in ns
 		*/
    static double GetPMTTransitTimeSpread(){return fPMTTransitTimeSpread;}

		/**
 		* \brief Set the surface separation gap 
 		*/
		void SetSurfaceSeparationGap(double value);
		/**
 		* \brief Get the surface separation gap 
 		*/
    double GetSurfaceSeparationGap(){return fSurfaceSeparationGap;}
	

		/**
 		* \brief Set the thickness of the superplane casing
 		*/
		void SetDetectorParametersFromXML(bool value){fUseDetectorParametersFromXML= value;}
		/**
 		* \brief Get the thickness of the superplane casing
 		*/
		double UseDetectorParametersFromXML(){return fUseDetectorParametersFromXML;}


  private:

		/**
		* \brief Define the materials (elements, materials, material property tables) used in the geometry
		*/
		void DefineMaterials();
		/**
		* \brief Build the detector geometry (re-called when geometry is modified, i.e. in the config file): it must returns the physical world
		*/
  	G4VPhysicalVolume* BuildDetectorGeometry();
		/**
		* \brief Set the detector parameters from XML config file
		*/
		void SetXMLParameters();
		/**
		* \brief Set the detector parameters 
		*/
		void SetParameters();


		DetectorMessenger* detMessenger;
 		bool updated;
		bool fUseDetectorParametersFromXML;

		//Elements
		G4Element* H;
    G4Element* B;
		G4Element* C;
		G4Element* N;
		G4Element* O;
		G4Element* Na;
		G4Element* Fe;
		G4Element* Mg;
		G4Element* Mn;
		G4Element* Ca;
		G4Element* K;
		G4Element* P;
		G4Element* Cl;
  	
		//Materials 
  	G4Material* LXe;
		G4Material* Air;
  	G4Material* Vacuum;
  	G4Material* Al;
		G4Material* Glass;
  	G4Material* Polystyrene;
    G4Material* Polystyrene_doped;
		G4Material* Polystyrene_TiO2mixed;
  	G4Material* PMMA;
    G4Material* PPO;
		G4Material* POPOP;
		G4Material* MalargueSoil;
  	G4Material* Pethylene;
  	G4Material* fPethylene;
		G4Material* PolyvinylToluene;
		G4Material* Silicone;
  	G4Material* TiO2;
  	G4Material* SiO2;
		G4Material* B2O3;
		G4Material* B2O2;
    G4Material* Na2O;
		G4Material* Al2O3;
		G4Material* Fe2O3;
		G4Material* MgO;
		G4Material* MnO;
		G4Material* CaO;
		G4Material* K2O;
		G4Material* P2O5;
		G4Material* Borosilicate;																											
    G4Material* Water;
    G4Material* PVC;
		G4Material* BisphenolAEpichlorohydrin;
		G4Material* NButylGlycidylEther;
		G4Material* EpoxyShell;
		G4Material* BisphenolADiglycidylEther;
		G4Material* ButanediolDiglycidylEther;
		G4Material* HexanediamineTrimethyl; 
		G4Material* EpoxyTek;

		G4Material* HDPE;
		G4Material* Pyrex;
		G4Material* Interface;
		G4Material* Lucite;

		//Material properties tables
		G4MaterialPropertiesTable* VacuumPT;	
		G4MaterialPropertiesTable* AirPT;	
		G4MaterialPropertiesTable* OpticalCouplingPT;
		G4MaterialPropertiesTable* ScintillatorPT;
		G4MaterialPropertiesTable* IdealScintCoatingPT;
		G4MaterialPropertiesTable* RealScintCoatingPT;
		
  	G4MaterialPropertiesTable* WLSFiberPT;
    G4MaterialPropertiesTable* WLSClad1PT;	 
	  G4MaterialPropertiesTable* WLSClad2PT;
		G4MaterialPropertiesTable* PMTPT;
		G4MaterialPropertiesTable* IdealPhotocathodePT;
		G4MaterialPropertiesTable* RealPhotocathodePT;



		//Volumes
		G4Box* expHall_box;
  	G4LogicalVolume* expHall_log;
  	G4VPhysicalVolume* expHall_phys;
		G4VisAttributes* expHall_va;

    double expHall_x;
    double expHall_y;
    double expHall_z;

		//scintillator
		double Xstrip;
		double Ystrip;
		double Zstrip;

		double Xcoat;
		double Ycoat;
		double Zcoat;

		double Xgroove;
		double Ygroove;
		double Zgroove;

		double fScintillationYield;
		double fScintillationWeight;
				
		double CoatingThickness;
    double CoatingReflectivity;
		int CoatingSurfaceType;
    bool UseRealReflectivity;

		double fSuperPlaneModuleThickness;		
		double fSuperPlaneModuleSizeZ;
		bool HasYPlaneInSuperModule;
		
	  //fiber variables
		double fFiberRadius;
		bool fUseFiber;
		
		double fiber_rmin;
		double fiber_rmax;
		double fiber_z;
		double fiber_sphi;
		double fiber_ephi;
		//clad 1 variables		
		double clad1_rmin;
		double clad1_rmax;
		double clad1_z;
		double clad1_sphi;
		double clad1_ephi;
		//clad 2 variables
		double clad2_rmin;
		double clad2_rmax;
		double clad2_z;
		double clad2_sphi;
		double clad2_ephi;
		
		double fiber_absLength;
		double fiber_captureLength;

		//PMT variables
		double pmt_sizeY;
    double pmt_sizeZ;
    double pmt_length;
		
    double photocathode_sizeY;
    double photocathode_sizeZ;
		double photocathode_length;

		double pmtOptCoupling_sizeY;
    double pmtOptCoupling_sizeZ;
    double pmtOptCoupling_length;
 
		//strip options
		int fStripDesignMode;
		int fStripReadoutMode;
		int fStripOptCouplingMode;
		int fFiberOptCouplingMode;

    bool fUseDoubleReadout;
		bool UseRealPhotocathode;
		
		int Nstrips;
		int Nplanes;
		int Nstrips_tot;
		std::vector<double> fPlaneDistance;
		double fXYPlaneDistance;
		double fPlaneTiltAngle;

		double fSurfaceSeparationGap;

		//World material name
		G4String fWorldMaterialName;
		int fWorldMaterialType;

		//Superplane material name
		G4String fSuperPlaneMaterialName;
		int fSuperPlaneMaterialType;
		
		//Sensitive Detectors
  	static G4ScintillatorSD* Scint_SD;
  	static G4PMTSD* PMT_SD;
		
		
		static double fStripVetoTime;
		static double fStripEnergyThreshold;
		static double fStripDecayTime;
		static double fFiberDecayTime;
		static double fFiberCriticalAngle;
		static double fFiberCoreRefractiveIndex;
		static double fFiberClad1RefractiveIndex;
		static double fFiberClad2RefractiveIndex;
		static double fPMTTransitTimeAverage;
		static double fPMTTransitTimeSpread;
		static double fPMTPhotoelectronYield;

		

		//friend class TStationSimData;

};

}//close namespace

#endif /*G4MuonCounterConstruction_h*/
