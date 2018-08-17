/**
* @file G4PMTSD.hh
* @class G4PMTSD
* @brief Define the PMT sensitive detector, create and handle the PMT hits
*
* When a photon hits the photocathode of a given PMT, a corresponding hit is created, or, if the hit already exists, the current information * is added to the previous hit. 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4PMTSD_h
#define _G4MuonCounterSimulatorUSC_G4PMTSD_h 1

#include "G4DataVector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4PMTHit.hh"

class G4Step;
class G4HCofThisEvent;

namespace G4MuonCounterSimulatorUSC {

class G4PMTSD : public G4VSensitiveDetector
{

public:

	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4PMTSD(G4String name);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~G4PMTSD();
  /** 
	\brief Initialize hit collections.
 	*/
  void Initialize(G4HCofThisEvent* HCE);
	/** 
	\brief Create scintillator hits
 	*/
  bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  
	/** 
	\brief Create scintillator hits, keeping the step constant
 	*/
  //A version of processHits that keeps aStep constant
  bool ProcessHits_constStep(const G4Step* aStep,G4TouchableHistory* ROhist);
	/** 
	\brief End of event operations
 	*/
  void EndOfEvent(G4HCofThisEvent* HCE);
  /** 
	\brief Clear operations
 	*/
  void clear();
	/** 
	\brief Drawing of scintillator hits
 	*/
  void DrawAll();
	/** 
	\brief Printing of scintillator hits
 	*/
  void PrintAll();
  
	/** 
	\brief Initialize the arrays to store PMT positions
 	*/
  //Initialize the arrays to store pmt positions
  inline void InitPMTs(int nPMTs){
    if(pmtPositionsX)delete pmtPositionsX;
    if(pmtPositionsY)delete pmtPositionsY;
    if(pmtPositionsZ)delete pmtPositionsZ;
    pmtPositionsX=new G4DataVector(nPMTs);
    pmtPositionsY=new G4DataVector(nPMTs);
    pmtPositionsZ=new G4DataVector(nPMTs);
  }
	/** 
	\brief Store a PMT position
 	*/
  //Store a pmt position
  inline void SetPMTPos(int n,double x,double y,double z){
    if(pmtPositionsX)pmtPositionsX->insertAt(n,x);
    if(pmtPositionsY)pmtPositionsY->insertAt(n,y);
    if(pmtPositionsZ)pmtPositionsZ->insertAt(n,z);
  }
  
private:
	/** 
	\brief Collection of PMT hits
 	*/
  PMTHitsCollection* pmtHitCollection;
  
	/** 
	\brief PMT X position
 	*/
  G4DataVector* pmtPositionsX;
	/** 
	\brief PMT Y position
 	*/
  G4DataVector* pmtPositionsY;
	/** 
	\brief PMT Z position
 	*/
  G4DataVector* pmtPositionsZ;
};

}//close namespace

#endif

