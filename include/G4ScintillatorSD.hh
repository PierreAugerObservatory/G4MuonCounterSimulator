/**
* @file G4ScintillatorSD.hh
* @class G4ScintillatorSD
* @brief Define the scintillator sensitive detector, create and handle the scintillator hits
*
* When a particle hits the scintillator, a corresponding hit is created only when the energy deposit is >0. If the hit already exists, the current information is added to the previous hit. 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4ScintillatorSD_h
#define _G4MuonCounterSimulatorUSC_G4ScintillatorSD_h 1

#include "G4ScintillatorHit.hh"

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

namespace G4MuonCounterSimulatorUSC {

class G4ScintillatorSD : public G4VSensitiveDetector{

public:

	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4ScintillatorSD(G4String name);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~G4ScintillatorSD();
  
	/**
 	* \brief Set the scintillator energy threshold in MeV
 	*/
	void SetScintillatorEnergyThreshold(double value) {fScintillatorEnergyThreshold= value;}

	/** 
	\brief Initialize hit collections.
 	*/
  void Initialize(G4HCofThisEvent* HCE);
	/** 
	\brief Create scintillator hits
 	*/
  bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
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
  
private:
	/** 
	\brief Collection of scintillator hits
 	*/
  ScintHitsCollection* scintHitCollection;
  
	double fScintillatorEnergyThreshold;
	

};

}//close namespace

#endif

