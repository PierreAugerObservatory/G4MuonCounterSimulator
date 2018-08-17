/**
* @file G4MuonCounterRunAction.hh
* @class G4MuonCounterRunAction
* @brief Handle all operations to be performed at run level
*
* Manage operations before the run starts (i.e. creation of output data structures, ...), actions to be performed for each event (access and storing of simulation information), operations after the run ends (i.e. saving of simulation outputs, ...).
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterRunAction_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterRunAction_h 1

#include "G4UserRunAction.hh"
#include "TScintHit.hh"
#include "TPMTHit.hh"
#include "TEventSimData.hh"
#include "TStationSimData.hh"
#include "TParticleSimData.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TVector3.h"

#include "vector"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterRecorderBase;
class RunMessenger;
class G4MuonCounterConstruction;
class G4MuonCounterPrimaryGenerator;

class G4MuonCounterSimulator;

class G4MuonCounterRunAction : public G4UserRunAction
{
public:
	/** 
	\brief Class constructor: initialize structures.
 	*/
  G4MuonCounterRunAction(G4MuonCounterRecorderBase*,G4MuonCounterConstruction*,G4MuonCounterPrimaryGenerator*);
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~G4MuonCounterRunAction();
	/**
	* \brief Generate new RUN object for MultiFunctionalDetector scheme
	*/
	virtual G4Run* GenerateRun();
  /**
	* \brief Operations to be done before run starts
	*/
  void BeginOfRunAction(const G4Run*);
	/**
	* \brief Operations to be done when run ends
	*/
  void EndOfRunAction(const G4Run*);
	/**
	* \brief Fill simulation information
	*/
  void FillEventInfo();
	/**
	* \brief Init data structures
	*/
  void Init();
	/**
	* \brief Manage run primitive scorers
	*/
	void HandlePrimitiveScorers(const G4Run*);

	/**
	* \brief Turn on/off the storage of huge vectors for scintillator and PMT hits
	*/
	void StoreFullInfo(bool choice){fStoreFullInfo=choice;}



private:
	/**
	* \brief Base recorder for this class
	*/
  G4MuonCounterRecorderBase* recorder;
	/**
	* \brief Run messenger for this class
	*/
	RunMessenger* runMessenger;
	/**
	* \brief Detector costruction for this class
	*/
	G4MuonCounterConstruction* Detector;
  /**
	* \brief Detector costruction for this class
	*/
	G4MuonCounterPrimaryGenerator* generator;


	/**
	* \brief Turn on/off the storage of huge vectors for scintillator and PMT hits
	*/
	bool fStoreFullInfo;


	//hit collections	
	/**
	* \brief Scintillator hit collection flag
	*/
	int scintCollID;
	/**
	* \brief PMT hit collection flag
	*/
  int pmtCollID;

	/**
	* \brief Collection of ROOT scintillator hits 
	*/
	std::vector<TScintHit> fScintHitCollection;
	std::vector<TScintHit> fSelectedScintHitCollection;
  
  
	/**
	* \brief Collection of ROOT PMT hits 
	*/
	std::vector<TPMTHit> fPMTHitCollection;


	int nStrips;
	int nPlanes;
	double StripSeparation;

	/**
	* \brief Vector of MultiFunctionalDetector names
	*/
  std::vector<G4String> theSDName; 
	/**
	* \brief nMuons hitting each detector plane 
	*/
	std::vector<double> fMuonInNumber;
	std::vector<double> fMuonOutNumber;
	std::vector<double> fMuonInNumberInEvent;
	std::vector<double> fMuonOutNumberInEvent;
	/**
	* \brief nEm hitting each detector plane 
	*/
	std::vector<double> fEmInNumber;
	std::vector<double> fEmOutNumber;
	std::vector<double> fEmInNumberInEvent;
	std::vector<double> fEmOutNumberInEvent;
	/**
	* \brief nHadrons hitting each detector plane 
	*/
	std::vector<double> fHadronInNumber;
	std::vector<double> fHadronOutNumber;
	std::vector<double> fHadronInNumberInEvent;
	std::vector<double> fHadronOutNumberInEvent;

	double fHitSigmaX;
	double fHitSigmaY;
	double fHitSigmaZ;

	/**
	* \brief Pointer to the G4MuonCounterSimulator
	*/
	G4MuonCounterSimulator* fG4MuonCounterSimulator;
	/**
	* \brief Pointer to current TEventSimData
	*/
	TEventSimData* currentEventSimData;	
	/**
	* \brief Pointer to current TStationSimData
	*/
	TStationSimData* fStationSimData;
	
	


};

}//close namespace

#endif
