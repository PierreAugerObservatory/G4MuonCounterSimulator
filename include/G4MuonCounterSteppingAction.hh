/**
* @file G4MuonCounterSteppingAction.hh
* @class G4MuonCounterSteppingAction
* @brief Handle all operations to be performed at step level
*
* Allows to access to step information, track information, ...
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterSteppingAction_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterSteppingAction_h 1

#include "UserTrackInformation.hh"
#include "UserEventInformation.hh"

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "G4OpBoundaryProcess.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterRecorderBase;
class G4MuonCounterEventAction;
class G4MuonCounterTrackingAction;
class SteppingMessenger;

class G4MuonCounterSteppingAction : public G4UserSteppingAction
{
	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
  	G4MuonCounterSteppingAction(G4MuonCounterRecorderBase*);
		/**
		* \brief Class destructor: free allocated memory
		*/
  	~G4MuonCounterSteppingAction();

		/**
		* \brief Absorption flag
		*/
		enum absorptionFlag{ kAbsInBoundary= 1, kAbsInScint= 2, kAbsInCoat= 3, kAbsInGroove= 4, kAbsInFiber= 5 };



		/**
		* \brief Handle operation to be performed for the current step
		*/
  	virtual void UserSteppingAction(const G4Step*);
		/**
		* \brief Enable/Disable primaries to go one step in the scintillator volume before being killed
		*/
  	void SetOneStepPrimaries(bool b){oneStepPrimaries=b;}
		/**
		* \brief Get the one step flag
		*/
  	bool GetOneStepPrimaries(){return oneStepPrimaries;}
  

	private:
		bool CountTracksOnSurface(const G4Step * theStep,G4String physVolName);	
		bool HandleOpticalPhotonTrack(const G4Step * theStep);
		bool HandlePrimaryTrack(const G4Step * theStep);


	private:
  	G4MuonCounterRecorderBase* recorder;
  	bool oneStepPrimaries;
  	SteppingMessenger* steppingMessenger;
		static G4OpBoundaryProcess* boundary;
		G4OpBoundaryProcessStatus boundaryStatus;
		bool recordAbsorptionInfo;

		UserTrackInformation* trackInformation;
		UserEventInformation* eventInformation;
};

}//close namespace

#endif
