/**
* @file G4MuonCounterSteppingVerbose.hh
* @class G4MuonCounterSteppingVerbose
* @brief Handle the verbosity of the stepping action
* @author S. Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterSteppingVerbose_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterSteppingVerbose_h 1


#include "G4SteppingVerbose.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterSteppingVerbose : public G4SteppingVerbose
{
 public:   
	/** 
	\brief Class constructor
 	*/
   G4MuonCounterSteppingVerbose();
	/**
	* \brief Class destructor
	*/
  ~G4MuonCounterSteppingVerbose();
	/**
	* \brief Set which information has to be printed for the current step
	*/
   void StepInfo();
	/**
	* \brief Set which information has to be printed for the track in the current step
	*/
   void TrackingStarted();

};

}//close namespace

#endif
