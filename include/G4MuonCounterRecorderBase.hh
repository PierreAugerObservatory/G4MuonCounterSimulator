
#ifndef _G4MuonCounterSimulatorUSC_RECORDER_BASE_H_
#define _G4MuonCounterSimulatorUSC_RECORDER_BASE_H_

// The following objects are the arguments to the methods
// invoked in the user action classes.  In other words, they
// contain the variables that we are normally able to record
// in Geant.

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4Step.hh"

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterRecorderBase {

public:

  virtual ~G4MuonCounterRecorderBase() {};

  // The following a list of methods that correspond to the available 
  // user action classes in Geant 4.0.1.   In this base class, the
  // methods are defined to do nothing.

  virtual void RecordBeginOfRun(const G4Run*) = 0;
  virtual void RecordEndOfRun(const G4Run*) = 0;
  virtual void RecordBeginOfEvent(const G4Event*) {};
  virtual void RecordEndOfEvent(const G4Event*) {};
  virtual void RecordTrack(const G4Track*) {};
  virtual void RecordStep(const G4Step*) {};

};

}//close namespace

#endif
