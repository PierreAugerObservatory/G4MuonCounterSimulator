/**
* @file G4MuonCounterRunScorer.hh
* @class G4MuonCounterRunScorer
* @brief Define a user G4Run
* Implementation for multi-functional-detector and primitive scorer.
* This G4MuonCounterRunScorer class has collections which accumulate an event information into a run information.
*
* @author S. Riggi
* @date 30/01/2011
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterRunScorer_h_
#define _G4MuonCounterSimulatorUSC_G4MuonCounterRunScorer_h_ 1

#include "G4Run.hh"
#include "G4Event.hh"

#include "G4THitsMap.hh"
#include <vector>

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterRunScorer : public G4Run {

public:
  // constructor and destructor.
  //  vector of multifunctionaldetector name has to given to constructor.
  G4MuonCounterRunScorer(const std::vector<G4String> mfdName);
  virtual ~G4MuonCounterRunScorer();

public:
  // virtual method from G4Run. 
  // The method is overriden in this class for scoring.
  virtual void RecordEvent(const G4Event*);

  // Access methods for scoring information.
  // - Number of HitsMap for this RUN. 
  //   This is equal to number of collections.
  int GetNumberOfHitsMap() const {return theRunMap.size();}
  // - Get HitsMap of this RUN.
  //   by sequential number, by multifucntional name and collection name,
  //   and by collection name with full path.
  G4THitsMap<double>* GetHitsMap(int i){return theRunMap[i];}
  G4THitsMap<double>* GetHitsMap(const G4String& detName, const G4String& colName);
  G4THitsMap<double>* GetHitsMap(const G4String& fullName);
  // - Dump All HitsMap of this RUN.
  //   This method calls G4THisMap::PrintAll() for individual HitsMap.
  void DumpAllScorer();

private:
  std::vector<G4String> theCollName;
  std::vector<int> theCollID;
  std::vector<G4THitsMap<double>*> theRunMap;
};

}//close namespace

#endif
