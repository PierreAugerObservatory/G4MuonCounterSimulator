/**
* @file G4MuonCounterRunScorer.hh
* @class G4MuonCounterRunScorer
* @brief Define a user G4Run
* G4MuonCounterRunScorer class is for accumulating scored quantities which is scored 
* using G4MultiFunctionalDetector and G4VPrimitiveScorer.
* Accumulation is done using G4THitsMap object.
*
* The constructor G4MuonCounterRunScorer(const std::vector<G4String> mfdName)
* needs a vector filled with MultiFunctionalDetector names which
* was assigned at instantiation of MultiFunctionalDetector(MFD).
* Then G4MuonCounterRunScorer constructor automatically scans primitive scorers
* in the MFD, and obtains collectionIDs of all collections associated
* to those primitive scorers. Futhermore, the G4THitsMap objects 
* for accumulating during a RUN are automatically created too.
* (*) Collection Name is same as primitive scorer name.
* 
* The resultant information is kept inside G4MuonCounterRunScorer objects as data members.
* std::vector<G4String> theCollName;            // Collection Name,
* std::vector<G4int> theCollID;                 // Collection ID,
* std::vector<G4THitsMap<G4double>*> theRunMap; // HitsMap for RUN.
*
* The resultant HitsMap objects are obtain using access method,
* GetHitsMap(..).
*
* @author S. Riggi
* @date 30/01/2011
*/


#include "G4MuonCounterRunScorer.hh"
#include "G4SDManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include <G4MuonCounterSimulator.hh>

using namespace std;
using namespace G4MuonCounterSimulatorUSC;


//
//  Constructor. 
//   (The vector of MultiFunctionalDetector name has to given.)
G4MuonCounterRunScorer::G4MuonCounterRunScorer(const std::vector<G4String> mfdName): G4Run(){

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //=================================================
  //  Initialize RunMaps for accumulation.
  //  Get CollectionIDs for HitCollections.
  //=================================================

	//## Looping over all MFD
  int Nmfd = mfdName.size();
  for(int idet=0;idet<Nmfd;idet++){ 
    G4String detName = mfdName[idet];

    // Seek and Obtain MFD objects from SDmanager.
    G4MultiFunctionalDetector* mfd = (G4MultiFunctionalDetector*)(SDman->FindSensitiveDetector(detName));
    
    if(mfd){
			//## Loop over the registered primitive scorers
			for(int icol=0;icol<mfd->GetNumberOfPrimitives();icol++){
	    
				// Get Primitive Scorer object.
	    	G4VPrimitiveScorer* scorer=mfd->GetPrimitive(icol);
	    	// collection name and collectionID for HitsCollection,
        // where type of HitsCollection is G4THitsMap in case of primitive scorer.
        // The collection name is given by <MFD name>/<Primitive Scorer name>.
	    	G4String collectionName = scorer->GetName();
	   	 	G4String fullCollectionName = detName+"/"+collectionName;
	    	int collectionID = SDman->GetCollectionID(fullCollectionName);
	    
	    	if(collectionID>= 0){
					//cout << "++ "<<fullCollectionName<< " id " << collectionID << endl;
					// Store obtained HitsCollection information into data members.
					// And, creates new G4THitsMap for accumulating quantities during RUN.
					theCollName.push_back(fullCollectionName);
					theCollID.push_back(collectionID);
					theRunMap.push_back(new G4THitsMap<G4double>(detName,collectionName));
	    	}
				else{
					cout << "** collection " << fullCollectionName << " not found. "<<G4endl;
	    	}
			}//end loop scorers
    }//if mfd
  }//end loop MFD

}//close constructor

//
// Destructor
//    clear all data members.
G4MuonCounterRunScorer::~G4MuonCounterRunScorer()
{
  //--- Clear HitsMap for RUN
  int Nmap = theRunMap.size();
  for(int i= 0;i<Nmap;i++){
    if(theRunMap[i]) theRunMap[i]->clear();
  }
  theCollName.clear();
  theCollID.clear();
  theRunMap.clear();
}//close destructor

//
//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a Run.
void G4MuonCounterRunScorer::RecordEvent(const G4Event* aEvent){

  numberOfEvent++;  // This is an original line.

  //=============================
  // HitsCollection of This Event
  //============================
  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  if (!HCE) return;

  //=======================================================
  // Sum up HitsMap of this Event  into HitsMap of this RUN
  //=======================================================
  int Ncol= theCollID.size();
  for(int i=0;i<Ncol;i++){  // Loop over HitsCollection
    G4THitsMap<double>* EvtMap=0;
    if(theCollID[i] >= 0 ){           // Collection is attached to HCE
      EvtMap = (G4THitsMap<double>*)(HCE->GetHC(theCollID[i]));
    }
		else{
      cout <<" Error EvtMap Not Found "<< i << endl;
    }
    if(EvtMap)  {
      //=== Sum up HitsMap of this event to HitsMap of RUN.===
      *theRunMap[i] += *EvtMap;
      //======================================================
    }
  }//end loop HitsCollection

}//close RecordEvent()


//=================================================================
//  Access method for HitsMap of the RUN
//
//-----
// Access HitsMap.
//  By  MultiFunctionalDetector name and Collection Name.
G4THitsMap<double>* G4MuonCounterRunScorer::GetHitsMap(const G4String& detName,const G4String& colName){
    
	G4String fullName = detName+"/"+colName;
  return GetHitsMap(fullName);

}//close GetHitsMap()

//-----
// Access HitsMap.
//  By full description of collection name, that is
//    <MultiFunctional Detector Name>/<Primitive Scorer Name>
G4THitsMap<G4double>* G4MuonCounterRunScorer::GetHitsMap(const G4String& fullName){
    
	int Ncol = theCollName.size();
  for(int i=0;i<Ncol;i++){
		if (theCollName[i] == fullName ){
	    return theRunMap[i];
		}
  }
  return NULL;
}//close GetHitsMap()

//-----
// - Dump All HitsMap of this RUN. (for debuging and monitoring of quantity).
//   This method calls G4THisMap::PrintAll() for individual HitsMap.
void G4MuonCounterRunScorer::DumpAllScorer(){

  // - Number of HitsMap in this RUN.
  int n = GetNumberOfHitsMap();
  //GetHitsMap and dump values.
  for(int i=0;i<n;i++){
    G4THitsMap<double>* RunMap =GetHitsMap(i);
    if (RunMap) {
      cout << " PrimitiveScorer RUN " 
	     	   << RunMap->GetSDname() <<","<< RunMap->GetName() << endl;
      cout << " Number of entries " << RunMap->entries() << endl;
      std::map<int,double*>::iterator itr = RunMap->GetMap()->begin();
      for(; itr != RunMap->GetMap()->end(); itr++) {
				cout << "  copy no.: " << itr->first
	       		 << "  Run Value : " << *(itr->second) 
	       		 << endl;
      }
    }
  }

}//close DumpAllScorer()


