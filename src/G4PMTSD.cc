/**
* @file G4PMTSD.cc
* @class G4PMTSD
* @brief Define the PMT sensitive detector, create and handle the PMT hits
*
* When a photon hits the photocathode of a given PMT, a corresponding hit is created, or, if the hit already exists, the current information * is added to the previous hit. 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4PMTSD.hh"
#include "G4PMTHit.hh"
#include "G4MuonCounterConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

#include "G4VProcess.hh"

using namespace G4MuonCounterSimulatorUSC;

G4PMTSD::G4PMTSD(G4String name)
  :G4VSensitiveDetector(name),pmtHitCollection(0),pmtPositionsX(0)
  ,pmtPositionsY(0),pmtPositionsZ(0)
{
  collectionName.insert("pmtHitCollection");
}

G4PMTSD::~G4PMTSD()
{}


void G4PMTSD::Initialize(G4HCofThisEvent* HCE){

  pmtHitCollection = new PMTHitsCollection(SensitiveDetectorName,collectionName[0]); 
  //Store collection with event and keep ID
  static int HCID = -1;
  if(HCID<0){ 
    HCID = GetCollectionID(0); 
  }
  HCE->AddHitsCollection( HCID, pmtHitCollection );
}


bool G4PMTSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist){
  return ProcessHits_constStep(aStep,ROhist);
}

//Generates a hit and uses the postStepPoint's mother volume replica number
//PostStepPoint because the hit is generated manually when the photon is
//absorbed by the photocathode
bool G4PMTSD::ProcessHits_constStep(const G4Step* aStep,G4TouchableHistory* ){

  //need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
 
  //#####################################################################
  //##  GetCopyNumber(int depth) o GetReplicaNumber(int depth) dovrebbero
  //##  essere la stessa cosa.
  //##     - depth=0 resituisce il numero di copia del volume corrente
  //##     - depth=1 restituisce il numero di copia del volume madre (il primo volume madre che contiene il volume corrente)
  //##   In questo caso pmtNumber è il numero di PMT in cui si trova il fotocatodo corrente
  //##    vedi FORUM http://hypernews.slac.stanford.edu:5090/HyperNews/geant4/get/geometry/69/1.html  
  //User replica number 1 since photocathode is a daughter volume
  //to the pmt which was replicated
  
	//get current PMT id - depth=1
  int pmtNumber= aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
  G4VPhysicalVolume* physVol= aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1);	
  
	//get current STRIP id - depth=2
  int stripNumber= aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(2);
  G4VPhysicalVolume* ModulePhysVol= aStep->GetPostStepPoint()->GetTouchable()->GetVolume(2);	

	//get current PLANE id - depth=3
	int planeNumber= aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(3);
  G4VPhysicalVolume* PlaneModulePhysVol= aStep->GetPostStepPoint()->GetTouchable()->GetVolume(3);	
	
	//get current SUPERPLANE id - depth=4 without inside casing, depth=5 without inside casing, 
	int superplaneNumber= aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(5);
  G4VPhysicalVolume* SuperPlaneModulePhysVol= aStep->GetPostStepPoint()->GetTouchable()->GetVolume(5);	

	/*
	cout<<"******** PMT *****************"<<endl;
  cout<<"current PMTVolName="<<physVol->GetName()<<"   copyNo="<<pmtNumber<<endl;
  cout<<"current ModuleVolName="<< ModulePhysVol->GetName()<<"   copyNo="<<stripNumber<<endl;
	cout<<"current PlaneModuleVolName="<< PlaneModulePhysVol->GetName()<<"   copyNo="<<planeNumber<<endl;
	cout<<"current SuperPlaneModuleVolName="<< SuperPlaneModulePhysVol->GetName()<<"   copyNo="<<superplaneNumber<<endl;
  cout<<"******************************"<<endl;
	*/

	//cout<<"current Hit PMT Number= "<<pmtNumber<<"  "<<physVol->GetCopyNo()<<endl;

	G4ThreeVector currentPrePos = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector currentPostPos = aStep->GetPostStepPoint()->GetPosition();
  //cout<<"currentPos ==> "<<" POST X,Y,Z "<<currentPostPos.x()<<"  "<<currentPostPos.y()<<"  "<<currentPostPos.z()<<endl;
  //cout<<"currentPos ==> "<<" PRE X,Y,Z "<<currentPrePos.x()<<"  "<<currentPrePos.y()<<"  "<<currentPrePos.z()<<endl;
   

  //Find the correct hit collection
  int n= pmtHitCollection->entries();
  G4PMTHit* hit=NULL;
  for(int i=0;i<n;i++){
    if((*pmtHitCollection)[i]->GetPMTNumber()==pmtNumber && (*pmtHitCollection)[i]->GetStripNumber()==stripNumber && (*pmtHitCollection)[i]->GetPlaneNumber()==planeNumber && (*pmtHitCollection)[i]->GetSuperPlaneNumber()==superplaneNumber){
      hit=(*pmtHitCollection)[i];
      break;
    }
  }
  
  if(hit==NULL){//this pmt wasnt previously hit in this event
    hit = new G4PMTHit(); //so create new hit
    hit->SetPMTNumber(pmtNumber);
    hit->SetStripNumber(stripNumber);
		hit->SetPlaneNumber(planeNumber);
		hit->SetSuperPlaneNumber(superplaneNumber);
    hit->SetPMTPhysVol(physVol);
        
 		hit->SetPos(currentPostPos);
  	
    pmtHitCollection->insert(hit);
    
    //Commented and replaced above by SetPos()
    //hit->SetPMTPos((*pmtPositionsX)[pmtNumber],(*pmtPositionsY)[pmtNumber],(*pmtPositionsZ)[pmtNumber]);
    //cout<<"HIT X,Y,Z ==> "<<(*pmtPositionsX)[pmtNumber]<<"  "<<(*pmtPositionsY)[pmtNumber]<<"  "<<(*pmtPositionsZ)[pmtNumber]<<endl;
   
  }

	//get photon arrival time and energy at the photocathode
  double time = aStep->GetTrack()->GetGlobalTime()/ns;
  hit->InsertPhotonTime(time); //insert this time in the photon time vector
	double energy = aStep->GetTrack()->GetTotalEnergy()/eV;
  hit->InsertPhotonEnergy(energy); //insert this time in the photon time vector

	hit->InsertPhotonPosition(currentPostPos); //insert this position in the photon position vector		

  hit->IncPhotonCount(); //increment counts for the selected pmt

	
	if(aStep->GetTrack()->GetParentID()>0){
    if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="Scintillation") {
			hit->IncPhotonCount_scint();
			hit->SetPhotonType(1);
		}
    if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="Cerenkov") {
			hit->IncPhotonCount_cerenk();	 
			hit->SetPhotonType(2);
		} 
		if(aStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="OpWLS") {
			hit->IncPhotonCount_wls();	 
			hit->SetPhotonType(3);
		}   
	}

  hit->SetDrawit(true);
  
  
  

  return true;
}


void G4PMTSD::EndOfEvent(G4HCofThisEvent* ){
}


void G4PMTSD::clear(){
}


void G4PMTSD::DrawAll(){
} 


void G4PMTSD::PrintAll(){
} 

