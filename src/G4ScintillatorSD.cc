/**
* @file G4ScintillatorSD.cc
* @class G4ScintillatorSD
* @brief Define the scintillator sensitive detector, create and handle the scintillator hits
*
* When a particle hits the scintillator, a corresponding hit is created only when the energy deposit is >0. If the hit already exists, the current information is added to the previous hit. 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4ScintillatorSD.hh"
#include "G4ScintillatorHit.hh"
#include "PMTSimulator.hh"
#include "G4MuonCounterConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

using namespace G4MuonCounterSimulatorUSC;
using namespace std;

G4ScintillatorSD::G4ScintillatorSD(G4String name)
  :G4VSensitiveDetector(name)  
{
  collectionName.insert("scintHitCollection");

}



G4ScintillatorSD::~G4ScintillatorSD()
{}


void G4ScintillatorSD::Initialize(G4HCofThisEvent* HCE){

  scintHitCollection = new ScintHitsCollection(SensitiveDetectorName,collectionName[0]); 
  //A way to keep all the hits of this event in one place if needed
  static int HCID = -1;
  if(HCID<0){ 
    HCID = GetCollectionID(0); 
  }
  HCE->AddHitsCollection( HCID, scintHitCollection );
}


bool G4ScintillatorSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

	//######################################################################################
	//##   Hit definition: track in strip with Edep!=0
	//##                       
	//###  - if current Hit already exists in the collection (same trackId&volume), simply add energy dep info
	//###  - otherwise create a new Hit, add all info and put in the collection
	//###  - tracks with same Id hitting different strips are considered as distinct hits
	//######################################################################################

  double currentEdep = aStep->GetTotalEnergyDeposit();
  //cout<<"currentEdep="<<currentEdep/MeV<<"  vol name="<<aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()<<endl;

  if(currentEdep == 0.) return false;//No edep so don't count as hit
	//if(currentEdep< fScintillatorEnergyThreshold) return false;//edep below threshold so don't count as hit


  G4ScintillatorHit* aHit;
  int nHit = scintHitCollection->entries();
  G4ThreeVector currentPrePos = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector currentPostPos = aStep->GetPostStepPoint()->GetPosition();
  
  
  // get volume of the current step
  G4VPhysicalVolume* currentVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	int currentCopyNo= currentVolume->GetCopyNo();
	G4String currentVolumeName= currentVolume->GetName();
	
	
	if(currentVolumeName!= "Scint_phys") return false;//consider only hit in scintillator, skip other volumes
	//if(currentVolumeName!= "Fiber_phys") return false;//consider only hit in fiber core, skip other volumes
 
	
	G4ThreeVector currentVolumeLocalOrigin(0,0,0);
  G4ThreeVector currentVolumeGlobalOrigin = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().Inverse().TransformPoint(currentVolumeLocalOrigin);

	//Correct precision problem, probably in the determination of the inverse transformation
	double precisionThreshold= 1.E-10*cm;
	if( fabs(currentVolumeGlobalOrigin.x())<precisionThreshold) currentVolumeGlobalOrigin.setX(0);
	if( fabs(currentVolumeGlobalOrigin.y())<precisionThreshold) currentVolumeGlobalOrigin.setY(0);	
	if( fabs(currentVolumeGlobalOrigin.z())<precisionThreshold) currentVolumeGlobalOrigin.setZ(0);


/*
	//alternative method
	//get traslation and rotation of current strip volume
	G4ThreeVector currentVolumeLocalToGlobalTraslation= aStep->GetPreStepPoint()->GetTouchableHandle()->GetTranslation();
	G4RotationMatrix* currentVolumeLocalToGlobalRotation= aStep->GetPreStepPoint()->GetTouchableHandle()->GetRotation();

  G4AffineTransform LocalToGlobalTransform = G4AffineTransform(currentVolumeLocalToGlobalRotation,currentVolumeLocalToGlobalTraslation);
  //LocalToGlobalTransform.Invert();
  G4ThreeVector localP = LocalToGlobalTransform.TransformPoint(position);
*/

	//## access to current strip module
	G4VPhysicalVolume* currentModuleVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1);
	int currentModuleCopyNo= currentModuleVolume->GetCopyNo();
	G4String currentModuleVolumeName= currentModuleVolume->GetName();
	//get transform up to its mother volume (the plane)
	G4ThreeVector ModuleToPlaneTranslation= currentModuleVolume->GetTranslation();
	G4RotationMatrix* ModuleToPlaneRotation= currentModuleVolume->GetRotation();
  G4AffineTransform ModuleToPlaneTransform = G4AffineTransform(ModuleToPlaneRotation,ModuleToPlaneTranslation);
  
	//## access to current plane module
	G4VPhysicalVolume* currentPlaneModuleVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2);
	int currentPlaneModuleCopyNo= currentPlaneModuleVolume->GetCopyNo();
	G4String currentPlaneModuleVolumeName= currentPlaneModuleVolume->GetName();
	//get transform up to its mother volume (the inside casing)
	G4ThreeVector PlaneToInsideCasingTranslation= currentPlaneModuleVolume->GetTranslation();
	G4RotationMatrix* PlaneToInsideCasingRotation= currentPlaneModuleVolume->GetRotation();
  G4AffineTransform PlaneToInsideCasingTransform = G4AffineTransform(PlaneToInsideCasingRotation,PlaneToInsideCasingTranslation);

	//## access to current inside casing
	G4VPhysicalVolume* currentInsideCasingVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(3);
	int currentInsideCasingCopyNo= currentInsideCasingVolume->GetCopyNo();
	G4String currentInsideCasingVolumeName= currentInsideCasingVolume->GetName();
	//get transform up to its mother volume (the superplane)
	G4ThreeVector InsideCasingToSuperPlaneTranslation= currentInsideCasingVolume->GetTranslation();
	G4RotationMatrix* InsideCasingToSuperPlaneRotation= currentInsideCasingVolume->GetRotation();
  G4AffineTransform InsideCasingToSuperPlaneTransform = G4AffineTransform(InsideCasingToSuperPlaneRotation,InsideCasingToSuperPlaneTranslation);

	//## access to current superplane module
	G4VPhysicalVolume* currentSuperPlaneModuleVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(4);//GetVolume(3) without inside superplane casing
	int currentSuperPlaneModuleCopyNo= currentSuperPlaneModuleVolume->GetCopyNo();
	G4String currentSuperPlaneModuleVolumeName= currentSuperPlaneModuleVolume->GetName();
	//get transform up to its mother volume (the world)
	G4ThreeVector SuperPlaneToWorldTranslation= currentSuperPlaneModuleVolume->GetTranslation();
	G4RotationMatrix* SuperPlaneToWorldRotation= currentSuperPlaneModuleVolume->GetRotation();
  G4AffineTransform SuperPlaneToWorldTransform = G4AffineTransform(SuperPlaneToWorldRotation,SuperPlaneToWorldTranslation);


	//## access to current PMTs
	int Ndaughters= currentModuleVolume->GetLogicalVolume()->GetNoDaughters();
	std::vector<G4ThreeVector> cathodeSurfacePos;	
	cathodeSurfacePos.clear();
	cathodeSurfacePos.resize(0);
	
	std::vector<double> hitCathodeSurfaceDistance;	
	hitCathodeSurfaceDistance.clear();
	hitCathodeSurfaceDistance.resize(0);

	//std::vector<double> hitExpPMTTime;	
	//hitExpPMTTime.clear();
	//hitExpPMTTime.resize(0);

	//search PMTs volumes among daughters	
	for(int k=0;k<Ndaughters;k++){
		G4VPhysicalVolume* thisDaughterPhysVol= currentModuleVolume->GetLogicalVolume()->GetDaughter(k);
		G4LogicalVolume* thisDaughterLogVol= thisDaughterPhysVol->GetLogicalVolume();
		G4String thisDaughterName= thisDaughterPhysVol->GetName();
		int thisDaughterCopyNo= thisDaughterPhysVol->GetCopyNo();
		
		if(thisDaughterName=="PMT_phys"){

			//define affine coordinates
			//## pmt->module
			G4ThreeVector PMTToModuleTranslation= thisDaughterPhysVol->GetTranslation();
			G4RotationMatrix* PMTToModuleRotation= thisDaughterPhysVol->GetRotation();
			G4AffineTransform PMTToModuleTransform = G4AffineTransform(PMTToModuleRotation,PMTToModuleTranslation);
			
			//get access to photocathode		
			int Ndaughters_pmt= thisDaughterLogVol->GetNoDaughters();
			for(int l=0;l<Ndaughters_pmt;l++){
				G4VPhysicalVolume* thisPMTDaughterPhysVol= thisDaughterLogVol->GetDaughter(l);
				G4LogicalVolume* thisPMTDaughterLogVol= thisPMTDaughterPhysVol->GetLogicalVolume();
				G4VSolid* thisPMTDaughterSolid= thisPMTDaughterLogVol->GetSolid();

				if(thisPMTDaughterPhysVol->GetName()=="Photocathode_phys"){
					//define affine coordinates
					//## photocathode->pmt
					G4ThreeVector PhotocathodeToPMTTranslation= thisPMTDaughterPhysVol->GetTranslation();
					G4RotationMatrix* PhotocathodeToPMTRotation= thisPMTDaughterPhysVol->GetRotation();
					//get Photocathode position in PMT coordinates
  				G4AffineTransform PhotocathodeToPMTTransform = G4AffineTransform(PhotocathodeToPMTRotation,PhotocathodeToPMTTranslation);

					double photocathode_offset;
					//cout<<"thisPMTDaughterSolid ==>"<<thisPMTDaughterSolid->GetEntityType()<<endl;
					G4Box* photocathode_box= (G4Box*)(thisPMTDaughterLogVol->GetSolid());
					
					double photocathode_halflength= photocathode_box->GetXHalfLength();
					G4ThreeVector photocathodeSurfacePosInPhotocathode;

					if(thisDaughterCopyNo==0 && currentPlaneModuleCopyNo==0){//right PMT && X plane
						photocathode_offset= +photocathode_halflength;
						photocathodeSurfacePosInPhotocathode= G4ThreeVector(photocathode_offset,0,0);
					}
					else if(thisDaughterCopyNo==1 && currentPlaneModuleCopyNo==0){//left PMT && X plane
						photocathode_offset= -photocathode_halflength;
						photocathodeSurfacePosInPhotocathode= G4ThreeVector(photocathode_offset,0,0);
					}
					else if(thisDaughterCopyNo==0 && currentPlaneModuleCopyNo==1){//right PMT && Y plane
						photocathode_offset= +photocathode_halflength;
						photocathodeSurfacePosInPhotocathode= G4ThreeVector(0,photocathode_offset,0);
					}
					else if(thisDaughterCopyNo==1 && currentPlaneModuleCopyNo==1){//left PMT && Y plane
						photocathode_offset= -photocathode_halflength;
						photocathodeSurfacePosInPhotocathode= G4ThreeVector(0,photocathode_offset,0);
					}
					else{
						cerr<<"ERROR in G4ScintillatorSD: Cannot calculate photocathode surface position...exit"<<endl;
						exit(1);
					}

					//cout<<"photocathodeSurfacePosInPhotocathode="<<photocathodeSurfacePosInPhotocathode<<endl;

					//get photocathode surface position in PMT coordinates
					G4ThreeVector photocathodeSurfacePosInPMT = PhotocathodeToPMTTransform.TransformPoint(photocathodeSurfacePosInPhotocathode);
					//cout<<"photocathodeSurfacePosInPMT="<<photocathodeSurfacePosInPMT/cm<<endl;
			
					//get photocathode surface position in module coordinates
  				G4ThreeVector photocathodeSurfacePosInModule = PMTToModuleTransform.TransformPoint(photocathodeSurfacePosInPMT);
					//cout<<"photocathodeSurfacePosInModule="<<photocathodeSurfacePosInModule/cm<<endl;

					//get photocathode surface position in plane coordinates
  				G4ThreeVector photocathodeSurfacePosInPlane = ModuleToPlaneTransform.TransformPoint(photocathodeSurfacePosInModule);
					//cout<<"photocathodeSurfacePosInPlane="<<photocathodeSurfacePosInPlane/cm<<endl;

					//get photocathode surface position in inside casing coordinates
  				G4ThreeVector photocathodeSurfacePosInCasing = PlaneToInsideCasingTransform.TransformPoint(photocathodeSurfacePosInPlane);
					//cout<<"photocathodeSurfacePosInCasing="<<photocathodeSurfacePosInCasing/cm<<endl;

					//get photocathode surface position in superplane coordinates
  				G4ThreeVector photocathodeSurfacePosInSuperPlane = InsideCasingToSuperPlaneTransform.TransformPoint(photocathodeSurfacePosInCasing);
					//cout<<"photocathodeSurfacePosInSuperPlane="<<photocathodeSurfacePosInSuperPlane/cm<<endl;

					//get photocathode surface position in world coordinates
  				G4ThreeVector photocathodeSurfacePosInWorld = SuperPlaneToWorldTransform.TransformPoint(photocathodeSurfacePosInSuperPlane);
					//cout<<"photocathodeSurfacePosInWorld="<<photocathodeSurfacePosInWorld/cm<<endl;

					cathodeSurfacePos.push_back(photocathodeSurfacePosInWorld);
		
					/*
					//######################################################
					//## Calculate expected hit time at the photocathode 
					//######################################################
					G4ThreeVector HitPos = currentPrePos+currentPostPos;
					G4ThreeVector HitPMTDistance= photocathodeSurfacePosInWorld-HitPos;
					double HitPMTHorizontalDistance= 0.;
					int SuperPlaneId= currentSuperPlaneModuleCopyNo;
					int PlaneId= currentPlaneModuleCopyNo;
					int AbsPlaneId= 2*SuperPlaneId+PlaneId;
					if(AbsPlaneId%2==0){//X plane: take distance in x
						HitPMTHorizontalDistance= fabs(HitPMTDistance.x());
					}
					else if(AbsPlaneId%2==1){//Y plane: take distance in y
						HitPMTHorizontalDistance= fabs(HitPMTDistance.y());
					}

					// add real distance to the PMT!!!
					//HitPMTHorizontalDistance+= ...

					PMTSimulator* fPMTSimulator= new PMTSimulator();
					fPMTSimulator->SetTransitTime(G4MuonCounterConstruction::GetPMTTransitTimeAverage());
					fPMTSimulator->SetTransitTimeSpread(G4MuonCounterConstruction::GetPMTTransitTimeSpread());
					fPMTSimulator->SetScintillatorDecayTime(G4MuonCounterConstruction::GetStripDecayTime());
					fPMTSimulator->SetFiberDecayTime(G4MuonCounterConstruction::GetFiberDecayTime());
					fPMTSimulator->SetFiberCriticalAngle(G4MuonCounterConstruction::GetFiberCriticalAngle());
					fPMTSimulator->SetFiberCoreRefractiveIndex(G4MuonCounterConstruction::GetFiberCoreRefractiveIndex());
					fPMTSimulator->SetPhotoelectronYield(G4MuonCounterConstruction::GetPMTPhotoelectronYield()); 
					double expPMTHitTime= aStep->GetPreStepPoint()->GetGlobalTime() +
																fPMTSimulator->GeneratePETime(HitPMTHorizontalDistance);
					cout<<"Generated PMTExpTime pmt "<<l+1<<"  ="<<expPMTHitTime<<endl;
					hitExpPMTTime.push_back(expPMTHitTime);
					delete fPMTSimulator;
					*/
				}//close if photocathode
			}//close loop PMT daughters
				
		}//close if PMT daughter
	}//close loop daughters

	/*
	cout<<"******** SCINT ****************"<<endl;
	cout<<"Scintillator Hit: currentCenter="<<currentVolumeGlobalOrigin.x()/cm<<"  "<<currentVolumeGlobalOrigin.y()/cm<<"  "<<currentVolumeGlobalOrigin.z()/cm<<endl;
	cout<<"Scintillator Hit: currentPrePos="<<currentPrePos/cm<<"  currentPostPos="<<currentPostPos/cm<<endl;
	cout<<"Scintillator Hit: currentCopyNo="<<currentCopyNo<<"  currentVolName="<<currentVolumeName<<endl;
	cout<<"Scintillator Hit: currentModuleCopyNo="<<currentModuleCopyNo<<"  currentModuleVolName="<<currentModuleVolumeName<<endl;
	cout<<"Scintillator Hit: currentPlaneModuleCopyNo="<<currentPlaneModuleCopyNo<<"  currentPlaneModuleVolName="<<currentPlaneModuleVolumeName<<endl;
	cout<<"Scintillator Hit: currentSuperPlaneModuleCopyNo="<<currentSuperPlaneModuleCopyNo<<"  currentSuperPlaneModuleVolName="<<currentSuperPlaneModuleVolumeName<<endl;
	cout<<"Scintillator Hit: currentTrackId="<<aStep->GetTrack()->GetTrackID()<<endl;
	cout<<"******************************"<<endl;
	*/

	// get id of the current track
	G4Track* currentTrack= aStep->GetTrack();
	int currentTrackId = currentTrack->GetTrackID();
	G4ParticleDefinition* ParticleDef= currentTrack->GetDefinition();
	G4String currentParticleName= ParticleDef->GetParticleName();
	
  for(int i=0;i<nHit;i++){
    aHit = (*scintHitCollection)[i];
    G4ThreeVector HitPos = aHit->GetPos();
    
    int HitCopyNo= aHit->GetSciCopyNo();
		int HitPlaneCopyNo= aHit->GetPlaneCopyNo();
		int HitSuperPlaneCopyNo= aHit->GetSuperPlaneCopyNo();
    int HitTrackId= aHit->GetTrackId();
    G4String HitVolName= aHit->GetPhysV()->GetName();
    
    //check if the current track has already be stored as hit (check trackid, volumename & vol copyNo)
    if(currentTrackId==HitTrackId && currentVolumeName==HitVolName && currentModuleCopyNo==HitCopyNo && currentPlaneModuleCopyNo==HitPlaneCopyNo && currentSuperPlaneModuleCopyNo==HitSuperPlaneCopyNo) {
    	//track found...append to previous hit info
    	aHit->AddEdep(currentEdep);
    	return true;
    }//close if
  }//end loop hits

	
	//Create a new hit
  aHit = new G4ScintillatorHit();
  aHit->SetEdep(currentEdep);
  G4ThreeVector currentPos = currentPrePos+currentPostPos;
  currentPos/=2.;//take average position between pre and post step positions
  aHit->SetPos(currentPos);
	aHit->SetStripPos(currentVolumeGlobalOrigin);
  aHit->SetSciCopyNo(currentModuleCopyNo);//insert here module copy number (strip copy number is 0)
	aHit->SetPlaneCopyNo(currentPlaneModuleCopyNo);//insert here plane module copy number
	aHit->SetSuperPlaneCopyNo(currentSuperPlaneModuleCopyNo);//insert here plane module copy number

  aHit->SetTimeHit(aStep->GetPreStepPoint()->GetGlobalTime());
  aHit->SetTrackId(aStep->GetTrack()->GetTrackID());
	aHit->SetTrackDirection(aStep->GetTrack()->GetMomentumDirection());

	G4ThreeVector currentTrackDir= aStep->GetTrack()->GetMomentumDirection();
	//cout<<"current TrackDir =("<<currentTrackDir.x()<<","<<currentTrackDir.y()<<","<<currentTrackDir.z()<<")"<<endl;

  //set kind of hit:  1 ==> e+,e-,gamma em hit
  //                  2 ==> mu+,mu- muon hit
  //									3 ==> optical photon
  //									4 ==> other particles??	
  if(currentParticleName=="mu-"||currentParticleName=="mu+") aHit->SetTrackPartType(2);
  else if(currentParticleName=="e-"||currentParticleName=="e+"||currentParticleName=="gamma") aHit->SetTrackPartType(1);
  else if(currentParticleName=="opticalphoton") aHit->SetTrackPartType(3);
  else{
  	//cout<<" Hit by other particles "<<endl;
  	aHit->SetTrackPartType(4);
  }
  aHit->SetPhysV(currentVolume);
	aHit->SetPhotocathodeSurfacePos(cathodeSurfacePos);
  //aHit->SetExpPMTTime(hitExpPMTTime);

  scintHitCollection->insert(aHit);

  return true;
}//close function


void G4ScintillatorSD::EndOfEvent(G4HCofThisEvent* ){
}


void G4ScintillatorSD::clear(){
} 

void G4ScintillatorSD::DrawAll(){
} 


void G4ScintillatorSD::PrintAll(){
} 

