/**
* @file G4MuonCounterSteppingAction.cc
* @class G4MuonCounterSteppingAction
* @brief Handle all operations to be performed at step level
*
* Allows to access to step information, track information, ...
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterSteppingAction.hh"
#include "G4MuonCounterEventAction.hh"
#include "G4MuonCounterTrackingAction.hh"
#include "G4MuonCounterTrajectory.hh"
#include "G4PMTSD.hh"
#include "UserTrackInformation.hh"
#include "UserEventInformation.hh"
#include "SteppingMessenger.hh"
#include "G4MuonCounterRecorderBase.hh"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4GeometryTolerance.hh"
#include "G4Box.hh"

#include <iostream>

using namespace G4MuonCounterSimulatorUSC;
using namespace std;

G4OpBoundaryProcess* G4MuonCounterSteppingAction::boundary; 

G4MuonCounterSteppingAction::G4MuonCounterSteppingAction(G4MuonCounterRecorderBase* r)
  :recorder(r),oneStepPrimaries(false)
{

  steppingMessenger = new SteppingMessenger(this);
	recordAbsorptionInfo= false;
  
}


G4MuonCounterSteppingAction::~G4MuonCounterSteppingAction(){
  
}


void G4MuonCounterSteppingAction::UserSteppingAction(const G4Step * theStep){

	
	//cout<<"G4MuonCounterSteppingAction::UserSteppingAction()"<<endl;
  //G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
  //G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  
	//## Exit if particle is out of world
	if(!theStep->GetPostStepPoint()->GetPhysicalVolume()){
    return;
  }

	//## Kill one-step tracks in scintillator
	G4Track* theTrack = theStep->GetTrack();
	trackInformation= (UserTrackInformation*)theTrack->GetUserInformation();
  eventInformation= (UserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();

	if(oneStepPrimaries && theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Scint_phys") {
		cout<<"G4MuonCounterSteppingAction::UserSteppingAction(): Killing one-step track in scintillator"<<endl;
   	theTrack->SetTrackStatus(fStopAndKill);   	
  }//close if
    

  //## Find the boundary process only once
	boundaryStatus= Undefined;
  boundary= NULL;
  if(!boundary){
    G4ProcessManager* pm = theStep->GetTrack()->GetDefinition()->GetProcessManager();
    int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    int i;
    for(i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
				boundary = (G4OpBoundaryProcess*)(*pv)[i];
				break;
      }
    }
  }//close if !boundary


	HandlePrimaryTrack(theStep);
	//CountTracksOnSurface(theStep,"PlaneHousing_phys");
	CountTracksOnSurface(theStep,"StripHousing_phys");
	
	HandleOpticalPhotonTrack(theStep);
  

  if(recorder)recorder->RecordStep(theStep);

  
}//close function




bool G4MuonCounterSteppingAction::CountTracksOnSurface(const G4Step * theStep,G4String physVolName){

	//## If current surface volume does not correspond to chosen one return false
	G4String volumeName= theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
	if(volumeName != physVolName || theStep->GetPreStepPoint()->GetStepStatus() != fGeomBoundary)	
		return false;
	
	G4String particleName= theStep->GetTrack()->GetDefinition()->GetParticleName();
	double particleWeight= theStep->GetPreStepPoint()->GetWeight();

	G4TouchableHandle theTouchable = theStep->GetPreStepPoint()->GetTouchableHandle();
	double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

	G4VSolid * solid = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
	if( solid->GetEntityType() != "G4Box" ){
		G4Exception("G4MuonCounterSteppingAction::CountTracksOnSurface(): Solid type is not supported. Use a G4Box");
		exit(1);
	}
	G4Box* boxSolid = (G4Box*)(solid);


	//## Find depth volume level in the geometry
	//## Determine plane & superplane id	
	int depth_level;
	if(volumeName=="SuperPlaneHousing_phys")
		depth_level= -2;
	else if(volumeName=="CasingInside_phys")
		depth_level= -1;
	else if(volumeName=="PlaneHousing_phys")
		depth_level= 0;
	else if(volumeName=="StripHousing_phys")
		depth_level= 1;
	else if(volumeName=="CoatUp_phys")
		depth_level= 2;
	else if(volumeName=="CoatSide_phys")
		depth_level= 2;
	else if(volumeName=="CoatBottom_phys")
		depth_level= 2;
	else if(volumeName=="Scint_phys")
		depth_level= 2;
	else if(volumeName=="Groove_phys")
		depth_level= 2;
	else if(volumeName=="FiberClad1_phys")
		depth_level= 2;
	else if(volumeName=="FiberClad2_phys")
		depth_level= 3;
	else if(volumeName=="Fiber_phys")
		depth_level= 4;
	else{
		G4Exception("G4MuonCounterSteppingAction::CountTracksOnSurface(): Invalid surface volume...exit");
		exit(1);	
	}
	

	int PlaneNo; 
	if(depth_level<0) PlaneNo= 0;
	else PlaneNo= theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(depth_level)->GetCopyNo();
	int SuperPlaneNo = theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(depth_level+2)->GetCopyNo();
	int AbsPlaneNo= 2*SuperPlaneNo + PlaneNo;			
						

	//## Counting on that surface
	//## Incoming particles   -->|   |
	G4ThreeVector GlobalPos= theStep->GetPreStepPoint()->GetPosition();
	G4ThreeVector LocalPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(GlobalPos);

	if( std::fabs(LocalPos.z()-boxSolid->GetZHalfLength()) < kCarTolerance ){
		//increment particle counts at surface		
		if(particleName=="mu-" || particleName=="mu+"){
			//cout<<"MU IN @ surface (SP,P,AP)=("<<SuperPlaneNo<<","<<PlaneNo<<","<<AbsPlaneNo<<")  nWeight="<<particleWeight<<endl;
			eventInformation->IncInMuonsAtPlaneSurface(AbsPlaneNo,+particleWeight);	
		}	
		else if(particleName=="e-" || particleName=="e+" || particleName=="gamma"){
			//cout<<"EM IN @ surface (SP,P,AP)=("<<SuperPlaneNo<<","<<PlaneNo<<","<<AbsPlaneNo<<")  nWeight="<<particleWeight<<endl;
			eventInformation->IncInEmAtPlaneSurface(AbsPlaneNo,+particleWeight);
		}
		else if(particleName=="proton" || particleName=="neutron"){
			//cout<<"HADR IN @ surface (SP,P,AP)=("<<SuperPlaneNo<<","<<PlaneNo<<","<<AbsPlaneNo<<")  nWeight="<<particleWeight<<endl;
			eventInformation->IncInHadronsAtPlaneSurface(AbsPlaneNo,+particleWeight);
		}

	}//close if incoming particles
	

	//## Outcoming particles  |<-- |
	GlobalPos= theStep->GetPostStepPoint()->GetPosition();
	LocalPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(GlobalPos);
	if(std::fabs( LocalPos.z()-boxSolid->GetZHalfLength()) < kCarTolerance ){
		//increment particle counts at surface
		if(particleName=="mu-" || particleName=="mu+"){
			//cout<<"MU OUT @ surface (SP,P,AP)=("<<SuperPlaneNo<<","<<PlaneNo<<","<<AbsPlaneNo<<")  nWeight="<<particleWeight<<endl;
			eventInformation->IncOutMuonsAtPlaneSurface(AbsPlaneNo,+particleWeight);	
		}	
		else if(particleName=="e-" || particleName=="e+" || particleName=="gamma"){
			//cout<<"EM OUT @ surface (SP,P,AP)=("<<SuperPlaneNo<<","<<PlaneNo<<","<<AbsPlaneNo<<")  nWeight="<<particleWeight<<endl;
			eventInformation->IncOutEmAtPlaneSurface(AbsPlaneNo,+particleWeight);
		}	
		else if(particleName=="proton" || particleName=="neutron"){
			//cout<<"HADR OUT @ surface (SP,P,AP)=("<<SuperPlaneNo<<","<<PlaneNo<<","<<AbsPlaneNo<<")  nWeight="<<particleWeight<<endl;
			eventInformation->IncOutHadronsAtPlaneSurface(AbsPlaneNo,+particleWeight);
		}
	}//close if outcoming particles

	
	return true;

}//close CountTracksOnSurface


bool G4MuonCounterSteppingAction::HandlePrimaryTrack(const G4Step * theStep){

	//## Check this is a primary track
	if(theStep->GetTrack()->GetParentID() != 0)	
		return false;

	G4Track* theTrack = theStep->GetTrack();
	G4TrackVector* fSecondary= fpSteppingManager->GetfSecondary();   
  int tN2ndariesTot =   fpSteppingManager->GetfN2ndariesAtRestDoIt()
                      + fpSteppingManager->GetfN2ndariesAlongStepDoIt()
                      + fpSteppingManager->GetfN2ndariesPostStepDoIt();
    
  //## If we havent already found the conversion position and there were 
  //## secondaries generated, then search for it
  if(!eventInformation->IsConvPosSet() && tN2ndariesTot>0 ){    
  	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; lp1<(*fSecondary).size(); lp1++){
			const G4VProcess* creator=(*fSecondary)[lp1]->GetCreatorProcess();
			if(creator){
	  		G4String creatorName= creator->GetProcessName();
	  		if(creatorName=="phot"||creatorName=="compt"||creatorName=="conv"){ 
	    		//since this is happening before the secondary is being tracked
	    		//the Vertex position has not been set yet(set in initial step)
	    		eventInformation->SetConvPos((*fSecondary)[lp1]->GetPosition());
	  		}//close if creatorName
			}//close if creator
     }//close for
   }//close if eventInformation

	
  
	return true;

}//close HandlePrimaryTrack()



bool G4MuonCounterSteppingAction::HandleOpticalPhotonTrack(const G4Step * theStep){

	if(theStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())	
		return false;

	//## Kill and count photons coming from the strip and entering the world or the plane container  	
	//## to speed-up the simulation
	G4Track* theTrack= theStep->GetTrack();
	G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
	G4String volumeName= thePostPV->GetName();
	if(volumeName=="expHall_phys" || volumeName=="PlaneHousing_phys" || volumeName=="SuperPlaneHousing_phys"){
		theTrack->SetTrackStatus(fStopAndKill);
		eventInformation->IncPhotonCount_killed();	
	}

	//## Count photons produced by different processes
	//## Rayleigh process?
	if(thePostPoint->GetProcessDefinedStep()->GetProcessName()=="OpRayleigh"){
		trackInformation->AddTrackStatusFlag(rayleighScattered);
	} 

	//## WLS process??
  //if(thePostPoint->GetProcessDefinedStep()->GetProcessName()=="OpWLS"){
		//cout<<"WLS emission!"<<endl;
	//}

  //## Absorption process?
  if(thePostPoint->GetProcessDefinedStep()->GetProcessName()=="OpAbsorption"){
  	eventInformation->IncAbsorption();
    trackInformation->AddTrackStatusFlag(absorbed);
			
		int depth_level=0;
		bool isAbsorbed= false;
		int whereAbsorbed= -1;
			
		if(thePostPV->GetName()=="CoatUp_phys"||thePostPV->GetName()=="CoatSide_phys"||thePostPV->GetName()=="CoatBottom_phys"){
			isAbsorbed=true;
			whereAbsorbed= kAbsInCoat;
			depth_level=1;
		}
		else if(thePostPV->GetName()=="Scint_phys"){
			isAbsorbed=true;
			whereAbsorbed= kAbsInScint;
			depth_level=1;
		}
		else if(thePostPV->GetName()=="Groove_phys"){
			isAbsorbed=true;
			whereAbsorbed= kAbsInGroove;
			depth_level=1;
		}
		else if(thePostPV->GetName()=="FiberClad2_phys"){
			isAbsorbed=true;
			whereAbsorbed= kAbsInFiber;
			depth_level=1;
		}
		else if(thePostPV->GetName()=="FiberClad1_phys"){
			isAbsorbed=true;
			whereAbsorbed= kAbsInFiber;
			depth_level=2;
		}
		else if(thePostPV->GetName()=="Fiber_phys"){
			isAbsorbed=true;
			whereAbsorbed= kAbsInFiber;
			depth_level=3;
		}
			
		// Keep count of absorption length in scintillator
		// Store track length of photons absorbed in the scintillator
		// Exclude photons absorbed at the border	
		if(isAbsorbed){
			//store track length
			double theTrackLength= theTrack->GetTrackLength();
			
			//store horizontal track length
			G4ThreeVector VertexPos= theTrack->GetVertexPosition();
			G4ThreeVector AbsPos= theTrack->GetPosition();
			double theTrackDist_x= fabs(AbsPos.x()-VertexPos.x());
        
			int currentStripNo= theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(depth_level)->GetCopyNo();
			int currentPlaneNo= theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(depth_level+1)->GetCopyNo();
			int currentSuperPlaneNo= theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(depth_level+3)->GetCopyNo();
				
				
			if(recordAbsorptionInfo) {
				eventInformation->InsertPhotonTrackLength(theTrackLength/CLHEP::cm);
				eventInformation->InsertPhotonHorizontalTrackLength(theTrackDist_x/CLHEP::cm);
				eventInformation->InsertStripNumber2(currentStripNo);
				eventInformation->InsertPlaneNumber2(currentPlaneNo);
				eventInformation->InsertSuperPlaneNumber2(currentSuperPlaneNo);
				eventInformation->InsertAbsorptionFlag(whereAbsorbed);		
			}

		}//close if absorbed
	}//close if OpAbsorption
   
		
  boundaryStatus= boundary->GetStatus();


  //## Check to see if the particle was actually at a boundary
  //## Otherwise the boundary status may not be valid
  //## Prior to Geant4.6.0-p1 this would not have been enough to check
  if(thePostPoint->GetStepStatus() == fGeomBoundary){
  	switch(boundaryStatus){
      case Absorption:
			{
				trackInformation->AddTrackStatusFlag(boundaryAbsorbed);
				eventInformation->IncBoundaryAbsorption();

				int depth_level=0;
				bool isAbsorbed= false;
				int whereAbsorbed= kAbsInBoundary;			

				if(thePostPV->GetName()=="CoatUp_phys"||thePostPV->GetName()=="CoatSide_phys"||thePostPV->GetName()=="CoatBottom_phys"){
					isAbsorbed=true;
					depth_level=1;
				}
				else if(thePostPV->GetName()=="Scint_phys"){
					isAbsorbed=true;
					depth_level=1;
				}
				else if(thePostPV->GetName()=="Groove_phys"){
					isAbsorbed=true;
					depth_level=1;
				}
				else if(thePostPV->GetName()=="FiberClad2_phys"){
					isAbsorbed=true;
					depth_level=1;
				}
				else if(thePostPV->GetName()=="FiberClad1_phys"){
					isAbsorbed=true;
					depth_level=2;
				}
				else if(thePostPV->GetName()=="Fiber_phys"){
					isAbsorbed=true;
					depth_level=3;
				}

				if(isAbsorbed){
					//store track length
					double theTrackLengthAtBoundary= theTrack->GetTrackLength();
						
					//store horizontal track length
					G4ThreeVector VertexAtBoundaryPos= theTrack->GetVertexPosition();
				 	G4ThreeVector AbsAtBoundaryPos= theTrack->GetPosition();
				 	double theTrackDistAtBoundary_x= fabs(AbsAtBoundaryPos.x()-VertexAtBoundaryPos.x());
          	
					int currentStripNo= theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(depth_level)->GetCopyNo();
					int currentPlaneNo= theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(depth_level+1)->GetCopyNo();
					int currentSuperPlaneNo= theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(depth_level+3)->GetCopyNo();
						
					if(recordAbsorptionInfo) {
						eventInformation->InsertPhotonTrackLength(theTrackLengthAtBoundary/CLHEP::cm);
						eventInformation->InsertPhotonHorizontalTrackLength(theTrackDistAtBoundary_x);
						eventInformation->InsertStripNumber2(currentStripNo);
						eventInformation->InsertPlaneNumber2(currentPlaneNo);	
						eventInformation->InsertSuperPlaneNumber2(currentSuperPlaneNo);
						eventInformation->InsertAbsorptionFlag(whereAbsorbed);	
					}	
				}//close if
					
				break;
			}//end case Absorption
      case Detection: 
      //## Note, this assumes that the volume causing detection
      //## is the photocathode because it is the only one with
	    //## non-zero efficiency
			{
	  		//## Trigger sensitive detector manually since photon is
	  		//## absorbed but status was Detection
	  		G4SDManager* SDman = G4SDManager::GetSDMpointer();
	  		G4String sdName="/Detector/PMTSD";
	  		G4PMTSD* PMT_SD = (G4PMTSD*)SDman->FindSensitiveDetector(sdName);
	  		if(PMT_SD) PMT_SD->ProcessHits_constStep(theStep,NULL);
	  		trackInformation->AddTrackStatusFlag(hitPMT);
						
	  		break;
			}//end case Detection
      case FresnelReflection:
				//cout<<"Fresnel reflection"<<endl;
			break;
			case FresnelRefraction:
				//cout<<"Fresnel refraction"<<endl;
			break;
      case TotalInternalReflection:
				//cout<<"Total internal reflection"<<endl;
				trackInformation->AddTrackStatusFlag(internalReflected);	
			break;				
      case SpikeReflection:
				//cout<<"Spike reflection"<<endl;
				trackInformation->IncReflections();
			break;
			case LobeReflection:
				//cout<<"Lobe reflection"<<endl;
			break;
			case BackScattering:
				//cout<<"Back scattering"<<endl;
			break;
			case StepTooSmall:
				//cout<<"StepTooSmall"<<endl;
				eventInformation->IncBoundaryStepTooSmall();
			break;
			case SameMaterial:
				//cout<<"SameMaterial"<<endl;
			break;
			case NoRINDEX:
				//cout<<"NoRINDEX"<<endl;
			break;
			case Undefined:
				//cout<<"Undefined"<<endl;
			break;
      default:
			break;
  	}//close switch
      
  }//close if at boundary
  

	return true;

}//close CountOpticalPhotons()


