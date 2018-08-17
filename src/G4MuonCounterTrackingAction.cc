/**
* @file G4MuonCounterTrackingAction.cc
* @class G4MuonCounterTrackingAction
* @brief Handle all operations to be performed at track level
*
* Manage operations before the track starts (i.e. access to primary and secondary track information, count given tracks,...), operations when the track stops.
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterTrajectory.hh"
#include "G4MuonCounterTrackingAction.hh"
#include "UserTrackInformation.hh"
#include "UserEventInformation.hh"
#include "G4MuonCounterConstruction.hh"
#include "G4MuonCounterPrimaryGenerator.hh"
#include "Utilities.hh"

#include "G4RunManager.hh"
#include "G4MuonCounterRecorderBase.hh" 

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"

#include <TVector3.h>
#include <TMath.h>

using namespace G4MuonCounterSimulatorUSC;

G4MuonCounterTrackingAction::G4MuonCounterTrackingAction(G4MuonCounterRecorderBase* r,G4MuonCounterPrimaryGenerator* gen)
  :recorder(r),generator(gen)
{

	Nscint=0;
	Ncerenk=0;
	recordEmissionInfo= true;
}


void G4MuonCounterTrackingAction::PreUserTrackingAction(const G4Track* aTrack){

	//cout<<"TrackingAction::PreUserTrackingAction"<<endl;

	UserEventInformation* eventInformation= (UserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();

  //Let this be up to the user via vis.mac
  //  fpTrackingManager->SetStoreTrajectory(true);
  
  //Use custom trajectory class
  fpTrackingManager->SetTrajectory(new G4MuonCounterTrajectory(aTrack));

  //This user track information is only relevant to the photons
  //fpTrackingManager->SetUserTrackInformation(new UserTrackInformation);


  //G4VPhysicalVolume* currentVolume = aTrack->GetStep()->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4VPhysicalVolume* currentVolume = aTrack->GetVolume();
	int currentCopyNo= currentVolume->GetCopyNo();
	G4String currentVolumeName= currentVolume->GetName();
	
	if(!aTrack->GetUserInformation()){
		fpTrackingManager->SetUserTrackInformation(new UserTrackInformation);
		//oppure 
		//UserTrackInformation* anInfo = new UserTrackInformation(aTrack);
   	//G4Track* theTrack = (G4Track*)aTrack;
    //theTrack->SetUserInformation(anInfo);

		if(aTrack->GetParentID()==0){
			G4ThreeVector direction= aTrack->GetMomentumDirection();
			G4ThreeVector momentum= aTrack->GetMomentum();

			TVector3 directionROOT= TVector3(direction.x(),direction.y(),direction.z());
			double Px= momentum.x();
			double Py= momentum.y();
			double Pz= momentum.z();
			G4ThreeVector vertex= aTrack->GetVertexPosition();
			double KinEnergy= aTrack->GetKineticEnergy();
      double Theta= generator->GetThetaOfGenParticles();
      double Phi= generator->GetPhiOfGenParticles();
           
			double time= aTrack->GetGlobalTime();
			double mass= aTrack->GetDefinition()->GetPDGMass();
			int pdgcode= aTrack->GetDefinition()->GetPDGEncoding();
			int trackId= aTrack->GetTrackID();

			eventInformation->InsertPrimaryParticleDirection(direction);
			eventInformation->InsertPrimaryParticleMomentum(momentum);
			eventInformation->InsertPrimaryParticlePosition(vertex);
			eventInformation->InsertPrimaryParticleEnergy(KinEnergy);

			eventInformation->InsertPrimaryParticleMass(mass);
			eventInformation->InsertPrimaryParticleTime(time);
			eventInformation->InsertPrimaryParticlePDGCode(pdgcode);

			
			cout<<"TRACK ID= "<<trackId<<"  PDGCode="<<pdgcode<<endl;
      cout<<"*** PRIMARY VERTEX ==> "<<"X[cm]="<<vertex.x()/cm<<"  Y[cm]="<<vertex.y()/cm<<"  Z[cm]="<<vertex.z()/cm<<endl;
      cout<<"*** PRIMARY DIRECTION ==> "<<" Theta[deg]="<<Theta/deg<<" Phi[deg]="<<Phi/deg<<"  Dx="<< direction.x()<<"  Dy="<< direction.y()<<"  Dz="<< direction.z()<<endl;
			cout<<"*** PRIMARY DIRECTION ==> "<<" Theta[deg]="<<directionROOT.Theta()*180./TMath::Pi()<<" Phi[deg]="<<directionROOT.Phi()*180./TMath::Pi()<<"  Dx="<< direction.x()<<"  Dy="<< direction.y()<<"  Dz="<< direction.z()<<endl;
      cout<<"*** PRIMARY MOMENTUM & ENERGY==> "<<" Px[GeV]="<<Px/GeV<<"  Py[GeV]="<<Py/GeV<<"  Pz[GeV]="<<Pz/GeV<<"  KinEnergy[GeV]="<<KinEnergy/GeV<<"  Mass[MeV]="<<mass/MeV<<"  PDGCode="<<pdgcode<<endl;
			
		}//if primary particle
		

	}//close if
	

	//if(aTrack->GetParentID()==1 && aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition() && currentVolumeName=="Fiber_phys"){
	if(aTrack->GetParentID()==1 && aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){	
		int vol_depth=0;
		
		if(currentVolumeName=="CoatUp_phys"||currentVolumeName=="CoatSide_phys"||currentVolumeName=="CoatBottom_phys") vol_depth=1;
		else if(currentVolumeName=="Scint_phys") vol_depth=1;
		else if(currentVolumeName=="Groove_phys") vol_depth=1;
		else if(currentVolumeName=="FiberClad2_phys") vol_depth=1;
		else if(currentVolumeName=="FiberClad1_phys") vol_depth=2;
		else if(currentVolumeName=="Fiber_phys") vol_depth=3;

		// get STRIP ID of the current track
		G4VPhysicalVolume* currentModuleVolume = aTrack->GetTouchable()->GetVolume(vol_depth);			
  	int currentModuleCopyNo= currentModuleVolume->GetCopyNo();
		
		// get PLANE ID of the current track
		G4VPhysicalVolume* currentPlaneModuleVolume = aTrack->GetTouchable()->GetVolume(vol_depth+1);			
  	int currentPlaneModuleCopyNo= currentPlaneModuleVolume->GetCopyNo();

		// get SUPER PLANE ID of the current track
		G4VPhysicalVolume* currentSuperPlaneModuleVolume = aTrack->GetTouchable()->GetVolume(vol_depth+3);//depth_level+2 without inside superplane casing			
  	int currentSuperPlaneModuleCopyNo= currentSuperPlaneModuleVolume->GetCopyNo();
			
		//get here the direction of the parent track (this must be a primary!) of this secondary
		// to be implemented...see this tips:
		//    http://geant4.slac.stanford.edu/Tips/event/1.html
		//    http://geant4.slac.stanford.edu/Tips/event/3.html
		std::vector<G4ThreeVector> primaryTrackDirectionVect= eventInformation->GetPrimaryParticleDirection();
		G4ThreeVector primaryTrackDirection= primaryTrackDirectionVect[0];//TO BE REPLACED!!!
		TVector3 primaryTrackDirection_ROOT;
		primaryTrackDirection_ROOT.SetXYZ(primaryTrackDirection.x(),primaryTrackDirection.y(),primaryTrackDirection.z());		

		//get the direction of the current secondary photon
		G4ThreeVector currentTrackDirection= aTrack->GetMomentumDirection();
		TVector3 currentTrackDirection_ROOT;
		currentTrackDirection_ROOT.SetXYZ(currentTrackDirection.x(),currentTrackDirection.y(),currentTrackDirection.z());		

		//calculate and set emission angle (angle between primary and secondary directions)
		double EmissionAngle= currentTrackDirection_ROOT.Angle(primaryTrackDirection_ROOT);
		EmissionAngle*= 180./TMath::Pi();
		
		//set emission time
		double EmissionTime= aTrack->GetGlobalTime();

		//set kin energy/wavelength
		double KinEnergy= aTrack->GetKineticEnergy();
		double EmissionWavelength=  Utilities::Energy2Wavelength(KinEnergy/eV);
		
		if(recordEmissionInfo){
			eventInformation->InsertPhotonEmissionAngle(EmissionAngle);
			eventInformation->InsertPhotonEmissionTime(EmissionTime/ns);
			eventInformation->InsertPhotonEmissionWavelength(EmissionWavelength);
		}

		//set process creator
		if(aTrack->GetCreatorProcess()->GetProcessName()=="Scintillation") {
			eventInformation->IncPhotonCount_Scint();
			if(recordEmissionInfo) eventInformation->InsertProcessType(kScintillation);
			Nscint++;
		}
		else if(aTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov"){
			eventInformation->IncPhotonCount_Ceren();
			if(recordEmissionInfo) eventInformation->InsertProcessType(kCerenkov);
			Ncerenk++;
		}
		else if(aTrack->GetCreatorProcess()->GetProcessName()=="OpWLS") {
			eventInformation->IncPhotonCount_WLS();			
			if(recordEmissionInfo) eventInformation->InsertProcessType(kWLS);
		}//close else if

		if(recordEmissionInfo){
			//insert current strip number
			eventInformation->InsertStripNumber(currentModuleCopyNo);

			//insert current plane number
			eventInformation->InsertPlaneNumber(currentPlaneModuleCopyNo);

			//insert current super plane number
			eventInformation->InsertSuperPlaneNumber(currentSuperPlaneModuleCopyNo);
		}
	
	}//close if		
	else if(aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){
		if(aTrack->GetCreatorProcess()->GetProcessName()=="OpWLS") {
			eventInformation->IncPhotonCount_WLS();			
		}//close if
	}

}


void G4MuonCounterTrackingAction::PostUserTrackingAction(const G4Track* aTrack){ 

	//cout<<"TrackingAction::PostUserTrackingAction"<<endl;

  G4MuonCounterTrajectory* trajectory=(G4MuonCounterTrajectory*)fpTrackingManager->GimmeTrajectory();
  UserTrackInformation* trackInformation=(UserTrackInformation*)aTrack->GetUserInformation();
	
	UserEventInformation* eventInformation= (UserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();



	G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries){
    UserTrackInformation* info = (UserTrackInformation*)(aTrack->GetUserInformation());
    size_t nSeco = secondaries->size();
    if(nSeco>0){
      for(size_t i=0;i<nSeco;i++){ 
        UserTrackInformation* infoNew = new UserTrackInformation(info);
        (*secondaries)[i]->SetUserInformation(infoNew);
      }//close for
    }//close if
  }//close if secondaries


  //Lets choose to draw only the photons that hit the pmt
  if(aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){
		
    const G4VProcess* creator=aTrack->GetCreatorProcess();
    if(creator && creator->GetProcessName()=="OpWLS"){
      trajectory->WLS();
			//cout<<"WLS trajectory"<<endl;
      //trajectory->SetDrawTrajectory(true);
    }
		if(creator && creator->GetProcessName()=="Cerenkov"){
      trajectory->CERENK();
      //trajectory->SetDrawTrajectory(true);
    }
		if(creator && creator->GetProcessName()=="Scintillation"){
      trajectory->SCINT();
      //trajectory->SetDrawTrajectory(true);
    }
    
		//cout<<"current track status "<<trackInformation->GetTrackStatus()<<endl;
		if(creator && creator->GetProcessName()=="OpWLS" && trackInformation->GetTrackStatus()&internalReflected) {
			trajectory->WLS(); 	
      //trajectory->SetDrawTrajectory(true);
		}
		
	/*
		if(trackInformation->GetTrackStatus()&rayleighScattered) {
			trajectory->SCATTERED(); 	
      trajectory->SetDrawTrajectory(true);
		}
	*/

    if(trackInformation->GetTrackStatus()&hitPMT) trajectory->SetDrawTrajectory(true);//draw only tracks hitting the PMTs
		//trajectory->SetDrawTrajectory(true);
    
  }
  else //draw all other trajectories
    trajectory->SetDrawTrajectory(true);
  
  

  if(trackInformation->GetForceDrawTrajectory()) trajectory->SetDrawTrajectory(true);

  if(recorder)recorder->RecordTrack(aTrack);

}


