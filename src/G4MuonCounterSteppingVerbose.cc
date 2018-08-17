/**
* @file G4MuonCounterSteppingVerbose.cc
* @class G4MuonCounterSteppingVerbose
* @brief Handle the verbosity of the stepping action
* @author S. Riggi
* @date 05/04/2010
*/


#include "G4MuonCounterSteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>

using namespace G4MuonCounterSimulatorUSC;
using namespace std;


G4MuonCounterSteppingVerbose::G4MuonCounterSteppingVerbose()
{}


G4MuonCounterSteppingVerbose::~G4MuonCounterSteppingVerbose()
{}


void G4MuonCounterSteppingVerbose::StepInfo()
{
	CopyState();
  
  int prec = cout.precision(3);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      cout << endl;    
      cout << std::setw( 5) << "#Step#"     << " "
	     << std::setw( 6) << "X"          << "    "
	     << std::setw( 6) << "Y"          << "    "  
	     << std::setw( 6) << "Z"          << "    "
	     << std::setw( 9) << "KineE"      << " "
	     << std::setw( 9) << "dEStep"     << " "  
	     << std::setw(10) << "StepLeng"     
	     << std::setw(10) << "TrakLeng" 
	     << std::setw(10) << "Volume"    << "  "
	     << std::setw(10) << "Process"   << endl;	          
    }

    cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	<< std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	<< std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
	<< "  ";

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
			int detectorId= fTrack->GetVolume()->GetCopyNo();
			G4String volumeName= fTrack->GetVolume()->GetName();
			if(volumeName=="PlaneHousing_phys" || volumeName=="SuperPlaneHousing_phys")
				cout << std::setw(10) << fTrack->GetVolume()->GetName() << " "<< detectorId;
			else
      	cout << std::setw(10) << fTrack->GetVolume()->GetName();
    } 
    else {
      cout << std::setw(10) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      cout << "  " 
        << std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    } 
    else {
      cout << "   UserLimit";
    }

    cout << endl;

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                          fN2ndariesAlongStepDoIt +
	                          fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot 
	       << "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << endl;

	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; 
                        lp1<(*fSecondary).size(); lp1++){
	  cout << "    : "
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << std::setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  cout << endl;
	}
              
	cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << endl;
      }
    }
    
  }
  cout.precision(prec);
}


void G4MuonCounterSteppingVerbose::TrackingStarted()
{

  CopyState();
	int prec = cout.precision(3);
  if( verboseLevel > 0 ){

    cout << std::setw( 5) << "Step#"      << " "
           << std::setw( 6) << "X"          << "    "
	   << std::setw( 6) << "Y"          << "    "  
	   << std::setw( 6) << "Z"          << "    "
	   << std::setw( 9) << "KineE"      << " "
	   << std::setw( 9) << "dEStep"     << " "  
	   << std::setw(10) << "StepLeng"  
	   << std::setw(10) << "TrakLeng"
	   << std::setw(10) << "Volume"     << "  "
	   << std::setw(10) << "Process"    << endl;	     

    cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	<< std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	<< std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
	<< "  ";

    if(fTrack->GetNextVolume()){
			int detectorId= fTrack->GetVolume()->GetCopyNo();
			G4String volumeName= fTrack->GetVolume()->GetName();
			if(volumeName=="PlaneHousing_phys" || volumeName=="SuperPlaneHousing_phys")
				cout << std::setw(10) << fTrack->GetVolume()->GetName() << " "<< detectorId;
			else
      	cout << std::setw(10) << fTrack->GetVolume()->GetName();
    } 
		else {
      cout << std::setw(10) << "OutOfWorld";
    }
    cout  << "    initStep" << endl;
  }
  cout.precision(prec);
}


