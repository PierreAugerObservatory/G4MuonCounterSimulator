/**
* @file G4MuonCounterTrajectory.cc
* @class G4MuonCounterTrajectory
* @brief Handle all trajectory tracks produced in the simulation 
*
* Set visualization attributes, colours, set drawing flags, ...
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterTrajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4Trajectory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ThreeVector.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4Polymarker.hh"

#include "Randomize.hh"

using namespace G4MuonCounterSimulatorUSC;
using namespace std;

G4Allocator<G4MuonCounterTrajectory> TrajectoryAllocator;


G4MuonCounterTrajectory::G4MuonCounterTrajectory()
  :G4Trajectory(),wls(false),cerenk(false),scint(false),scattered(false),drawit(false),forceNoDraw(false),forceDraw(false)
{
  particleDefinition=0;
  changeColorTrack= false;
  OptPhotonColor = G4Colour(0.,1.,0.);//green by default
  	
}


G4MuonCounterTrajectory::G4MuonCounterTrajectory(const G4Track* aTrack)
  :G4Trajectory(aTrack),wls(false),cerenk(false),scint(false),scattered(false),drawit(false)
{
  particleDefinition=aTrack->GetDefinition();
  changeColorTrack= false;	
  OptPhotonColor = G4Colour(0.,1.,0.);//green by default
}


G4MuonCounterTrajectory::G4MuonCounterTrajectory(G4MuonCounterTrajectory &right)
  :G4Trajectory(right),wls(right.wls),scint(right.scint),cerenk(right.cerenk),scattered(right.scattered),drawit(right.drawit)
{
  particleDefinition=right.particleDefinition;
}


G4MuonCounterTrajectory::~G4MuonCounterTrajectory()
{

}


void G4MuonCounterTrajectory::DrawTrajectory(int i_mode) const{
  //Taken from G4VTrajectory and modified to select colours based on particle
  //type and to selectively eliminate drawing of certain trajectories.
	
  if(!forceDraw && (!drawit || forceNoDraw)) return;
 
  // If i_mode>=0, draws a trajectory as a polyline and, if i_mode!=0,
  // adds markers - yellow circles for step points and magenta squares
  // for auxiliary points, if any - whose screen size in pixels is
  // given by std::abs(i_mode)/1000.  E.g: i_mode = 5000 gives easily
  // visible markers.
  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) return;
 
  const double markerSize = std::abs(i_mode)/1000;
  bool lineRequired (i_mode >= 0);
  bool markersRequired (markerSize > 0.);
  
  G4Polyline trajectoryLine;
  G4Polymarker stepPoints;
  G4Polymarker auxiliaryPoints;
  
  for (int i=0; i<GetPointEntries();i++) {
    G4VTrajectoryPoint* aTrajectoryPoint = GetPoint(i);
    const std::vector<G4ThreeVector>* auxiliaries= aTrajectoryPoint->GetAuxiliaryPoints();
    if(auxiliaries){
      for(size_t iAux=0;iAux<auxiliaries->size();++iAux) {
				const G4ThreeVector pos((*auxiliaries)[iAux]);
				if(lineRequired){
	  			trajectoryLine.push_back(pos);
				}
				if(markersRequired) {
	  			auxiliaryPoints.push_back(pos);
				}
      }//close for
    }//close if auxiliaries
    
    const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
    if (lineRequired) {
      trajectoryLine.push_back(pos);
    }
    if (markersRequired) {
      stepPoints.push_back(pos);
    }
  }//close for
  
  if(lineRequired){
    G4Colour colour;
    
    if(particleDefinition==G4OpticalPhoton::OpticalPhotonDefinition()){ 
			colour= OptPhotonColor;//default color...change then according to process
      if(wls) //WLS photons are red
				colour = G4Colour(1.,0.,0.);
			if(scint) //SCINT photons are green
				colour = G4Colour(0.,1.,0.);
			if(cerenk) //CERENKOV photons are magenta 
				colour = G4Colour(1.,0.,1.); 
			if(scattered) //SCATTERED photons are blue 
				colour = G4Colour(0.,0.,1.);    
    }//close if
    else if(particleDefinition==G4Electron::ElectronDefinition()||particleDefinition==G4Positron::PositronDefinition()){
      colour = G4Colour(0.,0.,1.);//blue
		}		
		else if(particleDefinition==G4MuonPlus::MuonPlusDefinition()||particleDefinition==G4MuonMinus::MuonMinusDefinition()){
      colour = G4Colour(1.,0.,0.);//red
		}
		else if(particleDefinition==G4Gamma::GammaDefinition()){
			colour = G4Colour(1.,0.,1.);//magenta 
		}
		else{
			colour = G4Colour(0.,1.,0.);//other particles in green
		}
    
    	G4VisAttributes trajectoryLineAttribs(colour);
    	trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
    	pVVisManager->Draw(trajectoryLine);
  	}
  	if(markersRequired) {
    	auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
    	auxiliaryPoints.SetScreenSize(markerSize);
    	auxiliaryPoints.SetFillStyle(G4VMarker::filled);
    	G4VisAttributes auxiliaryPointsAttribs(G4Colour(0.,1.,1.));  // Magenta
    	auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
    	pVVisManager->Draw(auxiliaryPoints);
    
    	stepPoints.SetMarkerType(G4Polymarker::circles);
    	stepPoints.SetScreenSize(markerSize);
    	stepPoints.SetFillStyle(G4VMarker::filled);
    	G4VisAttributes stepPointsAttribs(G4Colour(1.,1.,0.));  // Yellow.
    	stepPoints.SetVisAttributes(&stepPointsAttribs);
    	pVVisManager->Draw(stepPoints);
  	}
	}


