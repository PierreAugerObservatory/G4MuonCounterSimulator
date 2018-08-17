/**
* @file G4ScintillatorHit.cc
* @class G4ScintillatorHit
* @brief Define the scintillator hit structure
*
* A scintillator hit contains information about the Strip Id, Plane Id, number of photons per process type emitted in the scintillator, energy deposit, position and arrival time of particles hitting the strip.
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4ScintillatorHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4GeometryManager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableHandle.hh"

#include "G4Circle.hh"

using namespace G4MuonCounterSimulatorUSC;

G4Allocator<G4ScintillatorHit> ScintHitAllocator;


G4ScintillatorHit::G4ScintillatorHit()
  :physVol(0)
{
	stripNo= -1;
	planeNo= -1;
	superplaneNo= -1;
	edep=0.;
	pos= G4ThreeVector(0,0,0);
	strippos= G4ThreeVector(0,0,0);
	timehit= 0;
	trackid= -1;
	trackparttype= -1;
	direction= G4ThreeVector(0,0,0);

	photocathodeSurfacePos.clear();
	photocathodeSurfacePos.resize(0);

	//expPMTTimeHit.clear();
	//expPMTTimeHit.resize(0);

}


G4ScintillatorHit::G4ScintillatorHit(G4VPhysicalVolume* pVol)
  :physVol(pVol)
{
	stripNo= -1;
	planeNo= -1;
	superplaneNo= -1;
	edep= 0.;
	pos= G4ThreeVector(0,0,0);
	strippos= G4ThreeVector(0,0,0);
	timehit= 0;
	trackid= -1;
	trackparttype= -1;
	direction= G4ThreeVector(0,0,0);

	photocathodeSurfacePos.clear();
	photocathodeSurfacePos.resize(0);
	
	//expPMTTimeHit.clear();
	//expPMTTimeHit.resize(0);	
}


G4ScintillatorHit::~G4ScintillatorHit()
{}


G4ScintillatorHit::G4ScintillatorHit(const G4ScintillatorHit &right)
  : G4VHit()
{
  edep = right.edep;
  pos = right.pos;
	direction= right.direction;
	strippos= right.strippos;
  physVol = right.physVol;
  
  stripNo = right.stripNo;
	planeNo = right.planeNo;
	superplaneNo = right.superplaneNo;
  timehit= right.timehit;
  trackid= right.trackid;
  trackparttype= right.trackparttype;

	photocathodeSurfacePos= right.photocathodeSurfacePos;
	//expPMTTimeHit= right.expPMTTimeHit;
  
}


const G4ScintillatorHit& G4ScintillatorHit::operator=(const G4ScintillatorHit &right){
  edep = right.edep;
  pos = right.pos;	
	direction= right.direction;
	strippos= right.strippos;
  physVol = right.physVol;
  
  stripNo = right.stripNo;
	planeNo = right.planeNo;
	superplaneNo = right.superplaneNo;
  timehit= right.timehit;
  trackid= right.trackid;
  trackparttype= right.trackparttype;

	photocathodeSurfacePos= right.photocathodeSurfacePos;
	//expPMTTimeHit= right.expPMTTimeHit;
  
  return *this;
}


G4int G4ScintillatorHit::operator==(const G4ScintillatorHit& right) const{

	return (stripNo==right.stripNo && planeNo==right.planeNo && superplaneNo==right.superplaneNo && trackid==right.trackid);
}


void G4ScintillatorHit::Draw(){

	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager){
		/*
		//### Draw hit points
    G4Circle circle(pos);
    circle.SetScreenSize(1.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);//red
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
		*/

		//### Re-draw hit strips
		if(physVol){
			G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

			G4Colour stripNewColour(1.,0.,0.);//red
    	G4VisAttributes stripAttribs(stripNewColour);
			
			G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
			/*
      G4ThreeVector theLocalPoint = theNavigator->GetGlobalToLocalTransform().TransformPoint(theGlobalPoint);
      theGlobalVector = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalVector);
      theLocalVector = theNavigator->GetGlobalToLocalTransform().TransformAxis(theGlobalVector);
			*/

     	G4Navigator* aNavigator = new G4Navigator();
     	if(theNavigator->GetWorldVolume()) aNavigator->SetWorldVolume(theNavigator->GetWorldVolume());

     	G4GeometryManager* geomManager = G4GeometryManager::GetInstance();

     	if (!geomManager->IsGeometryClosed()) {
        geomManager->OpenGeometry();
        geomManager->CloseGeometry(true);
     	}

     	aNavigator->LocateGlobalPointAndSetup(strippos,0,false);

     	G4TouchableHistoryHandle fTouchable = aNavigator->CreateTouchableHistoryHandle();
			
			//G4Transform3D trans(*(physVol->GetRotation()),physVol->GetTranslation());
			G4Transform3D trans(*(fTouchable->GetRotation()),fTouchable->GetTranslation());
    	pVVisManager->Draw(*physVol,stripAttribs,trans);
		}

  }
}


void G4ScintillatorHit::Print(){
}


