/**
* @file G4PMTHit.cc
* @class G4PMTHit
* @brief Define the PMT hit structure
*
* A PMT hit contains information about the PMT Id, Strip Id, Plane Id, number of hit counts, arrival time and energy of * photons at the photocathode.
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "G4PMTHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Circle.hh"

using namespace G4MuonCounterSimulatorUSC;

G4Allocator<G4PMTHit> PMTHitAllocator;

G4PMTHit::G4PMTHit()
  :pmtNumber(-1),stripNumber(-1),planeNumber(-1),superplaneNumber(-1),photons(0),photons_scint(0),photons_cerenk(0),photons_wls(0),photons_type(-1),physVol(0),drawit(false)
{}


G4PMTHit::~G4PMTHit()
{}


G4PMTHit::G4PMTHit(const G4PMTHit &right)
  : G4VHit()
{
  pmtNumber=right.pmtNumber;
  stripNumber=right.stripNumber;
	planeNumber=right.planeNumber;
	superplaneNumber=right.superplaneNumber;
  photons=right.photons;
	photons_scint=right.photons_scint;
	photons_cerenk=right.photons_cerenk;
	photons_wls=right.photons_wls;
  photons_type= right.photons_type;
  photonTime=right.photonTime;
  physVol=right.physVol;
  drawit=right.drawit;
  
  pos = right.pos;
  
}


const G4PMTHit& G4PMTHit::operator=(const G4PMTHit &right){
  pmtNumber = right.pmtNumber;
  stripNumber = right.stripNumber;
	planeNumber=right.planeNumber;
	superplaneNumber=right.superplaneNumber;
  photons= right.photons;
  photons_scint= right.photons_scint;
	photons_cerenk= right.photons_cerenk;
	photons_wls=right.photons_wls;
  photons_type= right.photons_type;
  photonTime=right.photonTime;
  physVol= right.physVol;
  drawit= right.drawit;
  
  pos = right.pos;
  
  return *this;
}


int G4PMTHit::operator==(const G4PMTHit &right) const{
	
  return (pmtNumber==right.pmtNumber && stripNumber==right.stripNumber && planeNumber==right.planeNumber && superplaneNumber==right.superplaneNumber);
}

void G4PMTHit::Draw(){

  if(drawit&&physVol){ //ReDraw only the PMTs that have hit counts > 0
    //Also need a physical volume to be able to draw anything
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager){//Make sure that the VisManager exists
    	//draw PMT
      G4VisAttributes attribs(G4Colour(1.,0.,0.));//red
      //attribs.SetForceSolid(true);
      G4RotationMatrix rot;
      if(physVol->GetRotation())//If a rotation is defined use it
				rot=*(physVol->GetRotation());
      
      G4Transform3D trans(rot,physVol->GetTranslation());//Create transform
      //pVVisManager->Draw(*physVol,attribs,trans);//Draw it
      
      G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  
  		//Draw hit
  		G4ThreeVector hitpos= GetPMTPos();
  		G4cout<<"hitpos X,Y,Z ==> "<<hitpos.x()<<"  "<<hitpos.y()<<"  "<<hitpos.z()<<G4endl;
   	 	G4Circle circle(hitpos);
    	circle.SetScreenSize(1.);
    	circle.SetFillStyle(G4Circle::filled);
    	G4Colour colour(0.,1.,0.);//green
    	G4VisAttributes hitattribs(colour);
    	circle.SetVisAttributes(hitattribs);
    	pVVisManager->Draw(circle);
  
    }
  }
}

void G4PMTHit::Print(){
}

