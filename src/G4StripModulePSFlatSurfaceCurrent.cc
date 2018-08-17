/**
* @file G4StripModulePSFlatSurfaceCurrent.cc
* @class G4StripModulePSFlatSurfaceCurrent
* @brief Define a user G4PSFlatSurfaceCurrent
* This is a primitive scorer class for scoring Surface Current.
*  Current version assumes only for G4Box shape, and the surface
*  is defined at the -Z plane of the box.
*  The current is given in the unit of area.
*    e.g.  (Number of tracks)/mm2.
*
* Surface is defined at the -Z surface.
* Direction                  -Z   +Z
*   0  IN || OUT            ->|<-  |      fCurrent_InOut
*   1  IN                   ->|    |      fCurrent_In
*   2  OUT                    |<-  |      fCurrent_Out
*
* @author S. Riggi
* @date 30/01/2011
*/


#include "G4StripModulePSFlatSurfaceCurrent.hh"
#include <G4PSFlatSurfaceCurrent.hh>

#include <G4MuonCounterSimulator.hh>

using namespace std;
using namespace G4MuonCounterSimulatorUSC;


G4StripModulePSFlatSurfaceCurrent::G4StripModulePSFlatSurfaceCurrent(G4String name, int direction)
  :G4PSFlatSurfaceCurrent(name,direction)
{;}

G4StripModulePSFlatSurfaceCurrent::~G4StripModulePSFlatSurfaceCurrent()
{;}

int G4StripModulePSFlatSurfaceCurrent::GetIndex(G4Step* aStep){

	int StripId= aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetCopyNo();
	int PlaneId = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetCopyNo();
	int SuperPlaneId= aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(4)->GetCopyNo();
	int AbsPlaneId= 2*SuperPlaneId + PlaneId;
	int nStrips= (G4MuonCounterSimulator::GetCurrentMuonDetector())->GetNumberOfStrips();
  int AbsStripId= AbsPlaneId*nStrips+StripId;

	int index= AbsStripId;
	cout<<"GetIndex: (SP,P,S,AP,AS)=("<<SuperPlaneId<<","<<PlaneId<<","<<StripId<<","<<AbsPlaneId<<","<<AbsStripId<<")"<<endl;

  return index;
}


