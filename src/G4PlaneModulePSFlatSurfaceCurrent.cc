/**
* @file G4PlaneModulePSFlatSurfaceCurrent.cc
* @class G4PlaneModulePSFlatSurfaceCurrent
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


#include "G4PlaneModulePSFlatSurfaceCurrent.hh"
#include "G4PSFlatSurfaceCurrent.hh"

#include <G4MuonCounterSimulator.hh>

using namespace std;
using namespace G4MuonCounterSimulatorUSC;


G4PlaneModulePSFlatSurfaceCurrent::G4PlaneModulePSFlatSurfaceCurrent(G4String name, int direction)
  :G4PSFlatSurfaceCurrent(name,direction)
{;}

G4PlaneModulePSFlatSurfaceCurrent::~G4PlaneModulePSFlatSurfaceCurrent()
{;}

int G4PlaneModulePSFlatSurfaceCurrent::GetIndex(G4Step* aStep){

	int PlaneModuleCopyNo = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetCopyNo();
	int SuperPlaneModuleCopyNo= aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetCopyNo();
	int AbsPlaneModuleCopyNo= 2*SuperPlaneModuleCopyNo + PlaneModuleCopyNo;
  
	int index= AbsPlaneModuleCopyNo;

  return index;
}


