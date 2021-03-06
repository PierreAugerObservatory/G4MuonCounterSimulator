/**
* @file G4StripModulePSFlatSurfaceCurrent.hh
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

#ifndef _G4MuonCounterSimulatorUSC_G4StripModulePSFlatSurfaceCurrent_h_
#define _G4MuonCounterSimulatorUSC_G4StripModulePSFlatSurfaceCurrent_h_ 1

#include "G4PSFlatSurfaceCurrent.hh"

namespace G4MuonCounterSimulatorUSC {

class G4StripModulePSFlatSurfaceCurrent : public G4PSFlatSurfaceCurrent
{
   public: // with description
   	//G4StripModulePSFlatSurfaceCurrent(G4String name, int direction,int nx,int ny,int nz);
		G4StripModulePSFlatSurfaceCurrent(G4String name, int direction);
    virtual ~G4StripModulePSFlatSurfaceCurrent();

  protected: // with description
  	virtual int GetIndex(G4Step*);

  private:
      
			
};

}//close namespace

#endif

