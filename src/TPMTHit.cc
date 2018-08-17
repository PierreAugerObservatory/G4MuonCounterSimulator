/**
* @file TPMTHit.cc
* @class TPMTHit
* @brief Define the PMT hit structure for the output ROOT file
*
* A TPMTHit class is created for a pmt if a photon hits the photocathode
* A ROOT dictionary for this class is generated by the Makefile
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "TPMTHit.hh"

TPMTHit::TPMTHit()
: StripId(-1),PlaneId(-1),SuperPlaneId(-1),PMTId(-1),PhotonCounts(0),PhotonCounts_scint(0),PhotonCounts_cerenk(0),PhotonCounts_wls(0),PhotonType(-1)
{

	Position(0);

	PhotonTime.clear();
	PhotonTime.resize(0);

	PhotonEnergy.clear();
	PhotonEnergy.resize(0);

	PhotonPosition.clear();
	PhotonPosition.resize(0);
  
}//close constructor

TPMTHit::~TPMTHit(){

}//close destructor

ClassImp(TPMTHit)

#ifdef __MAKECINT__
#pragma link C++ class TPMTHit+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TPMTHit*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TPMTHit>+;
#endif
