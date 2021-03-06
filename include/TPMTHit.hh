/**
* @file TPMTHit.hh
* @class TPMTHit
* @brief Define the PMT hit structure for the output ROOT file
*
* A TPMTHit class is created for a pmt if a photon hits the photocathode
* A ROOT dictionary for this class is generated by the Makefile
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#ifndef TPMTHit_h
#define TPMTHit_h 1

#include <Rtypes.h>
#include <TObject.h>
#include <TVector3.h>

#include <iostream>
#include <vector>


class TPMTHit : public TObject {

	public:
   
		/** 
		\brief Class constructor: initialize structures.
 		*/
  	TPMTHit();
		/** 
		\brief Class destructor: free allocated memory.
 		*/
    //virtual ~TPMTHit();
		~TPMTHit();

		/** 
		\brief StripID of the current PMT
 		*/
  	int StripId;
		/** 
		\brief PlaneID of the current PMT
 		*/		
		int PlaneId;
		/** 
		\brief SuperPlaneID of the current PMT
 		*/		
		int SuperPlaneId;
		/** 
		\brief PMTID of the current PMT
 		*/
		int PMTId;
		/** 
		\brief Total number of photon counts for the current PMT
 		*/
		int PhotonCounts;
		/** 
		\brief Total number of scintillation photon counts for the current PMT
 		*/		
		int PhotonCounts_scint;
		/** 
		\brief Total number of cerenkov photon counts for the current PMT
 		*/
		int PhotonCounts_cerenk;
		/** 
		\brief Total number of WLS photon counts for the current PMT
 		*/
		int PhotonCounts_wls;
		/** 
		\brief Photon type for the current PMT
 		*/
    int PhotonType;
		/** 
		\brief Position of photon counts for the current PMT
 		*/
    TVector3 Position;
  
		/** 
		\brief Arrival time of photons hitting the photocathode for the current PMT
 		*/
  	std::vector<double> PhotonTime;
		/** 
		\brief Energy spectrum of photons hitting the photocathode for the current PMT
 		*/  		
		std::vector<double> PhotonEnergy;
		/** 
		\brief Position of photons hitting the photocathode for the current PMT
 		*/  		
		std::vector<TVector3> PhotonPosition;
		
  
  	ClassDef(TPMTHit,1)
};



#ifdef __MAKECINT__
#pragma link C++ class TPMTHit+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TPMTHit*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TPMTHit>+;
#endif


#endif
