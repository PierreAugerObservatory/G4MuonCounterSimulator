/**
* @file TScintHit.hh
* @class TScintHit
* @brief Define the scintillator hit structure for the output ROOT file
*
* A TScintHit class is created for a strip if the energy deposition is >0 for that strip
* A ROOT dictionary for this class is generated by the Makefile
* @author Dr. Simone Riggi
* @date 05/04/2010
*/


#ifndef TScintHit_h
#define TScintHit_h 1

#include <Rtypes.h>
#include <TObject.h>
#include <TVector3.h>

#include <iostream>
#include <vector>


class TScintHit : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/   
  	TScintHit();
		/** 
		\brief Class destructor: free allocated memory.
 		*/
    //virtual ~TScintHit();
		~TScintHit();

		/** 
		\brief Append info from other hit
 		*/
		void Append(TScintHit aNewHit);

		//these vectors loop over Nphotons produced in the strip
		/** 
		\brief Tracklength of photons absorbed in the strip
 		*/
  	std::vector<double> PhotonTrackLength;
		/** 
		\brief Horizontal tracklength of photons absorbed in the strip
 		*/		
		std::vector<double> PhotonTrackDistance;
		/** 
		\brief Emission angles of photons in the strip
 		*/
		std::vector<double> PhotonEmissionAngle;
		/** 
		\brief Emission times of photons in the strip
 		*/
    std::vector<double> PhotonEmissionTime;
		/** 
		\brief Emission wavelength of photons in the strip
 		*/
		std::vector<double> PhotonEmissionWavelength;
		/** 
		\brief Process responsible for photon creation in the strip
 		*/
		std::vector<int> PhotonProcessType;
		
		//these vectors loop over the total number of hits for the current strip
		/** 
		\brief Energy deposit of all hits produced in the current strip
 		*/		
		std::vector<double> Edep;
		/** 
		\brief Time of all hits produced in the current strip
 		*/
		std::vector<double> Time;
		/** 
		\brief Exp Time at the PMT anode of all hits produced in the current strip for each PMT
 		*/
		std::vector< std::vector<double> > ExpPMTTime;
		/** 
		\brief TrackID of all hits produced in the current strip
 		*/
  	std::vector<int> TrackId;
		/** 
		\brief Particle type of all hits produced in the current strip
 		*/
  	std::vector<int> ParticleType; 
		/** 
		\brief Position of all hits produced in the current strip
 		*/
		std::vector<TVector3> Position; 
		/** 
		\brief Track direction of all hits produced in the current strip
 		*/
		std::vector<TVector3> TrackDirection;
		

		//Summary vars
		/** 
		\brief Position of the current strip
 		*/
		TVector3 StripPosition;
		/** 
		\brief Position of the PMT photocathode surface for the current strip
 		*/
		std::vector<TVector3> PhotocathodeSurfacePosition;
		/** 
		\brief StripID of the current strip
 		*/
		int StripId;
		/** 
		\brief PlaneID of the current strip
 		*/
		int PlaneId;
		/** 
		\brief Super PlaneID of the current strip
 		*/
		int SuperPlaneId;
		/** 
		\brief Total energy deposit in the current strip
 		*/
		double Etot;

		/** 
		\brief Total number of photons produced in the current strip
 		*/
		int Nphotons;
		/** 
		\brief Total number of scintillation photons produced in the current strip
 		*/	  
		int Nphotons_scint;
		/** 
		\brief Total number of cerenkov photons produced in the current strip
 		*/
  	int Nphotons_cerenk;
		/** 
		\brief Total number of WLS photons produced in the current strip
 		*/
		int Nphotons_wls;
		/** 
		\brief Total number of photons absorbed in the current strip
 		*/
  	int Nphotons_abs;
		/** 
		\brief Total number of photons absorbed at boundaries in the current strip
 		*/
		int Nphotons_absBoundary;

  
  	ClassDef(TScintHit,2)
		//ClassDef(TScintHit,1)
};



#ifdef __MAKECINT__
#pragma link C++ class TScintHit+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TScintHit*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TScintHit>+;
#endif


#endif