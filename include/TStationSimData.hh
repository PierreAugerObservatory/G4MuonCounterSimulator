/**
* @file TStationSimData.hh
* @class TStationSimData
* @brief Define the sim data structure for each muon counter station for the output ROOT file
*
* A ROOT dictionary for this class is generated by the Makefile
* @author S. Riggi
* @date 22/12/2010
*/

#ifndef TStationSimData_h
#define TStationSimData_h 1

#include "TParticleSimData.hh"
#include "TScintHit.hh"
#include "TPMTHit.hh"

#include <Rtypes.h>
#include <TObject.h>
#include <TVector3.h>

#include <iostream>
#include <vector>


class TStationSimData : public TObject {

	public:
   
		/** 
		\brief Class constructor: initialize structures.
 		*/
  	TStationSimData();
		/** 
		\brief Class destructor: free allocated memory.
 		*/
    //virtual ~TStationSimData();
		~TStationSimData();

		/** 
		\brief Search hit mode
 		*/
		enum HitSearchMode { eInFullCollection=1, eInSelectionCollection=2};
		
		enum StationSimEventTag { eNull = 0, eSinglePlane = -1, eOneFold = 1,  eTwoFold= 2, eThreeFold = 3};


		/**
		* \brief Add a sim particle data to the list 
		*/
		void AddSimParticle(TParticleSimData part) {fParticleSimDataCollection.push_back(part);}
		void AddSimParticleInTimeBin(TParticleSimData part) {fParticleSimDataCollectionInTimeBin.push_back(part);}

		/**
		* \brief Add a PMT hit to the list 
		*/
		void AddPMTHit(TPMTHit hit) {fPMTHitCollection.push_back(hit);}
		/**
		* \brief Add a Scint hit to the list 
		*/
		bool AddScintHit(TScintHit hit,HitSearchMode whereToAdd);
		/**
		* \brief Get a selected hit from given hit
		*/
		bool SelectScintHit(TScintHit hit,TScintHit& selHit);
		/**
		* \brief Is the given hit present in the chosen collection?
		*/
		bool HasScintHit(TScintHit hit, HitSearchMode whereToSearch, int& index);

		/**
		* \brief Fill track point collection
		*/
		void CreateTrackPoints();
    /** 
		\brief Tag the events according to the number of involved detector planes
 		*/
		StationSimEventTag StationEventTagger(std::vector<int> hitPlaneList);	


		/**
		* \brief Add the current number of muons hitting each plane to the list
		*/
		void AddMuonCount(std::vector<double> vect,std::string dir);
		/**
		* \brief Add the current number of em particles hitting each plane to the list
		*/
		void AddEmCount(std::vector<double> vect,std::string dir);
		/**
		* \brief Add the current number of hadrons hitting each plane to the list
		*/
		void AddHadronCount(std::vector<double> vect,std::string dir);
		

		/**
		* \brief Make PMT simulations for triggered strips of this station
		*/
		void MakePMTSimulation();


public:

		/** 
		\brief Strip number
 		*/
  	int nStrips;
		/** 
		\brief Plane number
 		*/
  	int nPlanes;

		/** 
		\brief StationID 
 		*/
  	int fId;	
		/** 
		\brief Station name 
 		*/
  	std::string fName;
		/** 
		\brief Station position in global Auger coordinates [m]
 		*/
  	TVector3 fPosition;
		/** 
		\brief Station radius in shower plane system [m]
 		*/
  	double fRadius;	
		/** 
		\brief Station phi in shower plane system [degrees]
 		*/
  	double fPhi;
		/** 
		\brief Station radius in core system [m]
 		*/
  	double fRadiusGrd;	
		/** 
		\brief Station phi in in core system [degrees]
 		*/
  	double fPhiGrd;
		
		/** 
		\brief Station sim data event tag
 		*/
		StationSimEventTag fStationSimDataTag;

		/**
		* \brief Collection of simulated particles for this station
		*/
		std::vector<TParticleSimData> fParticleSimDataCollection;
		std::vector<TParticleSimData> fParticleSimDataCollectionInTimeBin;

		/**
		* \brief Collection of ROOT scintillator hits 
		*/
		std::vector<TScintHit> fScintHitCollection;
		std::vector<TScintHit> fSelScintHitCollection;
		std::vector<TScintHit> fScintHitCollectionInTimeBin;
  
  
		/**
		* \brief Collection of ROOT PMT hits 
		*/
		std::vector<TPMTHit> fPMTHitCollection;


		/**
		* \brief Collection of track points
		*/
		std::vector<TrackPoint> fTrackPointCollection;
		/**
		* \brief Multiplicity of hit strips per detector plane
		*/
		std::vector<int> fStripMultiplicity;
		/**
		* \brief Multiplicity of hit strips per detector plane X
		*/
		std::vector<int> fStripMultiplicityX;
		/**
		* \brief Average X of hit strips per detector plane X
		*/
		std::vector<double> fAverageXHit;
		/**
		* \brief Adiacent hit strip weight per detector plane X
		*/
		std::vector<double> fRX;
		/**
		* \brief Multiplicity of hit strips per detector plane Y
		*/
		std::vector<int> fStripMultiplicityY;
		/**
		* \brief Average Y of hit strips per detector plane Y
		*/
		std::vector<double> fAverageYHit;
		/**
		* \brief Adiacent hit strip weight per detector plane Y
		*/
		std::vector<double> fRY;

		/**
		* \brief nMuons hitting each detector plane 
		*/
		std::vector<double> fMuonNumber;
		std::vector<double> fMuonAlbedoNumber;
		/**
		* \brief nEm hitting each detector plane 
		*/
		std::vector<double> fEmNumber;
		std::vector<double> fEmAlbedoNumber;
		/**
		* \brief nHadrons hitting each detector plane 
		*/
		std::vector<double> fHadronNumber;
		std::vector<double> fHadronAlbedoNumber;
		
  	double fHitSigmaX;
		double fHitSigmaY;	
		double fHitSigmaZ;

		std::vector<int> fListOfHitAbsStripId;
		std::vector<int> fSelListOfHitAbsStripId;

  	//ClassDef(TStationSimData,1)
		ClassDef(TStationSimData,2)
};

//## PMT HIT
#ifdef __MAKECINT__
#pragma link C++ class TPMTHit+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TPMTHit*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TPMTHit>+;
#endif

//## SCINT HIT
#ifdef __MAKECINT__
#pragma link C++ class TScintHit+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TScintHit*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TScintHit>+;
#endif

//## TRACK POINT
#ifdef __MAKECINT__
#pragma link C++ class TrackPoint+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TrackPoint*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TrackPoint>+;
#endif

//## TParticleSimData
#ifdef __MAKECINT__
#pragma link C++ class TParticleSimData+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TParticleSimData*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TParticleSimData>+;
#endif

//## TStationSimData
#ifdef __MAKECINT__
#pragma link C++ class TStationSimData+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TStationSimData*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TStationSimData>+;
#endif


#endif
