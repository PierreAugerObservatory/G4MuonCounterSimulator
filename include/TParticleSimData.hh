/**
* @file TParticleSimData.hh
* @class TParticleSimData
* @brief Define the sim data structure for each muon counter station for the output ROOT file
*
* A ROOT dictionary for this class is generated by the Makefile
* @author Simone Riggi
* @date 22/12/2010
*/

#ifndef TParticleSimData_h
#define TParticleSimData_h 1

#include "TScintHit.hh"
#include "TPMTHit.hh"
#include "TrackPoint.hh"
#include "MuonDetector.hh"


#include <Rtypes.h>
#include <TObject.h>
#include <TVector3.h>

#include <iostream>
#include <vector>


class TParticleSimData : public TObject {

	public:
   
		/** 
		\brief Class constructor: initialize structures.
 		*/
  	TParticleSimData();
		/** 
		\brief Class destructor: free allocated memory.
 		*/
    //virtual ~TParticleSimData();
		~TParticleSimData();

		enum ParticleSimEventTag { eNull = 0, eSinglePlane = -1, eOneFold = 1,  eTwoFold= 2, eThreeFold = 3};

		/**
		* \brief Add a PMT hit to the list 
		*/
		void AddPMTHit(TPMTHit hit) {fPMTHitCollection.push_back(hit);}
		/**
		* \brief Add a Scint hit to the list 
		*/
		void AddScintHit(TScintHit hit) {fScintHitCollection.push_back(hit);}
		/**
		* \brief Add a Scint hit to the selected list 
		*/
		void AddSelScintHit(TScintHit hit) {fSelScintHitCollection.push_back(hit);}
    /**
		* \brief Fill track point collection
		*/
		void CreateTrackPoints();
    /** 
		\brief Tag the events according to the number of involved detector planes
 		*/
		ParticleSimEventTag ParticleEventTagger(std::vector<int> hitPlaneList);	
		/** 
		\brief Tag the events according to the geometrical trajectory
 		*/
		int GeomEventTagger();

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
		\brief Y Separation among strips
 		*/
  	double StripSeparation;

		/** 
		\brief Particle ID 
 		*/
  	int fId;
		/** 
		\brief Particle name
 		*/
  	std::string fName;
		/** 
		\brief Particle PDGCode 
 		*/
  	int fPDGCode;
		/** 
		\brief Particle position in station coordinate system [m]
 		*/
  	TVector3 fPosition;
		/** 
		\brief Particle position in Auger coordinate system [m]
 		*/
  	TVector3 fGlobalPosition;
		/** 
		\brief Particle position radius in shower coordinate system [m]
 		*/
		double fRadius;
		/** 
		\brief Particle position radius in core coordinate system [m]
 		*/
		double fRadiusGrd;
		/** 
		\brief Particle position phi in shower coordinate system [deg]
 		*/
		double fPsi;
		/** 
		\brief Particle position phi in core coordinate system [deg]
 		*/
		double fPsiGrd;

		/** 
		\brief Particle direction
 		*/
  	TVector3 fDirection;
		/** 
		\brief Particle momentum [MeV]
 		*/
  	TVector3 fMomentum;
		/** 
		\brief Particle kinetic energy [MeV]
 		*/
  	double fEnergy;
		/** 
		\brief Particle theta [degrees]
 		*/
  	double fTheta;
		/** 
		\brief Particle phi [degrees]
 		*/
  	double fPhi;
		/** 
		\brief Particle time [ns]
 		*/
  	double fTime;
		/** 
		\brief Particle sim data event tag
 		*/
		ParticleSimEventTag fParticleSimDataTag;
		int fParticleSimDataGeomTag;
		int fHasGeomCrossedFirstPlane;

		double fHitSigmaX;
		double fHitSigmaY;	
		double fHitSigmaZ;

		MuonDetector* fDetectorData;


		/**
		* \brief Collection of hits produced by this particle in the scintillators 
		*/
		std::vector<TScintHit> fScintHitCollection;//all hits
		std::vector<TScintHit> fSelScintHitCollection;//selected hits (E>threshold, strip in veto, ...)
  
  
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

		
  	ClassDef(TParticleSimData,2)
		//ClassDef(TParticleSimData,1)
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


#endif
