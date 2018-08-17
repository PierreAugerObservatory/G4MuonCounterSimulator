/**
* @file TrackReco.hh
* @class TrackReco
* @brief Reconstruct the particle tracks in the detector
*
* Reconstruct the particle tracks in the detector
* @author Dr. Simone Riggi
* @date 07/09/2010
*/

#ifndef TrackReco_h
#define TrackReco_h 1

#include "TrackPoint.hh"
#include "Track.hh"

#include <vector>
#include <map>
#include <string>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TGraph.h>
#include <TVector3.h>

using namespace std;

class TrackReco{

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    TrackReco();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~TrackReco();

		/**
		* \brief Reconstruct tracks from collection of hit points
		*/
   	bool TrackReconstruction();


	public:

		/**
		* \brief Set the collection of track points
		*/
   	void SetTrackPoints(std::vector<TrackPoint> vect){fTrackPointCollection= vect; }
		/**
		* \brief Set the collection of track points
		*/
   	std::vector<TrackPoint> GetTrackPoints(){return fTrackPointCollection;}

		/**
		* \brief Get the collection of reconstructed tracks
		*/
   	std::vector<Track*> GetRecTrackCollection(){return fRecTrackCollection;}
		/**
		* \brief Get the collection of cluster points
		*/
   	std::vector<TrackPoint> GetClusterTrackPoints(){return fClusterTrackPointCollection;}

		/**
		* \brief Set the collection of track points
		*/
   	void SetHitSeparationInCluster(double value){fHitSeparationInCluster= value; }
		/**
		* \brief Set the collection of track points
		*/
   	double GetHitSeparationInCluster(){return fHitSeparationInCluster;}
		/**
		* \brief Set the track point default X uncertainty
		*/
   	void SetHitSigmaX(double value){fHitSigmaX= value; }
		/**
		* \brief Set the track point default X uncertainty
		*/
   	double GetHitSigmaX(){return fHitSigmaX;}
		/**
		* \brief Set the track point default Y uncertainty
		*/
   	void SetHitSigmaY(double value){fHitSigmaY= value; }
		/**
		* \brief Set the track point default Y uncertainty
		*/
   	double GetHitSigmaY(){return fHitSigmaY;}
		/**
		* \brief Set the track point default Z uncertainty
		*/
   	void SetHitSigmaZ(double value){fHitSigmaZ= value; }
		/**
		* \brief Set the track point default Z uncertainty
		*/
   	double GetHitSigmaZ(){return fHitSigmaZ;}

		/**
		* \brief Set the superplane size Z
		*/
   	void SetSuperPlaneSizeZ(double value){fSuperPlaneSizeZ= value; }
		/**
		* \brief Set the superplane size Z
		*/
   	double GetSuperPlaneSizeZ(){return fSuperPlaneSizeZ;}

		/**
		* \brief Set the superplane Z depths
		*/
   	void SetSuperPlaneDepth(std::vector<double> value){fSuperPlaneDepth= value; }
		/**
		* \brief Set the superplane Z depths
		*/
   	std::vector<double> GetSuperPlaneDepth(){return fSuperPlaneDepth;}

	private:

		
		int fVerbosity;

		std::vector<TrackPoint> fTrackPointCollection;
		std::vector<TrackPoint> fClusterTrackPointCollection;
		std::vector<Track*> fRecTrackCollection;
		
		double fHitSigmaX;
		double fHitSigmaY;
		double fHitSigmaZ;
		double fHitSeparationInCluster;
		std::vector<double> fSuperPlaneDepth;
		double fSuperPlaneSizeZ;

};

#endif /*TrackReco_h*/
