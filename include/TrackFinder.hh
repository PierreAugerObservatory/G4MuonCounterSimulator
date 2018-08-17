/**
* @file TrackFinder.hh
* @class TrackFinder
* @brief Reconstruct the particle tracks in the detector
*
* Reconstruct the particle tracks in the detector
* @author S. Riggi
* @date 07/09/2010
*/

#ifndef TrackFinder_h
#define TrackFinder_h 1

#include "TrackPoint.hh"
#include "Track.hh"
#include "Cluster.hh"

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

class TrackFinder{

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    TrackFinder();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~TrackFinder();

		/**
		* \brief Reconstruct tracks from collection of hit points
		*/
   	bool TrackReconstructor();
		/**
		* \brief Find tracks from collection of hit points
		*/
   	bool FindTrack();
		/**
		* \brief Find clusters in each plane starting from the given hits
		*/
		bool FindCluster();
 

		/**
		* \brief Set verbosity
		*/
   	void SetVerbosity(int value){fVerbosity=value;}


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

		/**
		* \brief Set the superplane trigger times
		*/
   	void SetTriggerTime(std::vector<double> value){fTriggerTime= value; }
		/**
		* \brief Set the superplane Z depths
		*/
   	std::vector<double> GetTriggerTime(){return fTriggerTime;}

		/**
		* \brief Set the true energy
		*/
   	void SetGenEnergy(double value){fGenEnergy= value; }
		/**
		* \brief Get the true energy
		*/
   	double GetGenEnergy(){return fGenEnergy;}	
		/**
		* \brief Set the true theta
		*/
   	void SetGenTheta(double value){fGenTheta= value; }
		/**
		* \brief Get the true theta
		*/
   	double GetGenTheta(){return fGenTheta;}	
		/**
		* \brief Set the true phi
		*/
   	void SetGenPhi(double value){fGenPhi= value; }
		/**
		* \brief Set the true phi
		*/
   	double GetGenPhi(){return fGenPhi;}	
		/**
		* \brief Set the true tx
		*/
   	void SetGenTx(double value){fGenTx= value; }
		/**
		* \brief Get the true tx
		*/
   	double GetGenTx(){return fGenTx;}	
		/**
		* \brief Set the true ty
		*/
   	void SetGenTy(double value){fGenTy= value; }
		/**
		* \brief Get the true ty
		*/
   	double GetGenTy(){return fGenTy;}		
		/**
		* \brief Set the true vertex X
		*/
   	void SetGenVertexX(double value){fGenVertexX= value; }
		/**
		* \brief Get the true vertex X
		*/
   	double GetGenVertexX(){return fGenVertexX;}			
		/**
		* \brief Set the true vertex Y
		*/
   	void SetGenVertexY(double value){fGenVertexY= value; }
		/**
		* \brief Get the true vertex Y
		*/
   	double GetGenVertexY(){return fGenVertexY;}	
		/**
		* \brief Set the run number
		*/
   	void SetRunNumber(int value){fRunNumber= value;}	


			
		/**
		* \brief Turn on/off the fitting with MS correction
		*/
   	void TrackWithMSCorrection(bool choice){fTrackWithMSCorrection= choice;}	
		/**
		* \brief Turn on/off the energy loss correction in MS treatment
		*/
   	void UseEnergyLossCorrectionInMS(bool choice){fUseEnergyLossCorrectionInMS= choice;}	
		/**
		* \brief Use an average energy in MS treatment 
		*/
   	void UseAverageEnergyInMS(bool choice){fUseAverageEnergyInMS= choice;}
		/**
		* \brief Use an average energy in MS treatment 
		*/
   	void SetAverageEnergyInMS(double choice){fAverageEnergyInMS= choice;}	
		/**
		* \brief Use the true energy in MS treatment 
		*/
   	void UseTrueEnergyInMS(bool choice){fUseTrueEnergyInMS= choice;}	
		/**
		* \brief Use the true energy in MS treatment 
		*/
   	void IncludeMomentumInTracking(bool choice){fIncludeMomentumInTracking= choice;}			

		/**
		* \brief Use clusters in tracking and not single hits
		*/
   	void UseClusterInTracking(bool choice){fUseClusterInTracking= choice;}	
		/**
		* \brief Split clusters according to hits multiplicity threshold
		*/
   	void SplitClustersInTracking(bool choice){fSplitClustersInTracking= choice;}
		/**
		* \brief Split clusters according to hits multiplicity threshold
		*/
   	void SetSplitClusterThreshold(int value){fSplitClusterThreshold= value;}
		/**
		* \brief Turns on/off fake hits in tracking
		*/
   	void RemoveFakeHitsInTracking(bool choice){fRemoveFakeHitsInTracking= choice;}
		/**
		* \brief Remove fake hits in tracking
		*/
   	void RemoveFakeHits();

	private:
	
		/**
		* \brief Define straight line equation in 3D in parametric form
		*/
		void GeomLine3D(double t, double *p, double &x, double &y, double &z);
		/**
		* \brief Define distance in 3D between a point and a straight line
		*/
		static double GeomDistance3D(double x,double y,double z, double* p);
		
		/**
		* \brief Fit a candidate track with a straight 3D line
		*/		
		void TrackFitter(Track* aTrack);
		/**
		* \brief Fit a candidate track with a straight 3D line
		*/		
		void TrackFitterWithMS(Track* aTrack);
		/**
		* \brief Fit a candidate track with a straight 3D line - analytic solution
		*/		
		int FastTrackFitter(Track* aTrack);
		/**
		* \brief Fit a candidate track with a straight 3D line - analytic solution + MultipleScattering
		*/		
		int FastTrackFitterWithMS(Track* aTrack);
		/**
		* \brief Fit function definition for track fitter method
		*/		
		static void TrackChi2Fcn(int& nPar, double* const grad,double& value, double* par,const int iFlag);
		/**
		* \brief Fit function definition for track fitter method
		*/		
		static void TrackChi2WithMSFcn(int& nPar, double* const grad,double& value, double* par,const int iFlag);


		/**
		* \brief Get the number of all possible combinations obtained with the track points
		*/
		int GetCombinationNo();
		/**
		* \brief Calculate the combination matrix
		*/
		void CalculateCombinationMatrix();
		/**
		* \brief Check if tracks can be reconstructed
		*/
		bool Checker();

		/**
		* \brief Set track point uncertainty
		*/
		void SetTrackPointUncertainty();
		/**
		* \brief Calculate elements of covariance matrix due to multiple scattering
		*/
		static double MSCovariance(double Zi, double Zj, double Zs, double ZsStart, double L, double H);

	private:


		int fRunNumber;
		int fVerbosity;

		static bool fIncludeMomentumInTracking;
		static bool fUseAverageEnergyInMS;
		static double fAverageEnergyInMS;
		static bool fUseTrueEnergyInMS;
		static bool fUseEnergyLossCorrectionInMS;
		static bool fTrackWithMSCorrection;
		bool fUseClusterInTracking;
		bool fSplitClustersInTracking;	
		int fSplitClusterThreshold;
		bool fRemoveFakeHitsInTracking;

		static int fFitStatus;
		static double fFitFCNMin;

		std::vector<TrackPoint> fTrackPointCollection;
		std::vector<TrackPoint> fCandidateTrackPointCollection;
		std::vector<TrackPoint> fClusterTrackPointCollection;
		std::vector< std::vector<TrackPoint> > fTrackPointCollectionInPlane;
			
		std::vector< std::vector<int> > fCombinationMatrix;	
		std::vector<Track*> fRecTrackCollection;

		double fHitSigmaX;
		double fHitSigmaY;
		double fHitSigmaZ;
		double fHitSeparationInCluster;
		static std::vector<double> fSuperPlaneDepth;
		static double fSuperPlaneSizeZ;
		std::vector<double> fTriggerTime;

		static double fTheta;
		static double fPhi;
		static double fVertexX;
		static double fVertexY;
		static double fTx;
		static double fTy;

		static double fGenEnergy;
		static double fGenTheta;
		static double fGenPhi;
		static double fGenVertexX;
		static double fGenVertexY;
		static double fGenTx;
		static double fGenTy;
		static double fThetaStart;
		static double fPhiStart;
		static double fVertexXStart;
		static double fVertexYStart;
		static double fTxStart;
		static double fTyStart;

		static double fThetaErr;
		static double fPhiErr;
		static double fVertexXErr;
		static double fVertexYErr;
		static double fTxErr;
		static double fTyErr;

		static double fThetaStartErr;
		static double fPhiStartErr;
		static double fVertexXStartErr;
		static double fVertexYStartErr;
		static double fTxStartErr;
		static double fTyStartErr;

		static const int MAXTRACKNO= 100;
		int nTrack;
		double RecTheta[MAXTRACKNO];
		double RecPhi[MAXTRACKNO];
		double RecX[MAXTRACKNO];
		double RecY[MAXTRACKNO];
		int TrackFitStatus[MAXTRACKNO];
		double TrackFitChi2[MAXTRACKNO];

		static const int MAXHITNO= 5000; 
		static double lenz[MAXHITNO]; 
		int nHit;	
		double HitX[MAXHITNO];
		double HitY[MAXHITNO];
		double HitZ[MAXHITNO];
		double HitEdepX[MAXHITNO];
		double HitEdepY[MAXHITNO];
		double HitTimeX[MAXHITNO];
		double HitTimeY[MAXHITNO];
		int IsMuonHit[MAXHITNO];
		int nClusterHit;
		double ClusterHitX[MAXHITNO];
		double ClusterHitY[MAXHITNO];
		double ClusterHitZ[MAXHITNO];
	
		int nTrackHit;
		double TrackHitX[MAXHITNO];
		double TrackHitY[MAXHITNO];
		double TrackHitZ[MAXHITNO];
		double TrackHitEdepX[MAXHITNO];
		double TrackHitEdepY[MAXHITNO];
		double TrackHitTimeX[MAXHITNO];
		double TrackHitTimeY[MAXHITNO];
		int TrackId[MAXHITNO];
	
		static int nHits;
		static int nPar;

		static double fKinEnergy;
		static double fMomentum;

		static const int MAXSCATNO= 10;
		static double Ls[MAXSCATNO];
		static double Zs[MAXSCATNO];
		static double ZsStart[MAXSCATNO];

};

#endif /*TrackFinder_h*/
