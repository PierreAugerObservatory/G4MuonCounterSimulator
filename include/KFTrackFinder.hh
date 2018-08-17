/**
* @file KFTrackFinder.hh
* @class KFTrackFinder
* @brief Reconstruct the particle tracks in the detector
*
* Reconstruct the particle tracks in the detector
* @author S. Riggi
* @date 07/09/2010
*/

#ifndef KFTrackFinder_h
#define KFTrackFinder_h 1

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
#include <TCanvas.h>

using namespace std;

class KFTrackFinder{

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    KFTrackFinder();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~KFTrackFinder();

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
		* \brief Get the collection of track points
		*/
   	double GetHitSeparationInCluster(){return fHitSeparationInCluster;}
		/**
		* \brief Set the track point default X uncertainty
		*/
   	void SetHitSigmaX(double value){fHitSigmaX= value; }
		/**
		* \brief Get the track point default X uncertainty
		*/
   	double GetHitSigmaX(){return fHitSigmaX;}
		/**
		* \brief Set the track point default Y uncertainty
		*/
   	void SetHitSigmaY(double value){fHitSigmaY= value; }
		/**
		* \brief Get the track point default Y uncertainty
		*/
   	double GetHitSigmaY(){return fHitSigmaY;}
		/**
		* \brief Set the track point default Z uncertainty
		*/
   	void SetHitSigmaZ(double value){fHitSigmaZ= value; }
		/**
		* \brief Get the track point default Z uncertainty
		*/
   	double GetHitSigmaZ(){return fHitSigmaZ;}
		/**
		* \brief Set the superplane size Z
		*/
   	void SetSuperPlaneSizeZ(double value){fSuperPlaneSizeZ= value; }
		/**
		* \brief Get the superplane size Z
		*/
   	double GetSuperPlaneSizeZ(){return fSuperPlaneSizeZ;}

		/**
		* \brief Set the superplane Z depths
		*/
   	void SetSuperPlaneDepth(std::vector<double> value){fSuperPlaneDepth= value; }
		/**
		* \brief Get the superplane Z depths
		*/
   	std::vector<double> GetSuperPlaneDepth(){return fSuperPlaneDepth;}

		/**
		* \brief Set the superplane trigger times
		*/
   	void SetTriggerTime(std::vector<double> value){fTriggerTime= value; }
		/**
		* \brief Get the superplane Z depths
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
		* \brief Set the hits multiplicity threshold to split clusters
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
		* \brief Fit a candidate track with a Kalman Filter 
		*/
		int KFTrackFitter(Track* aTrack);
		/**
		* \brief Fit a candidate track with a fast track fitter
		*/
		static int LSTrackFitter(Track* aTrack);
	
		/**
		* \brief Fit function definition for track fitter method
		*/		
		static void TrackChi2Fcn(int& nPar, double* const grad,double& value, double* par,const int iFlag);
		/**
		* \brief Fit a candidate track with a straight 3D line - analytic solution + MultipleScattering
		*/		
		int KFTrackFitterWithMomentum(Track* aTrack);
		

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
		* \brief Accept candidate track?
		*/
		bool AcceptCandidateTrack(Track* aCandidateTrack);	
		/**
		* \brief Accept hit combination?
		*/
		bool AcceptHitCombination(int* aCandidateCombination);

		/**
		* \brief Debug printing
		*/
		void DebugPrint();
		/**
		* \brief Kalman Filter Prediction Step
		*/
		static void KFPredictionStep(int k);
		/**
		* \brief Kalman Filter filtering Step
		*/
		static void KFFilterStep(int k);
		/**
		* \brief Kalman Filter smoothing Step
		*/
		static void KFSmoothingStep(int k);
		/**
		* \brief Calculate elements of covariance matrix due to multiple scattering
		*/
		static void MSCovariance(int k);

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
		double fMaxAllowedAngle;

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

		static TVector3 fDirection;
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
		static double TrackHitChi2[MAXHITNO];
		static double TrackHitPredictedChi2[MAXHITNO];
	
		int TrackId[MAXHITNO];

		static int nHits;
		static int nPar;
		static TMatrixD fKFInitState;
		static TMatrixD fKFInitCovariance;
		static std::vector<TMatrixD> fKFState;
		static std::vector<TMatrixD> fKFCovariance;
		static std::vector<TMatrixD> fKFPredictedState;
		static std::vector<TMatrixD> fKFPredictedCovariance;
		static std::vector<TMatrixD> fKFMeas;	
		static std::vector<TMatrixD> fKFMeasVariance;
		static std::vector<TMatrixD> fKFPropagator;
		static std::vector<TMatrixD> fKFProjector;
		static std::vector<TMatrixD> fKFCovarianceMS;
		static std::vector<TMatrixD> fKFMeasResidual;
		static std::vector<TMatrixD> fKFCovarianceResidual;
		static std::vector<TMatrixD> fKFPredictedMeasResidual;
		static std::vector<TMatrixD> fKFPredictedCovarianceResidual;
		static std::vector<TMatrixD> fKFSmoothedState;
		static std::vector<TMatrixD> fKFSmoothedCovariance;
		static std::vector<TMatrixD> fKFSmoothedMeasResidual;
		static std::vector<TMatrixD> fKFSmoothedCovarianceResidual;

		static std::vector<double> fKFHitChi2;
		static std::vector<double> fKFHitPredictedChi2;
		static double fKinEnergy;
		static double fMomentum;

		static const int MAXSCATNO= 10;
		static double Ls[MAXSCATNO];
		static double Zs[MAXSCATNO];
		static double ZsStart[MAXSCATNO];


};

#endif /*KFTrackFinder_h*/
