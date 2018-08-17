
/**
* @file G4MuonCounterReconstructor.hh
* @class module
* @brief module used to handle the reconstruction of simulated data
*
* @author S. Riggi
* @date 05/04/2010
*/

#ifndef _G4MuonCounterReconstructorUSC_G4MuonCounterReconstructor_h_
#define _G4MuonCounterReconstructorUSC_G4MuonCounterReconstructor_h_

#include <fwk/VModule.h>

#include <G4MuonCounterSimulator.hh>
#include <TScintHit.hh>
#include <TPMTHit.hh>
#include <PMTSimulator.hh>
#include <TrackPoint.hh>
#include <TStationSimData.hh>
#include <TEventSimData.hh>
#include <MuonDetector.hh>


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
#include <TGraphErrors.h>
#include <TVector3.h>
#include <TFile.h>
#include <TTree.h>

//using namespace G4MuonCounterSimulatorUSC;
namespace G4MuonCounterSimulatorUSC{
	class G4MuonCounterSimulator;
	class PMTSimulator;
}


namespace G4MuonCounterReconstructorUSC {
	
	//class G4MuonCounterSimulator;

  class G4MuonCounterReconstructor : public fwk::VModule {

  public:
    G4MuonCounterReconstructor();
    virtual ~G4MuonCounterReconstructor();

    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& theEvent);
		fwk::VModule::ResultFlag RunFromSimulator(evt::Event& theEvent);
		fwk::VModule::ResultFlag RunFromFile();
    fwk::VModule::ResultFlag Finish();

    std::string GetSVNId() const
    { return std::string("$Id$"); }


		enum TrackAlgo { eLeastSquare = 1, eKalmanFilter = 2};


		/** 
		\brief Set the output filename where to store processed info (i.e. histograms, graph, ...)
 		*/
		void SetInputFileName(const char * fileName) { fInputFileName = fileName;}
		/** 
		\brief Set the output filename where to store processed info (i.e. histograms, graph, ...)
 		*/
		void SetOutputFileName(const char * fileName) { fOutputFileName = fileName;}

		/** 
		\brief Create data structures where to store processed information
 		*/
		void InitInfoFromFile();//allocate (just once!) all histograms
		

		/** 
		\brief Get needed information for detector design
 		*/
		void FillDetectorInfo();
		/** 
		\brief Get needed information for sim particle
 		*/
		//fwk::VModule::ResultFlag FillSimParticleInfo();	
		void FillSimParticleInfo();	
		/** 
		\brief Get needed information for sim station
 		*/
		void FillSimStationInfo();



		/** 
		\brief Set the threshold in deposited energy
 		*/
		void SetEnergyThreshold(double value){fEnergyThreshold= value;}

		/** 
		\brief Set printing verbosity
 		*/
		void SetVerbosity(int choice){fVerbosity= choice;}
		/** 
		\brief Process particle hits 
 		*/
		void ProcessParticleHits();
		/** 
		\brief Init station hits 
 		*/
		void InitStationHits();
		/** 
		\brief Merge particle hits 
 		*/
		void MergeParticleHits();
		/** 
		\brief Process station hits 
 		*/
		void ProcessStationHits();
		/** 
		\brief Find station start time
 		*/
		void FindStationStartTime();

		/** 
		\brief Calculate max relative angle among hits for a given particle track 
 		*/
		double CalculateMaxRelAngle(TParticleSimData* aSimParticle);

		/** 
		\brief Check if a particle is inside the geometrical FOV
 		*/
		int GeomEventTagger(TParticleSimData* aSimParticle);

		/** 
		\brief Identify events according to NN response
 		*/
		bool IsMuonEvent(std::vector<double> NNInputPars);


  private:
    TFile* fInputFile;
		std::string fInputFileName;
		TTree* DetTree;
		TTree* SimTree;
		TTree* GenTree;
		
		TFile* fOutputFile;
		std::string fOutputFileName;
		TTree* RecTree;
		TTree* TrackTree;

		bool fSaveRecInfo;
		bool fReadSimInfoFromFile;
		bool fProcessHits;
		bool fSaveRecTree;
		bool fSaveTrackTree;
		bool fSaveHisto;

		bool fProcessSelectedHits;
		double fEnergyThreshold;	
		double fPEThreshold;
		double fPhotonThreshold;

		bool fMergeSimParticles;
		double fMergeTimeInterval;
		double fSamplingTime;
		double fCurrentSamplingTime;
		int fNSamplingBins;	
		

		bool fTrackMuons;
		int fTrackingAlgo;
		bool fTrackWithMSCorrection;
		bool fUseAverageEnergyInMS;
		double fAverageEnergyInMS;
		bool fUseTrueEnergyInMS;
		bool fUseEnergyLossCorrectionInMS;
		bool fIncludeMomentumInTracking;
		
		bool fUseClusterInTracking;
		bool fSplitClustersInTracking;	
		int fSplitClusterThreshold;
		bool fRemoveFakeHitsInTracking;
		bool fUsePID;
		double fNNCutInPID;

		G4MuonCounterSimulatorUSC::G4MuonCounterSimulator* fG4MuonCounterSimulator;
		MuonDetector* fDetectorData;
		std::vector<TStationSimData> fStationSimDataCollection;
		TEventSimData* fEventSimData;
		TStationSimData* currentStationSimData;
		TParticleSimData* currentSimParticle;
		TScintHit* currentScintHit;
		
		G4MuonCounterSimulatorUSC::PMTSimulator* fPMTSim;
		std::vector<G4MuonCounterSimulatorUSC::PMTSimulator*> fPMTSimCollection;

	
		int fNevents;	
		int fVerbosity;
		int fCurrentEventNo;
		int fCurrentStation;


		//## Output histograms
		
		static const int nMaxPlanes= 10;
		TH1D* fEnergyDepositInPlaneHisto[nMaxPlanes];
		TH1D* fElectronEnergyDepositInPlaneHisto[nMaxPlanes];
		TH1D* fMuonEnergyDepositInPlaneHisto[nMaxPlanes];
		TH1D* fGammaEnergyDepositInPlaneHisto[nMaxPlanes];	
		TH1D* fProtonEnergyDepositInPlaneHisto[nMaxPlanes];
		TH1D* fStripMultiplicityInPlaneHisto[nMaxPlanes];
		TGraphErrors* fMuonDensityInPlaneGraph[nMaxPlanes];
		TGraphErrors* fMuonDensityAlbedoInPlaneGraph[nMaxPlanes];
		int nMuonDensityPointsCounter[nMaxPlanes];
		int nMuonDensityAlbedoPointsCounter[nMaxPlanes];
		TGraphErrors* fEmDensityInPlaneGraph[nMaxPlanes];	
		TGraphErrors* fEmDensityAlbedoInPlaneGraph[nMaxPlanes];
		int nEmDensityPointsCounter[nMaxPlanes];
		int nEmDensityAlbedoPointsCounter[nMaxPlanes];

		TH1D* fTrigger3FoldEnergyHisto;
		TH1D* fTrigger2FoldEnergyHisto;

		//## GEN VARS 
		double lgE;
		double GPSTimeSec, GPSTimeNanoSec;
		double Theta, Phi;
		double Tx, Ty;
		double Rp, Psi;
		double CoreX, CoreY;
		double nMu, nMuPlus, nMuMinus;
		double nProton, nNeutron;
		double nE, nEPlus, nEMinus, nGamma;
		int HadronicModel, PrimaryParticle;
		double rCut;
		int ShowerNumber;
		const char* ShowerRunId;

		//## DET VARS
		int nStrips;
		int nPlanes;
		double StripSeparation;

		//## REC VARS
		int Id;
		double E;
		double ThetaRec;
		double PhiRec;
		double Time;
		double X;
		double Y;
		int StationId;
		double RStat;
		double RStatGrd;
		double R;
		double RGrd;
		double PsiGrd;
		int EventTag;
		int GeomEventTag;
		int HasGeomCrossedFirstPlane;
		static const int nMaxHitPlanes= 10;
		int nHitPlanes;
		double AverageXHit[nMaxHitPlanes];
		double AverageYHit[nMaxHitPlanes];
		double MultiplicityX[nMaxHitPlanes];
		double MultiplicityY[nMaxHitPlanes];
		double RX[nMaxHitPlanes];
		double RY[nMaxHitPlanes];
		double TriggerTimeX[nMaxHitPlanes];
		double TriggerTimeY[nMaxHitPlanes];

			
		
		double MuonTrackTxAtPlaneX[nMaxHitPlanes];
		double MuonTrackTxAtPlaneY[nMaxHitPlanes];
		double MuonTrackTyAtPlaneX[nMaxHitPlanes];
		double MuonTrackTyAtPlaneY[nMaxHitPlanes];

		double nMuonsX[nMaxHitPlanes];
		double nMuonsAlbedoX[nMaxHitPlanes];
		double nMuonsY[nMaxHitPlanes];
		double nMuonsAlbedoY[nMaxHitPlanes];
		double nEmX[nMaxHitPlanes];
		double nEmAlbedoX[nMaxHitPlanes];
		double nEmY[nMaxHitPlanes];
		double nEmAlbedoY[nMaxHitPlanes];
		double nHadronsX[nMaxHitPlanes];
		double nHadronsAlbedoX[nMaxHitPlanes];
		double nHadronsY[nMaxHitPlanes];
		double nHadronsAlbedoY[nMaxHitPlanes];

		double fHitSeparationInCluster;
		double fHitSigmaX;
		double fHitSigmaY;
		double fHitSigmaZ;
		double fDetectorAreaSize;
		std::vector<double> fSuperPlaneDepth;
		double fSuperPlaneSizeZ;
	
		

		//## Muon Tracking Info
		static const int MAXTRACKNO= 10000;
		int nTrack;
		int nTrueTrack;
		int nTrueMuonTrack;
		
		double GenE[MAXTRACKNO];
		double GenTheta[MAXTRACKNO];
		double GenPhi[MAXTRACKNO];
		double GenX[MAXTRACKNO];
		double GenY[MAXTRACKNO];
		double GenTx[MAXTRACKNO];
		double GenTy[MAXTRACKNO];
		double GenR[MAXTRACKNO];
		double GenMaxRelAngle[MAXTRACKNO];

		double RecOpeningAngle[MAXTRACKNO];
		double RecTheta[MAXTRACKNO];
		double RecPhi[MAXTRACKNO];
		double RecX[MAXTRACKNO];
		double RecY[MAXTRACKNO];
		double RecTx[MAXTRACKNO];
		double RecTy[MAXTRACKNO];
		int TrackFitStatus[MAXTRACKNO];
		int TrackRecoStatus;
		double TrackFitChi2[MAXTRACKNO];
		double RecThetaStart[MAXTRACKNO];
		double RecPhiStart[MAXTRACKNO];
		double RecXStart[MAXTRACKNO];
		double RecYStart[MAXTRACKNO];
		double RecTxStart[MAXTRACKNO];
		double RecTyStart[MAXTRACKNO];

		double RecThetaErr[MAXTRACKNO];
		double RecPhiErr[MAXTRACKNO];
		double RecXErr[MAXTRACKNO];
		double RecYErr[MAXTRACKNO];
		double RecTxErr[MAXTRACKNO];
		double RecTyErr[MAXTRACKNO];
		double RecThetaStartErr[MAXTRACKNO];
		double RecPhiStartErr[MAXTRACKNO];
		double RecXStartErr[MAXTRACKNO];
		double RecYStartErr[MAXTRACKNO];
		double RecTxStartErr[MAXTRACKNO];
		double RecTyStartErr[MAXTRACKNO];

		static const int MAXHITNO= 10000; 
		int nHit;
		double HitX[MAXHITNO];
		double HitY[MAXHITNO];
		double HitZ[MAXHITNO];
		double HitTimeX[MAXHITNO];
		double HitEdepX[MAXHITNO];
		double HitTimeY[MAXHITNO];
		double HitEdepY[MAXHITNO];
		int IsMuonHit[MAXHITNO];
		int HitPlaneId[MAXHITNO];
		int nClusterHit;
		double ClusterHitX[MAXHITNO];
		double ClusterHitY[MAXHITNO];
		double ClusterHitZ[MAXHITNO];
		int nTrackHit;
		double TrackHitX[MAXHITNO];
		double TrackHitY[MAXHITNO];
		double TrackHitZ[MAXHITNO];
		double TrackHitTimeX[MAXHITNO];
		double TrackHitEdepX[MAXHITNO];
		double TrackHitTimeY[MAXHITNO];
		double TrackHitEdepY[MAXHITNO];
		double TrackHitChi2[MAXHITNO];
		double TrackHitExpChi2[MAXHITNO];
		int TrackId[MAXHITNO];
		int nMuonTracks;
		int IsMuon;
		
		double fStartTime;


    REGISTER_MODULE("G4MuonCounterReconstructorUSC", G4MuonCounterReconstructor);

  };

}//close namespace

#endif
