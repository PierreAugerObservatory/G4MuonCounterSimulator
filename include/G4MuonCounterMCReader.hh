/**
* @file G4MuonCounterMCReader.hh
* @class G4MuonCounterMCReader
* @brief Read the output of the GEANT4 simulation
*
* Useful to read the output of the ROOT file produced in the simulation, get and draw sim info
* @author Dr. Simone Riggi, Dr.ssa Enrica Trovato
* @date 14/04/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4MuonCounterMCReader_h
#define _G4MuonCounterSimulatorUSC_G4MuonCounterMCReader_h 1

#include "TScintHit.hh"
#include "TPMTHit.hh"
#include "PMTSimulator.hh"
#include "TrackPoint.hh"

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

namespace G4MuonCounterSimulatorUSC {

class G4MuonCounterMCReader{

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    G4MuonCounterMCReader();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~G4MuonCounterMCReader();


		enum eventTag { kEMPTY = 0, kSINGLEPLANE = -1, kONEFOLD = 1,  kTWOFOLD= 2, kTHREEFOLD = 3};


	public:

		/** 
		\brief Set the input MC filename to be read and processed
 		*/
		void SetMCFileName(const char * fileName) { fMCFileName = fileName;}
		/** 
		\brief Set the output filename where to store processed info (i.e. histograms, graph, ...)
 		*/
		void SetOutputFileName(const char * fileName) { fOutputFileName = fileName;}


		/** 
		\brief Create data structures where to store processed information
 		*/
		void Init();//allocate (just once!) all histograms
		/** 
		\brief Reset allocated data structures (i.e. reset histograms, ...)
 		*/
 		void Reset();
		/** 
		\brief Read MC file and get needed information
 		*/
		void ReadInfo();
		/** 
		\brief Get needed information for detector design
 		*/
		void ReadDetectorInfo(int EventId);
		/** 
		\brief Get needed information for generated particles
 		*/
		void ReadGeneratedInfo(int EventId);
		/** 
		\brief Get needed information for scintillation hits
 		*/
		void ReadScintillatorInfo(int EventId);
		/** 
		\brief Get needed information for scintillation hits
 		*/
		void ReadPMTInfo(int EventId);

		/** 
		\brief Fill data structures (histograms) with information
 		*/
		void Fill();
		/** 
		\brief Draw data structures
 		*/
		void DrawInfo();
		/** 
		\brief Draw data structures
 		*/
		void StoreInfo();

		/** 
		\brief Set the threshold in deposited energy
 		*/
		void SetEnergyThreshold(double value){fEnergyThreshold= value;}
		/** 
		\brief Set the threshold in number of photoelectrons
 		*/
		void SetPEThreshold(double value){fPEThreshold= value;}
		/** 
		\brief Set the threshold in number of photons
 		*/
		void SetPhotonThreshold(double value){fPhotonThreshold= value;}
		/** 
		\brief Turn on/off the generation of PMT pulses
 		*/
		void GeneratePMTPulses(bool choice){fGeneratePMTPulses= choice;}
		/** 
		\brief Turn on/off the generation of PMT photons, according to a simple parametrization
 		*/
		void GeneratePMTPhotonsFast(bool choice){fGeneratePMTPhotonsFast= choice;}
		/** 
		\brief Turn on/off the generation of PMT pe, according to a simple parametrization
 		*/
		void GeneratePMTPEFast(bool choice){fGeneratePMTPEFast= choice;}
		/** 
		\brief Turn on/off the simple simulation (without photons tracking, only strip hits)
 		*/
		void IsFastSimulation(bool choice){fIsFastSimulation= choice;}
		/** 
		\brief Turn on/off the simple simulation (without photons tracking, only strip hits)
 		*/
		void SetFastSimulationSmearing(double choice){fFastSimulationSmearing= choice;}

		/** 
		\brief Set printing verbosity
 		*/
		void SetVerbosity(int choice){fVerbosity= choice;}

	private:
	
		/** 
		\brief Get quantum efficiency probability from QE table according to the energy entry E
 		*/
		double GetQEFromTable(double E);
		/** 
		\brief Get strip multiplicity
 		*/
		int GetStripMultiplicity(int AbsPlaneId);
		/** 
		\brief Create track points from scintillator hits
 		*/
		void CreateTrackPoints();
		

	private:
	
		TFile* fInputFile;
		TFile* fOutputFile;
		TTree* SimTree_pmt;
		TTree* SimTree_scint;
		TTree* SimTree_gen;
		TTree* SimTree_det;

		PMTSimulator* fPMTSim;
		std::vector<PMTSimulator*> fPMTSimCollection;

		/** 
		\brief Calculate geometrical hits in each planes
 		*/
		std::vector<int> CalculateGeomHits();
		/** 
		\brief Tag the events according to the number of involved detector planes
 		*/
		eventTag EventTagger(std::vector<int> hitPlaneList);


		const char * fMCFileName;
		const char * fOutputFileName;

		TDirectory* genDir;
		TDirectory* recDir;
		TDirectory* stripEdepDir;
		TDirectory* planeEdepDir;

		int fNevents;	
		int fVerbosity;

		//## Output histograms
		//## generation info
		TH1D* fGenEnergyHisto;
		TH1D* fGenMomentumXHisto;
		TH1D* fGenMomentumYHisto;
		TH1D* fGenMomentumZHisto;			
		TH1D* fGenVertexXHisto;
		TH1D* fGenVertexYHisto;
		TH1D* fGenVertexZHisto;
		TH1D* fGenThetaHisto;
		TH1D* fGenPhiHisto;

		//## rec info
		//scintillator
		TH1D* fEnergyDepositHisto;
		std::vector<TH1D*> fEnergyDepositHistoList;
		TH1D* fEnergyDepositInPlaneHisto;
		std::vector<TH1D*> fEnergyDepositInPlaneHistoList;
		TH1D* fStripMultiplicityInPlaneHisto;
		std::vector<TH1D*> fStripMultiplicityInPlaneHistoList;

		//pmt
		TH1D* fNumberOfEventsHisto;
		std::vector<TH1D*> fNumberOfEventsHistoList;
		TH1D* fPhotonCountsHisto;
		std::vector<TH1D*> fPhotonCountsHistoList;
		TH1D* fAveragePhotonCountsHisto;
		std::vector<TH1D*> fAveragePhotonCountsHistoList;		
		TH1D* fPhotonEfficiencyHisto;
		std::vector<TH1D*> fPhotonEfficiencyHistoList;
		
		TH1D* fPECountsHisto;
		std::vector<TH1D*> fPECountsHistoList;
		TH1D* fAveragePECountsHisto;
		std::vector<TH1D*> fAveragePECountsHistoList;
		TH1D* fPEEfficiencyHisto;
		std::vector<TH1D*> fPEEfficiencyHistoList;

		TH2D* fPhotonCountsMap;
		TH2D* fPECountsMap;
		TTree* fCountsMapTree;
		TTree* fTimeMapTree;
		double x;
		double t;
		double t_hit;
		double Edep;
		double d;
		int phCounts;
		int peCounts;

		//## GEN VARS 
		//vector loops on number of events in file and number of generated particles
		double X0;
		double Y0;
		double Z0;
		double Px;
		double Py;
		double Pz;
		double KinEnergy;
		double Mass;
    double Theta;
		double Phi;
		double Time;
		int PDGCode;

		//vector loops on number of events in file and number of planes
		std::vector< vector<int> > fGeomHitStripId;

		//## DET VARS 
		double Xstrip;
		double Ystrip;
		double Zstrip;	
		double Xcoat;
		double Ycoat;
		double Zcoat;	
		int Nstrips;
		int Nplanes;		
		static const int MAXNOPLANES=100;
		double XYPlaneDistance;
		double PlaneDistance[MAXNOPLANES];
		double PlaneDepth[MAXNOPLANES];
	

		//## SCINT VARS 
		std::vector<TScintHit>* fScintHitCollection;
		
		double fEnergyThreshold;

		struct HitStripStruct { 
			int superplaneID;
      int planeID;
   		int stripID;
			int absplaneID;
			int absstripID;
			double etot;
			double time;
			TVector3 stripPosition;
			TVector3 stripHitPosition;
			std::vector<TVector3> stripPMTPosition;
			std::vector<double> stripHitPMTDistance;
			std::vector<int> stripPMTPECounts;
			std::vector<int> stripPMTPhotonCounts;
		};
		
		struct order_by_superplaneID { 
    	bool operator()(HitStripStruct const &a, HitStripStruct const &b) { 
        return a.superplaneID < b.superplaneID;
    	}
		};
		struct order_by_absstripID { 
    	bool operator()(HitStripStruct const &a, HitStripStruct const &b) { 
        return a.absstripID < b.absstripID;
    	}
		};
		
		//vector looping over the events and number of hit strips
		std::vector<HitStripStruct> fHitStripList;
		std::vector<int> fSelectedStripList;
		std::vector<int> fSelectedPlaneList;

		std::vector<eventTag> fEventType;
		std::vector<TrackPoint> fTrackPointCollection;
		

		//## PMT VARS
		//vector looping over the events and number of hit strips
		std::vector<TPMTHit>* fPMTHitCollection;
		std::vector<int> fSelectedPMTStripList;
		std::vector<int> fSelectedPMTList;

		
		struct HitPMTStruct { 
			int superplaneID;
      int planeID;
   		int stripID;
			int pmtID;
			int absplaneID;
			int absstripID;
			int photonCounts;
			int peCounts;
			TVector3 pmtPosition;
			TVector3 pmtHitPosition;
			PMTSimulator* pmtSim;
		};
		std::vector<HitPMTStruct> fHitPMTList;

		
		double fPEThreshold;
		double fPhotonThreshold;
		bool fGeneratePMTPulses;
		bool fGeneratePMTPhotonsFast;
		bool fGeneratePMTPEFast;
		bool fIsFastSimulation;
		double fFastSimulationSmearing;
		
		std::map<double,double> QETable;

};

}//close namespace
#endif /*MCReader_h*/
