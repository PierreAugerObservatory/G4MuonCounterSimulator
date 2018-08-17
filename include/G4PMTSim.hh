/**
* @file G4PMTSim.hh
* @class G4PMTSim
* @brief Simulate the PMT signals (spe pulses, ...)
*
* Simulation of the PMT signals on the basis of the information produced in the G4 simulation
* @author Dr. Simone Riggi
* @date 20/08/2010
*/

#ifndef _G4MuonCounterSimulatorUSC_G4PMTSIM_h
#define _G4MuonCounterSimulatorUSC_G4PMTSIM_h 1

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TVector3.h"
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

using namespace std;

namespace G4MuonCounterSimulatorUSC {

class G4PMTSim {

	public:

		/** 
		\brief Class constructor: initialize structures.
 		*/
		G4PMTSim();
		/**
		* \brief Class destructor: free allocated memory
		*/
  	~G4PMTSim();

		void SetDataFileName(const char * fileName) { fDataFileName = fileName;}
		void SetOutputFileName(const char * fileName) { fOutputFileName = fileName;}
		
		/** 
		\brief Create data structures where to store processed information
 		*/
		void Init();
		/** 
		\brief Generate the photoelectrons at cathode
 		*/                                    
		bool GeneratePE();
		/** 
		\brief Generate the photons at cathode, according to a parametrization, given the hit position x and the energy deposit in strip
 		*/                                    
		bool GeneratePhotonsFast(double x,double Edep);
		/** 
		\brief Generate the pe at cathode, according to a parametrization, given the hit position x and the energy deposit in strip
 		*/                                    
		bool GeneratePEFast(double x,double Edep);
		/** 
		\brief Generate the PMT pulse 
 		*/
		void PulseFinder();
		/** 
		\brief Generate the PMT single photoelectron spectrum
 		*/
		void SPEPulseSpectrum(double time,int PulseModel);
		/** 
		\brief Generate the PMT noise pulse
 		*/
    void NoisePulseSpectrum(); 
		/** 
		\brief Execute the PMT simulation (signal generation)
 		*/                                    
		void DoPMTSim();

		/** 
		\brief PMT gain generator
 		*/
		double GainGenerator(int GainModel);
		/** 
		\brief Polya model for gain generation
 		*/
	  static double PolyaFcn(double* x,double* par);

		/** 
		\brief Available gain models
 		*/
		enum GainModel { kGauss=1, kPolya=2, kFixed=3};
		/** 
		\brief Available pulse models
 		*/
    enum PulseModel { kGaussTimePulse=1, kTruncGaussTimePulse=2, kLandauTimePulse=3};
		/** 
		\brief Available noise pulse models
 		*/
    enum NoisePulseModel { kGaussNoisePulse=1, kPoissonNoisePulse=2};


		//### PMT characteristics ###
		/** 
		\brief Set the PMT quantum efficiency 
 		*/
		void SetQE(double value) {fQE= value;}
		/** 
		\brief Get the PMT quantum efficiency
 		*/
		double GetQE() {return fQE;}
		/** 
		\brief Set the PMT quantum efficiency, according to a spectrum
 		*/
		void SetQETable(std::map<double,double> table) {fQETable= table;}
		/** 
		\brief Get the PMT quantum efficiency, according to a spectrum
 		*/
    std::map<double,double> GetQETable() {return fQETable;}

		//RiseTime
		/** 
		\brief Set the PMT signal rise time
 		*/
		void SetRiseTime(double value) {fRiseTime= value;}
		/** 
		\brief Get the PMT signal rise time
 		*/
		double GetRiseTime() {return fRiseTime;}

		//Transit time
		/** 
		\brief Set the PMT transit time (TTS)
 		*/
		void SetTransTime(double value) {fTransTime= value;}
		/** 
		\brief Get the PMT transit time (TTS)
 		*/
		double GetTransTime() {return fTransTime;}
		/** 
		\brief Set the PMT transit time spread (TTS)
 		*/
		void SetTransTimeSpread(double value) {fTransTimeSpread= value;}
		/** 
		\brief Get the PMT transit time spread (TTS)
 		*/
		double GetTransTimeSpread() {return fTransTimeSpread;}
		/** 
		\brief Set the PMT pulse model
 		*/
		void SetPulseModel(int value) {fPulseModel= value;}


    //Gain
		/** 
		\brief Set the PMT gain
 		*/
		void SetGain(double value) {fGain= value;}
		/** 
		\brief Get the PMT gain
 		*/
		double GetGain() {return fGain;}
		/** 
		\brief Set the PMT gain spread
 		*/
		void SetGainSpread(double value) {fGainSpread= value;}
		/** 
		\brief Get the PMT gain spread
 		*/
		double GetGainSpread() {return fGainSpread;}
		/** 
		\brief Set the PMT gain model
 		*/
		void SetGainModel(int value) {fGainModel= value;}
		/** 
		\brief Set the PMT noise baseline
 		*/
		void SetNoiseBaseline(double value) {fNoiseBaseline= value;}
		/** 
		\brief Get the PMT noise baseline
 		*/
		double GetNoiseBaseline() {return fNoiseBaseline;}
		/** 
		\brief Set the PMT noise RMS
 		*/
		void SetNoiseRMS(double value) {fNoiseRMS= value;}
		/** 
		\brief Get the PMT noise RMS
 		*/
		double GetNoiseRMS() {return fNoiseRMS;}
	
		/** 
		\brief Get the list of photon times at cathode
 		*/
		std::vector<double> GetPhotonTimeAtCathode() {return fPhotonTime;}
		/** 
		\brief Set the list of photon times at cathode
 		*/
		void SetPhotonTimeAtCathode(std::vector<double> vect) {fPhotonTime= vect;}	
		/** 
		\brief Get the list of photon energies at cathode
 		*/
		std::vector<double> GetPhotonEnergyAtCathode() {return fPhotonEnergy;}
		/** 
		\brief Set the list of photon energies at cathode
 		*/
		void SetPhotonEnergyAtCathode(std::vector<double> vect) {fPhotonEnergy= vect;}	
		/** 
		\brief Get the list of photon positions at cathode
 		*/
		std::vector<TVector3> GetPhotonPositionAtCathode() {return fPhotonPosition;}
		/** 
		\brief Set the list of photon positions at cathode
 		*/
		void SetPhotonPositionAtCathode(std::vector<TVector3> vect) {fPhotonPosition= vect;}	


		/** 
		\brief Get the list of photoelectrons times at cathode
 		*/
		std::vector<double> GetPETimeAtCathode() {return fPETime;}
		/** 
		\brief Set the list of photoelectrons times at cathode
 		*/
		void SetPETimeAtCathode(std::vector<double> vect) {fPETime= vect;}
		/** 
		\brief Get the list of photoelectrons energies at cathode
 		*/
		std::vector<double> GetPEEnergyAtCathode() {return fPEEnergy;}
		/** 
		\brief Set the list of photoelectrons energies at cathode
 		*/
		void SetPEEnergyAtCathode(std::vector<double> vect) {fPEEnergy= vect;}


		/** 
		\brief Set the PE counts
 		*/
		void SetPECounts(int value) {fPECounts= value;}
		/** 
		\brief Get the PE counts
 		*/
		int GetPECounts() {return fPECounts;}
		/** 
		\brief Set the Photon counts
 		*/
		void SetPhotonCounts(int value) {fPhotonCounts= value;}
		/** 
		\brief Get the Photon counts
 		*/
		int GetPhotonCounts() {return fPhotonCounts;}
		/** 
		\brief Set the parametrization smearing width
 		*/
		void SetParametrizationSmearingWidth(double value) {fParametrizationSmearingWidth= value;}
		/** 
		\brief Get the parametrization smearing width
 		*/
		double GetParametrizationSmearingWidth() {return fParametrizationSmearingWidth;}


		/** 
		\brief Set the QDC counts
 		*/
		void SetQDCCounts(double value) {fQDCCounts= value;}
		/** 
		\brief Get the QDC counts
 		*/
		double GetQDCCounts() {return fQDCCounts;}
	
		/** 
		\brief Set the arrival time
 		*/
		void SetArrivalTime(double value) {fArrivalTime= value;}
		/** 
		\brief Get the arrival time
 		*/
		double GetArrivalTime(double value) {return fArrivalTime;}

		/** 
		\brief Get the time trace histo
 		*/
		TH1D* GetTimeTrace(){return fTimeTraceHisto;}	
		/** 
		\brief Get the signal time trace histo
 		*/
		TH1D* GetSignalTimeTrace(){return fSignalTimeTraceHisto;}	
		/** 
		\brief Get the noise time trace histo
 		*/
		TH1D* GetNoiseTimeTrace(){return fNoiseTimeTraceHisto;}	
	


	private:


		//## FILE VARS
		const char* fDataFileName;
		TFile* fDataFile;
		TTree* SimTree;
		const char* fOutputFileName;
		TFile* fOutputFile;
		
		
		static TH1D* fTimeTraceHisto;
		static TH1D* fSignalTimeTraceHisto;
		static TH1D* fNoiseTimeTraceHisto;
		static TH1D* fPETimeDistributionHisto;
		static TH1D* fSPETimeTraceHisto;
		static std::vector<TH1D*> fSPETimeTrace;	

		
		int fPulseModel;
		int fGainModel;
	
		double fQE;//integrated QE
		std::map<double,double> fQETable;//spectrum of QE 

		double fFWHM;
		double fRiseTime;
		double fTransTime;
		double fTransTimeSpread;
		double fGain;
		double fGainSpread;
		double fDarkCurrent;
		double fDarkCurrentMax;
		double fDarkCurrentSpread;
		double fNoiseBaseline;
		double fNoiseRMS;

		double fQDCCounts;
		double fQDCPedestal;
		double fCharge2QDCCounts;

		double fDiscriminatorThreshold;
		double fArrivalTime;

		int fPECounts;
		int fPhotonCounts;
		double fParametrizationSmearingWidth;

		double fPMTTimeTraceBinning;
	  double fPMTTimeTraceMin;
		double fPMTTimeTraceMax;
		int fNTimeBins;
	
		static const double ElectronCharge= 1.60217646e-19;

		std::vector<double> fPhotonTime; 
		std::vector<double> fPhotonEnergy;
		std::vector<TVector3> fPhotonPosition;
  	std::vector<double> fPETime;
		std::vector<double> fPEEnergy;

		std::vector<double> fSPETrace;
		std::vector<double> fSignalTrace;
		std::vector<double> fNoiseTrace;

		std::vector<double> fSPETraceContainer;
		std::vector<double> fNoiseTraceContainer;
		std::vector<double> fTimeBinCenter;

		/** 
		\brief Get quantum efficiency probability from QE table according to the energy entry E
 		*/
		double GetQEFromTable(double E);

};

}//close namespace

#endif
