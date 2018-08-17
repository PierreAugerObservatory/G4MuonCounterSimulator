/**
* @file MuonDetector.hh
* @class MuonDetector
* @brief Container for muon detector info
*
* @author S. Riggi
* @date 31/12/2010
*/

#ifndef MuonDetector_h
#define MuonDetector_h 1

#include <Rtypes.h>
#include <TObject.h>

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

class MuonDetector: public TObject{

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MuonDetector();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~MuonDetector();


	public:
			
		/**
		* \brief Get the number of stations
		*/
		int GetNumberOfStations() {return fNStations;}
		/**
		* \brief Set the number of stations
		*/
		void SetNumberOfStations(int value) {fNStations= value;}
		/**
		* \brief Add a station to the total number of stations
		*/
		void AddStation() {fNStations++;}
		/**
		* \brief Get the number of strips per plane
		*/
		int GetNumberOfStrips() {return fNStrips;}
		/**
		* \brief Set the number of strips per plane
		*/
		void SetNumberOfStrips(int value) {fNStrips= value;}
		/**
		* \brief Get the number of planes 
		*/
		int GetNumberOfPlanes() {return fNPlanes;}
		/**
		* \brief Set the number of planes
		*/
		void SetNumberOfPlanes(int value) {fNPlanes= value;}

		
		/**
		* \brief Get the superplane thickness
		*/
		double GetSuperPlaneThickness() {return fSuperPlaneThickness;}
		/**
		* \brief Set the superplane thickness
		*/
		void SetSuperPlaneThickness(double value) {fSuperPlaneThickness= value;}
		/**
		* \brief Get the superplane height
		*/
		double GetSuperPlaneHeight() {return fSuperPlaneHeight;}
		/**
		* \brief Set the superplane height
		*/
		void SetSuperPlaneHeight(double value) {fSuperPlaneHeight= value;}
	
		/**
		* \brief Get the superplane distance
		*/
		std::vector<double> GetSuperPlaneDistance() {return fSuperPlaneDistance;}
		/**
		* \brief Set the superplane distance
		*/
		void SetSuperPlaneDistance(std::vector<double> vect) {fSuperPlaneDistance= vect;}

		/**
		* \brief Get the superplane depths
		*/
		std::vector<double> GetSuperPlaneDepth() {return fSuperPlaneDepth;}
		/**
		* \brief Set the superplane depths
		*/
		void SetSuperPlaneDepth(std::vector<double> vect) {fSuperPlaneDepth= vect;}
		/**
		* \brief Calculate the superplane depths according to the plane relative distance
		*/
		void CalculateSuperPlaneDepth();
		
		/**
		* \brief Get the plane XY distance
		*/
		double GetXYPlaneDistance() {return fPlaneXYDistance;}
		/**
		* \brief Set the plane XY distance
		*/
		void SetXYPlaneDistance(double value) {fPlaneXYDistance= value;}
		/**
		* \brief Get the plane XY title angle
		*/
		double GetXYPlaneTiltAngle() {return fPlaneXYTiltAngle;}
		/**
		* \brief Set the plane XY title angle
		*/
		void SetXYPlaneTiltAngle(double value) {fPlaneXYTiltAngle= value;}
		
		/**
		* \brief Get the strip size
		*/
		TVector3 GetStripSize() {return fStripSize;}
		/**
		* \brief Set the strip size
		*/
		void SetStripSize(TVector3 value) {fStripSize= value;}
		/**
		* \brief Get the strip groove size
		*/
		TVector3 GetStripGrooveSize() {return fStripGrooveSize;}
		/**
		* \brief Set the strip groove size
		*/
		void SetStripGrooveSize(TVector3 value) {fStripGrooveSize= value;}
		/**
		* \brief Get the strip coating size
		*/
		TVector3 GetStripCoatingSize() {return fStripCoatingSize;}
		/**
		* \brief Set the strip coating size
		*/
		void SetStripCoatingSize(TVector3 value) {fStripCoatingSize= value;}
		/**
		* \brief Get the strip coating thickness
		*/
		double GetStripCoatingThickness() {return fStripHousingThickness;}
		/**
		* \brief Set the strip coating thickness
		*/
		void SetStripCoatingThickness(double value) {fStripHousingThickness= value;}

		/**
		* \brief Get the fiber length
		*/
		double GetFiberLength() {return fFiberLength;}
		/**
		* \brief Set the fiber length
		*/
		void SetFiberLength(double value) {fFiberLength= value;}
		/**
		* \brief Get the fiber radius
		*/
		double GetFiberRadius() {return fFiberRadius;}
		/**
		* \brief Set the fiber radius
		*/
		void SetFiberRadius(double value) {fFiberRadius= value;}
		
		/**
		* \brief Get the pmt size
		*/
		TVector3 GetPMTSize() {return fPMTSize;}
		/**
		* \brief Set the pmt size
		*/
		void SetPMTSize(TVector3 value) {fPMTSize= value;}
		/**
		* \brief Get the pmt photocathode size
		*/
		TVector3 GetPhotocathodeSize() {return fPhotocathodeSize;}
		/**
		* \brief Set the pmt size
		*/
		void SetPhotocathodeSize(TVector3 value) {fPhotocathodeSize= value;}
		/**
		* \brief Get the optical coupling size
		*/
		TVector3 GetOptCouplingSize() {return fOptCouplingSize;}
		/**
		* \brief Set the optical coupling size
		*/
		void SetOptCouplingSize(TVector3 value) {fOptCouplingSize= value;}

	private:
	
		int fNStations;
		int fNStrips;
		int fNPlanes;

		double fSuperPlaneThickness;
		double fSuperPlaneHeight;
		std::vector<double> fSuperPlaneDistance;
		std::vector<double> fSuperPlaneDepth;
		double fPlaneXYDistance;
		double fPlaneXYTiltAngle;

		TVector3 fStripSize;
		TVector3 fStripGrooveSize;
		TVector3 fStripCoatingSize;
		double fStripHousingThickness;

		double fFiberLength;
		double fFiberRadius;
	
		TVector3 fPMTSize;
		TVector3 fPhotocathodeSize;
		TVector3 fOptCouplingSize;
	
		ClassDef(MuonDetector,1)
};

#ifdef __MAKECINT__
#pragma link C++ class MuonDetector+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<MuonDetector*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<MuonDetector>+;
#endif

#endif 
