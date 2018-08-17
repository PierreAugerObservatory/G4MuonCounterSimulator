/**
* @file MuonDetector.cc
* @class MuonDetector
* @brief Container for muon detector info
*
* @author S. Riggi
* @date 31/12/2010
*/

#include "MuonDetector.hh"
#include <G4UnitsTable.hh>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <TGeoTrack.h>
#include <TGraph.h>
#include <TGraph2D.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <assert.h>

using namespace std;

MuonDetector::MuonDetector(){

	fNStations= 0;
	fNStrips= 0;
	fNPlanes= 0;

	fSuperPlaneDepth.clear();
	fSuperPlaneDepth.resize(0);

}

MuonDetector::~MuonDetector(){

}


void MuonDetector::CalculateSuperPlaneDepth(){

	double d= 0.;

	fSuperPlaneDepth.clear();
	fSuperPlaneDepth.resize(0);

	for(unsigned int i=0;i<fSuperPlaneDistance.size();i++){
		d+= fSuperPlaneDistance[i]/cm;
		fSuperPlaneDepth.push_back(d);	
	}
	
}//end SetSuperPlaneDepth()

ClassImp(MuonDetector)

#ifdef __MAKECINT__
#pragma link C++ class MuonDetector+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<MuonDetector*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<MuonDetector>+;
#endif

