/**
* @file TrackReco.cc
* @class TrackReco
* @brief Reconstruct the particle tracks in the detector
*
* Reconstruct the particle tracks in the detector
* @author Dr. Simone Riggi
* @date 07/09/2010
*/

#include "TrackReco.hh"
#include "TrackFinder.hh"
#include "KFTrackFinder.hh"
#include "TrackPoint.hh"
#include "Track.hh"
#include "AnalysisConsts.hh"
#include "G4UnitsTable.hh"

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

using namespace std;

TrackReco::TrackReco(){

	fVerbosity= 1;

	fHitSigmaX= 0.;
	fHitSigmaY= 0.;
	fHitSigmaZ= 0.;

}

TrackReco::~TrackReco(){

}


bool TrackReco::TrackReconstruction(){

	cout<<"*******************************"<<endl;
	cout<<"***  TRACK RECONSTRUCTION   ***"<<endl;
	cout<<"*******************************"<<endl;
	//## Associate hit points to a track
	//## according to the least track Chi2
	/*
	TrackFinder* theTrackFinder= new TrackFinder;
	theTrackFinder->SetTrackPoints(fTrackPointCollection);
	theTrackFinder->SetHitSeparationInCluster(fHitSeparationInCluster);
	theTrackFinder->SetHitSigmaX(fHitSigmaX);
	theTrackFinder->SetHitSigmaY(fHitSigmaY);
	theTrackFinder->SetHitSigmaZ(fHitSigmaZ);
	theTrackFinder->SetSuperPlaneDepth(fSuperPlaneDepth);
	theTrackFinder->SetSuperPlaneSizeZ(fSuperPlaneSizeZ);
	theTrackFinder->SetVerbosity(1);
	//theTrackFinder->FindCluster();
	//theTrackFinder->FindTrack();
	theTrackFinder->TrackReconstructor();

	fRecTrackCollection= theTrackFinder->GetRecTrackCollection();
	fClusterTrackPointCollection= theTrackFinder->GetClusterTrackPoints();

	delete theTrackFinder;
	*/

	KFTrackFinder* theKFTrackFinder= new KFTrackFinder;
	theKFTrackFinder->SetTrackPoints(fTrackPointCollection);
	theKFTrackFinder->SetHitSeparationInCluster(fHitSeparationInCluster);
	theKFTrackFinder->SetHitSigmaX(fHitSigmaX);
	theKFTrackFinder->SetHitSigmaY(fHitSigmaY);
	theKFTrackFinder->SetHitSigmaZ(fHitSigmaZ);
	theKFTrackFinder->SetSuperPlaneDepth(fSuperPlaneDepth);
	theKFTrackFinder->SetSuperPlaneSizeZ(fSuperPlaneSizeZ);
	theKFTrackFinder->SetVerbosity(1);
	theKFTrackFinder->TrackReconstructor();

	fRecTrackCollection= theKFTrackFinder->GetRecTrackCollection();
	fClusterTrackPointCollection= theKFTrackFinder->GetClusterTrackPoints();

	delete theKFTrackFinder;

	return true;

}//close TrackReco::TrackFinder()





