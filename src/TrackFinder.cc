/**
* @file TrackFinder.cc
* @class TrackFinder
* @brief Reconstruct the particle tracks in the detector
*
* Reconstruct the particle tracks in the detector
* @author S. Riggi
* @date 07/09/2010
*/

#include "TrackFinder.hh"
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
#include <assert.h>

using namespace std;


bool TrackFinder::fIncludeMomentumInTracking;
bool TrackFinder::fTrackWithMSCorrection;
bool TrackFinder::fUseEnergyLossCorrectionInMS;
bool TrackFinder::fUseAverageEnergyInMS;
double TrackFinder::fAverageEnergyInMS;
bool TrackFinder::fUseTrueEnergyInMS;

double TrackFinder::fSuperPlaneSizeZ;
std::vector<double> TrackFinder::fSuperPlaneDepth;
double TrackFinder::fTheta;
double TrackFinder::fPhi;
int TrackFinder::fFitStatus;
double TrackFinder::fFitFCNMin;
double TrackFinder::fVertexX;
double TrackFinder::fVertexY;
double TrackFinder::fTx;
double TrackFinder::fTy;
double TrackFinder::fThetaStart;
double TrackFinder::fPhiStart;
double TrackFinder::fVertexXStart;
double TrackFinder::fVertexYStart;
double TrackFinder::fTxStart;
double TrackFinder::fTyStart;

double TrackFinder::fGenEnergy;
double TrackFinder::fGenTheta;
double TrackFinder::fGenPhi;
double TrackFinder::fGenVertexX;
double TrackFinder::fGenVertexY;
double TrackFinder::fGenTx;
double TrackFinder::fGenTy;

double TrackFinder::fThetaErr;
double TrackFinder::fPhiErr;
double TrackFinder::fVertexXErr;
double TrackFinder::fVertexYErr;
double TrackFinder::fTxErr;
double TrackFinder::fTyErr;
double TrackFinder::fThetaStartErr;
double TrackFinder::fPhiStartErr;
double TrackFinder::fVertexXStartErr;
double TrackFinder::fVertexYStartErr;
double TrackFinder::fTxStartErr;
double TrackFinder::fTyStartErr;

int TrackFinder::nHits;
int TrackFinder::nPar;
double TrackFinder::fMomentum;
double TrackFinder::fKinEnergy;
double TrackFinder::lenz[MAXHITNO];
double TrackFinder::Ls[MAXSCATNO];
double TrackFinder::Zs[MAXSCATNO];
double TrackFinder::ZsStart[MAXSCATNO];


TrackFinder::TrackFinder(){

	fRunNumber= 0;

	fVerbosity= 1;

	fHitSigmaX= 0.;
	fHitSigmaY= 0.;
	fHitSigmaZ= 0.;

	nTrack= 0;
	nHit= 0;
	nClusterHit= 0;

	fFitStatus= -1;
	fFitFCNMin= -999;
	fTheta= -999;
	fPhi= -999;
	fVertexX= -999;
	fVertexY= -999;
	fTx= -999;
	fTy= -999;
	fThetaStart= -999;
	fPhiStart= -999;
	fVertexXStart= -999;
	fVertexYStart= -999;
	fTxStart= -999;
	fTyStart= -999;

	fGenEnergy= -999;
	fGenTheta= -999;
	fGenPhi= -999;
	fGenVertexX= -999;
	fGenVertexY= -999;
	fGenTx= -999;
	fGenTy= -999;

	fThetaErr= -999;
	fPhiErr= -999;
	fVertexXErr= -999;
	fVertexYErr= -999;
	fTxErr= -999;
	fTyErr= -999;
	fThetaStartErr= -999;
	fPhiStartErr= -999;
	fVertexXStartErr= -999;
	fVertexYStartErr= -999;
	fTxStartErr= -999;
	fTyStartErr= -999;
	
	fKinEnergy= 5000.;//MeV

	fIncludeMomentumInTracking= false;
	fTrackWithMSCorrection= false;
	fUseEnergyLossCorrectionInMS= false;
	fUseAverageEnergyInMS= true;
	fAverageEnergyInMS= 5000.;//MeV
	fUseTrueEnergyInMS= false;

	fUseClusterInTracking= true;
	fSplitClustersInTracking= false;	
	fSplitClusterThreshold= 2;
	fRemoveFakeHitsInTracking= false;
}

TrackFinder::~TrackFinder(){

}


bool TrackFinder::TrackReconstructor(){

	//## Set track point uncertainty
	SetTrackPointUncertainty();

	//## Remove fake hits?
	if(fRemoveFakeHitsInTracking) 
		RemoveFakeHits();

	//## Find clusters
	bool clusterSearchStatus= FindCluster();
	if(!clusterSearchStatus) return false;

	//## Find track
	bool trackSearchStatus= FindTrack();
	if(!trackSearchStatus) return false;

	return true;

}//close TrackFinder::FindTrack()


bool TrackFinder::FindTrack(){

	//## Associate hit points to a track
	//## according to the least track Chi2
	//copy full track point collection to candidate track point collection
	if(fUseClusterInTracking){ 
		fCandidateTrackPointCollection= std::vector<TrackPoint>(fClusterTrackPointCollection.begin(),fClusterTrackPointCollection.end()); 
	}
	else{
		fCandidateTrackPointCollection= std::vector<TrackPoint>(fTrackPointCollection.begin(),fTrackPointCollection.end());	
	}


	fRecTrackCollection.clear();
	fRecTrackCollection.resize(0);	

	int rectrackCounter= 0;		
	bool IsTrackSearchFinished=false;
	while(!IsTrackSearchFinished){
		IsTrackSearchFinished= !Checker();
		if(IsTrackSearchFinished) {
			cout<<"TrackFinder::FindTrack(): End track search"<<endl;	
			break;
		}

		//## Create a combination of hit points
		CalculateCombinationMatrix();

		std::vector<Track*> CandidateTrackList;
		CandidateTrackList.clear();
		CandidateTrackList.resize(0);

		std::vector<int> TrackFitStatus;
		TrackFitStatus.clear();
		TrackFitStatus.resize(0);
	
		std::vector<double> TrackFitFCNMin;
		TrackFitFCNMin.clear();
		TrackFitFCNMin.resize(0);

		std::vector<double> ThetaFit;
		ThetaFit.clear();
		ThetaFit.resize(0);

		std::vector<double> PhiFit;
		PhiFit.clear();
		PhiFit.resize(0);

		std::vector<double> VertexXFit;
		VertexXFit.clear();
		VertexXFit.resize(0);

		std::vector<double> VertexYFit;
		VertexYFit.clear();
		VertexYFit.resize(0);

		std::vector<double> TxFit;
		TxFit.clear();
		TxFit.resize(0);

		std::vector<double> TyFit;
		TyFit.clear();
		TyFit.resize(0);

		std::vector<double> ThetaFitErr;
		ThetaFitErr.clear();
		ThetaFitErr.resize(0);

		std::vector<double> PhiFitErr;
		PhiFitErr.clear();
		PhiFitErr.resize(0);

		std::vector<double> VertexXFitErr;
		VertexXFitErr.clear();
		VertexXFitErr.resize(0);

		std::vector<double> VertexYFitErr;
		VertexYFitErr.clear();
		VertexYFitErr.resize(0);

		std::vector<double> TxFitErr;
		TxFitErr.clear();
		TxFitErr.resize(0);

		std::vector<double> TyFitErr;
		TyFitErr.clear();
		TyFitErr.resize(0);

  	int trackCounter=0;	

		//## Loop and fit all possible combinations
		for(int i = 0; i < fCombinationMatrix.size(); i++) {
			int rowIndex = 0;
   		int colIndex = 0;
			std::vector<int> row (fCombinationMatrix[i].begin(),fCombinationMatrix[i].end());

			Track* aCandidateTrack= new Track;
			
    	while(colIndex < fTrackPointCollectionInPlane.size()) {
				TrackPoint aPoint= fTrackPointCollectionInPlane[rowIndex++][row[colIndex++]];
				TVector3 aPointPos= aPoint.fPosition;
				TVector3 aPointPosErr= aPoint.fPositionErr;
				int aPointDetId= aPoint.fDetectorPlaneId;
				double aPointEdepX= aPoint.fEdepX;
				double aPointEdepY= aPoint.fEdepY;	
				double aPointTimeX= aPoint.fTimeX;
				double aPointTimeY= aPoint.fTimeY;
				//cout<<"  P=("<<(aPoint.fPosition).X()<<","<<(aPoint.fPosition).Y()<<","<<(aPoint.fPosition).Z()<<")"<<endl;
				//cout<<"  PErr=("<<(aPoint.fPositionErr).X()<<","<<(aPoint.fPositionErr).Y()<<","<<(aPoint.fPositionErr).Z()<<")"<<endl;
				//cout<<"  P EdepX="<<aPointEdepX<<"  EdepY="<<aPointEdepY<<"  TimeX="<<aPointTimeX<<"  TimeY="<<aPointTimeY<<endl;
				aCandidateTrack->AddPoint(aPointPos.X(),aPointPos.Y(),aPointPos.Z(),0.);//add first point to cand track
				aCandidateTrack->AddPointError(aPointPosErr.X(),aPointPosErr.Y(),aPointPosErr.Z(),0.);//add first point to cand track
				aCandidateTrack->AddPointDetId(aPointDetId);
				aCandidateTrack->AddTrackPoint(aPoint);
    	}
			CandidateTrackList.push_back(aCandidateTrack);
	
			//## fit this candidate track and get fit results	
			//TrackFitter(aCandidateTrack);
			//TrackFitterWithMS(aCandidateTrack);

			//## fast track fitter
			if(fTrackWithMSCorrection){
				fFitStatus= FastTrackFitterWithMS(aCandidateTrack);
			}
			else {
				fFitStatus= FastTrackFitter(aCandidateTrack);
			}
			
			
			//cout<<"Track "<<trackCounter<<"  FitFCN="<<fFitFCNMin<<"  FitStatus="<<fFitStatus<<endl;
			TrackFitStatus.push_back(fFitStatus);
			TrackFitFCNMin.push_back(fFitFCNMin);
			ThetaFit.push_back(fTheta);	
			PhiFit.push_back(fPhi);
			VertexXFit.push_back(fVertexX);	
			VertexYFit.push_back(fVertexY);	
			TxFit.push_back(fTx);
			TyFit.push_back(fTy);

			ThetaFitErr.push_back(fThetaErr);	
			PhiFitErr.push_back(fPhiErr);
			VertexXFitErr.push_back(fVertexXErr);	
			VertexYFitErr.push_back(fVertexYErr);	
			TxFitErr.push_back(fTxErr);
			TyFitErr.push_back(fTyErr);
					
			trackCounter++;
			
    	cout<<endl;
		}//close for combination matrix

		//## search best fit among candidate tracks
		double min=1.e+99;
		int bestTrackIndex= -1;
		for(unsigned int k=0;k<TrackFitFCNMin.size();k++){
			//cout<<"Track "<<k<<"  FitFCN="<<TrackFitFCNMin[k]<<"  FitStatus="<<TrackFitStatus[k]<<endl;
			if(TrackFitFCNMin[k]<min && TrackFitStatus[k]==0) {
				min= TrackFitFCNMin[k];
				bestTrackIndex= k;
			}
		}

		//## fill the list of reconstructed tracks
		CandidateTrackList[bestTrackIndex]->SetId(rectrackCounter);
		CandidateTrackList[bestTrackIndex]->SetTheta(ThetaFit[bestTrackIndex]);
		CandidateTrackList[bestTrackIndex]->SetPhi(PhiFit[bestTrackIndex]);
		CandidateTrackList[bestTrackIndex]->SetFitStatus(TrackFitStatus[bestTrackIndex]);
		CandidateTrackList[bestTrackIndex]->SetFitChi2(TrackFitFCNMin[bestTrackIndex]);
		CandidateTrackList[bestTrackIndex]->SetVertexPos(TVector3(VertexXFit[bestTrackIndex],VertexYFit[bestTrackIndex],0.));
		CandidateTrackList[bestTrackIndex]->SetTx(TxFit[bestTrackIndex]);
		CandidateTrackList[bestTrackIndex]->SetTy(TyFit[bestTrackIndex]);

		CandidateTrackList[bestTrackIndex]->SetThetaErr(ThetaFitErr[bestTrackIndex]);
		CandidateTrackList[bestTrackIndex]->SetPhiErr(PhiFitErr[bestTrackIndex]);
		CandidateTrackList[bestTrackIndex]->SetVertexPosErr(TVector3(VertexXFitErr[bestTrackIndex],VertexYFitErr[bestTrackIndex],0.));
		CandidateTrackList[bestTrackIndex]->SetTxErr(TxFitErr[bestTrackIndex]);
		CandidateTrackList[bestTrackIndex]->SetTyErr(TyFitErr[bestTrackIndex]);

		fRecTrackCollection.push_back(CandidateTrackList[bestTrackIndex]);
		rectrackCounter++;

		//## remove assigned track points from the list of candidate points
		//get track points of reconstructed track
		std::vector<int> RecTrackPointDetIdList= CandidateTrackList[bestTrackIndex]->GetPointDetId();
		double x,y,z,t;
		double xErr,yErr,zErr,tErr;
		std::vector<TrackPoint> RecTrackPointList;
		TrackPoint aRecTrackPoint;
		for(int k=0;k<CandidateTrackList[bestTrackIndex]->GetNpoints();k++){
			CandidateTrackList[bestTrackIndex]->GetPoint(k,x,y,z,t);
			CandidateTrackList[bestTrackIndex]->GetPointError(k,xErr,yErr,zErr,tErr);
			CandidateTrackList[bestTrackIndex]->GetTrackPoint(k,aRecTrackPoint);	
	
			TVector3 PointPos(x,y,z); 
			TVector3 PointPosErr(xErr,yErr,zErr); 
			TrackPoint aTrackPoint;
			aTrackPoint.fPosition= PointPos;
			aTrackPoint.fPositionErr= PointPosErr;

			RecTrackPointList.push_back(aTrackPoint);
		}
	
		std::vector<int> pointToBeErased;
		for(int k=0;k<fCandidateTrackPointCollection.size();k++){
			TrackPoint currentTrackPoint= fCandidateTrackPointCollection[k];
			for(unsigned l=0;l<RecTrackPointDetIdList.size();l++){
				TrackPoint recTrackPoint;
				recTrackPoint.fPosition= RecTrackPointList[l].fPosition;
				recTrackPoint.fPositionErr= RecTrackPointList[l].fPositionErr;
				recTrackPoint.fDetectorPlaneId= RecTrackPointDetIdList[l];
				
				if(currentTrackPoint==recTrackPoint) {
					pointToBeErased.push_back(k);
					break;
				}
			}
		}

		
		std::vector<TrackPoint>::iterator it= fCandidateTrackPointCollection.begin(); 
		for(unsigned k=0;k<pointToBeErased.size();k++){
			//cout<<"pointToBeErased "<<pointToBeErased[k]<<endl;
			//fCandidateTrackPointCollection.erase(fCandidateTrackPointCollection.begin()+pointToBeErased[k]);
			fCandidateTrackPointCollection.erase(it+pointToBeErased[k]);
			it--;
		}

	
		
		cout<<"*** TRACK FINDER RESULTS ***"<<endl;
		cout<<CandidateTrackList.size()<<" candidate tracks fitted"<<endl;
		cout<<" Track "<<bestTrackIndex<<" best candidate ("<<min<<")"<<endl;		
		cout<<fCandidateTrackPointCollection.size()<<" candidate points remaining"<<endl;
		cout<<"****************************"<<endl;
	}//close while search


	/*
	nTrack= (int)(fRecTrackCollection.size());
	nTrackHit= 0;
	for(int s=0;s<nTrack;s++){
		Track* currentRecTrack= fRecTrackCollection[s];
		RecTheta[s]= currentRecTrack->GetTheta();
		RecPhi[s]= currentRecTrack->GetPhi();
		RecX[s]= (currentRecTrack->GetVertexPos()).X();
		RecY[s]= (currentRecTrack->GetVertexPos()).Y();
		TrackFitStatus[s]= currentRecTrack->GetFitStatus();	
		TrackFitChi2[s]= currentRecTrack->GetFitChi2();
	
		double x,y,z,t;
		int id= currentRecTrack->GetId();
		for(int p=0;p<currentRecTrack->GetNpoints();p++){
			currentRecTrack->GetPoint(p,x,y,z,t);
			TrackHitX[nTrackHit]= x;
			TrackHitY[nTrackHit]= y;
			TrackHitZ[nTrackHit]= z;
			TrackId[nTrackHit]= id;	 	
			nTrackHit++;
		}
	}		
	*/

	return true;

}//close TrackFinder::TrackFinder()


void TrackFinder::GeomLine3D(double t, double *p, double &x, double &y, double &z) { 
	
	//##  Define the parameteric line equation in 3D
  //## 	A parameteric line is defined by 6 parameters (x0,y0,z0,z1,y1,z1) 
	//##  which are the coordinates of two points on the line) but only 4 are independent
  //##  can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
  x = p[0] + p[1]*t; 
  y = p[2] + p[3]*t;
  z = t; 

}//close TrackFinder::3DGeomLine()



double TrackFinder::GeomDistance3D(double x,double y,double z, double* p) {

	//## Calculate distance line-point  
  //## The distance line-point is D= | (xp-x0) cross  ux | 
  //## where ux is direction of line and x0 is a point in the line (like t = 0) 
  TVector3 xp(x,y,z); 
  TVector3 x0(p[0], p[2], 0. ); 
  TVector3 x1(p[0] + p[1], p[2] + p[3], 1. ); 
  TVector3 u = (x1-x0).Unit(); 
  double d2 = ((xp-x0).Cross(u)) .Mag2(); 
  
	return d2; 

}//close TrackFinder::3DGeomDistance()

void TrackFinder::TrackChi2Fcn(int& nPar, double* const grad,double& value, double* par,const int iFlag){
    
	// the TGraph must be a global variable
  TGraph2D* trackPointGraph = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
  assert(trackPointGraph != 0);
  
	double * x = trackPointGraph->GetX();
  double * y = trackPointGraph->GetY();
  double * z = trackPointGraph->GetZ();
  int Npoints = trackPointGraph->GetN();
   
	value = 0;
  for(int i= 0; i<Npoints;++i) { 
		double d= GeomDistance3D(x[i],y[i],z[i],par); 
    value += d;
  }

	//cout<<"## SumGeomDistance for this track ="<< value <<" ##"<<endl;

}//close TrackFinder::SumGeomDistance3D()

void TrackFinder::TrackFitter(Track* aTrack){

	double x,y,z,time;
	int Npoints= aTrack->GetNpoints();

	//## graph with candidate track points 	
	TGraph2D* trackPointGraph = new TGraph2D();

  //## generate graph with the 3d points of the candidate track
  for(int i=0; i<Npoints; i++) {
  	aTrack->GetPoint(i,x,y,z,time);
		trackPointGraph->SetPoint(i,x,y,z);
  }
   
	//## fit the graph now 
  TVirtualFitter* fit = TVirtualFitter::Fitter(0,4);
  fit->SetObjectFit(trackPointGraph);
  fit->SetFCN(TrackFinder::TrackChi2Fcn);
   
  double arglist[10];
  arglist[0] = 3;
  fit->ExecuteCommand("SET PRINT",arglist,1);
  
  double pStart[4] = {1,1,1,1};
  fit->SetParameter(0,"x0",pStart[0],0.01,0,0);
  fit->SetParameter(1,"Ax",pStart[1],0.01,0,0);
  fit->SetParameter(2,"y0",pStart[2],0.01,0,0);
  fit->SetParameter(3,"Ay",pStart[3],0.01,0,0);
    
  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  fFitStatus= fit->ExecuteCommand("MIGRAD",arglist,2);

  //if (minos) min->ExecuteCommand("MINOS",arglist,0);
  int nvpar,nparx; 
  double amin,edm, errdef;
  fit->GetStats(amin,edm,errdef,nvpar,nparx);
	fFitFCNMin= amin;
	cout<<"FCN min="<<fFitFCNMin<<" FitStatus="<<fFitStatus<<endl;
  fit->PrintResults(1,amin);
  
	//trackPointGraph->Draw("p0");


  //## get fit parameters
  double parFit[4];
  for(int i=0;i<4;++i) parFit[i] = fit->GetParameter(i);

	//## get track direction
	TVector3 x0(parFit[0], parFit[2], 0.); 
  TVector3 x1(parFit[0] + parFit[1], parFit[2] + parFit[3], 1.); 
  TVector3 direction = (x1-x0).Unit();
	TVector3 oppositeDirection = TVector3(-direction.X(),-direction.Y(),-direction.Z());  
  //double theta= oppositeDirection.Theta()*180./TMath::Pi();
	//double phi= oppositeDirection.Phi()*180./TMath::Pi();
	fTheta= direction.Theta()*180./TMath::Pi();
	fPhi= direction.Phi()*180./TMath::Pi();
	fVertexX= x0.X();
	fVertexY= x0.Y();	
	fTx= parFit[2];
	fTy= parFit[3];

	cout<<"RecTheta[deg]="<<fTheta<<"  RecPhi[deg]="<<fPhi<<endl;
	   
  //## draw the fitted line
	/*
  int n = 1000;
  double t0 = 0;
  double dt = 10;
  TPolyLine3D* l= new TPolyLine3D(n);
  for(int i=0;i<n;++i) {
		double t = t0+ dt*i/n;
   	double x,y,z;
    GeomLine3D(t,parFit,x,y,z);
    l->SetPoint(i,x,y,z);
  }
  l->SetLineColor(kRed);
  //l->Draw("same");
  */
  	
}//close TrackFinder::TrackFitter()


void TrackFinder::TrackFitterWithMS(Track* aTrack){

	double x,y,z,time;
	double xErr,yErr,zErr,timeErr;
	int Npoints= aTrack->GetNpoints();
	//if (Npoints<3) return false;

	//## fit the graph now 
  TVirtualFitter* fit = TVirtualFitter::Fitter(0,4);
  fit->SetObjectFit(aTrack);
  fit->SetFCN(TrackFinder::TrackChi2WithMSFcn);
   
  double arglist[10];
  arglist[0] = 3;
  fit->ExecuteCommand("SET PRINT",arglist,1);
  
  double pStart[4] = {1,1,1,1};
  fit->SetParameter(0,"x0",pStart[0],0.01,0,0);
  fit->SetParameter(1,"Ax",pStart[1],0.01,0,0);
  fit->SetParameter(2,"y0",pStart[2],0.01,0,0);
  fit->SetParameter(3,"Ay",pStart[3],0.01,0,0);
    
  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  fFitStatus= fit->ExecuteCommand("MIGRAD",arglist,2);

  //if (minos) min->ExecuteCommand("MINOS",arglist,0);
  int nvpar,nparx; 
  double amin,edm, errdef;
  fit->GetStats(amin,edm,errdef,nvpar,nparx);
	fFitFCNMin= amin;
	cout<<"FCN min="<<fFitFCNMin<<" FitStatus="<<fFitStatus<<endl;
  fit->PrintResults(1,amin);
  

  //## get fit parameters
  double parFit[4];
  for(int i=0;i<4;++i) parFit[i] = fit->GetParameter(i);

	//## get track direction
	TVector3 x0(parFit[0], parFit[1], 0.); 
  TVector3 x1(parFit[0] + parFit[2], parFit[1] + parFit[3], 1.); 
  TVector3 direction = (x1-x0).Unit();
	
	//fTheta= Theta*180./TMath::Pi();
	//fPhi= Phi*180./TMath::Pi();
	fTheta= direction.Theta()*180./TMath::Pi();
	fPhi= direction.Phi()*180./TMath::Pi();

	fVertexX= x0.X();
	fVertexY= x0.Y();
	fTx= parFit[2];
	fTy= parFit[3];

	cout<<"*** FAST TRACK RECO WITH MS ***"<<endl;
	cout<<"Theta [deg]="<<fTheta<<"  Phi [deg]="<<fPhi<<endl;
	cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 	cout<<"*******************************"<<endl;
  
	//if(fFitStatus==0) return true;
	//else return false;

	//return true;
	
}//close TrackFinder::TrackFitterWithMS()



int TrackFinder::FastTrackFitter(Track* aTrack){

	//cout<<"**************************"<<endl;
	//cout<<"*** FAST TRACK FITTER  ***"<<endl;
	//cout<<"**************************"<<endl;
	double x,y,z,time;
	double xErr,yErr,zErr,timeErr;
	//int nHits= aTrack->GetNpoints();
	nHits= aTrack->GetNpoints();

	//## Get hits and set uncertainties
	double hits[nHits][3];
	double sigma[nHits][3];
	for(int i=0; i<nHits; i++) {
  	aTrack->GetPoint(i,x,y,z,time);
		hits[i][0] = x;
		hits[i][1] = y;
		hits[i][2] = z;
		
		aTrack->GetPointError(i,xErr,yErr,zErr,timeErr);
		sigma[i][0] = xErr;
		sigma[i][1] = yErr;
		sigma[i][2] = zErr;
	}	


	//## Reset fit values
	double Chi2 = 0;
	double Ndof = 2*nHits-4;
	double Theta = 0.0;
	double Phi = 0.0;
	double U[3];
	U[0] = 0.0;
	U[1] = 0.0;
	U[2] = 1.0;
	double P0[3];
	P0[0] = 0.0;
	P0[1] = 0.0;
	P0[2] = 0.0;
	
	//## Consistency check on the number of hits
	if (nHits<3) return -1;
	P0[2] = hits[0][2];

	//## Lenghts (Zi-Z0)
	//double lenz[nHits];
	for(int i=0;i<nHits;i++) lenz[i] = hits[i][2] - P0[2];
 
  //## Calculate derivative matrix F & G
	//## Derivative of predictions f with respect to the four track parameter
	const int idim = 4;
	double d[2*nHits][idim];
	for(int i=0;i<nHits;i++) {
		int ix = i;
		int iy = i+nHits;
		for(int j=0;j<idim;j++) { d[ix][j] = 0; d[iy][j] = 0;}
		d[ix][0] = 1.;
		d[iy][1] = 1.;
		d[ix][2] = lenz[i];
		d[iy][3] = lenz[i];
	}
 
  //## Calculate matrix product F*S_x*x + G*S_y*y
  double dx[idim];
	for(int j=0;j<idim;j++) {
		dx[j] = 0.;
		for(int l=0;l<nHits;l++) {
			dx[j] += d[l][j]/sigma[l][0]/sigma[l][0]*hits[l][0];
			dx[j] += d[l+nHits][j]/sigma[l][1]/sigma[l][1]*hits[l][1];
		}
	}
 
	//## Calculate matrix product (F*S_x*F + G*S_y*G): Covariance Matrix
	double par[idim];
	double parErr[idim];
	double InvCov[idim][idim];
	for (int j=0;j<idim;j++) {
		for (int k=0;k<idim;k++) {
			InvCov[j][k] = 0.;
			for (int l=0;l<nHits;l++) {
      	InvCov[j][k] += d[l][j]/sigma[l][0]/sigma[l][0]*d[l][k];
				InvCov[j][k] += d[l+nHits][j]/sigma[l][1]/sigma[l][1]*d[l+nHits][k];
			} 
		}
	}
       
  //## Calculate matrix product (F*S_x*F + G*S_y*G)**{-1}: Invert of Covariance Matrix
  double determ = 0.0;
  TMatrixD ParaCovariance = TMatrixD(idim,idim,(double*)InvCov," ");
  ParaCovariance = ParaCovariance.Invert(&determ);
  if (determ<=0) return -1;
 
  //## Find track solution & errors
 	for (int k=0;k<idim;k++) {
		par[k] = 0.;
		parErr[k]= sqrt( ParaCovariance(k,k) );
   	for (int i=0;i<idim;i++) {
	  	par[k] += ParaCovariance(k,i)*dx[i];
   	}
 	}
 
	//## Chi2
	Chi2 = 0.;
	for (int l=0;l<nHits;l++) {
		double xl = hits[l][0];
		double yl = hits[l][1];
		for (int k=0;k<idim;k++) {
			xl -= d[l][k]*par[k];
			yl -= d[l+nHits][k]*par[k];
		}
		Chi2 += xl/sigma[l][0]/sigma[l][0]*xl + yl/sigma[l][1]/sigma[l][1]*yl;
	}
 
  //## Final result
  P0[0] = par[0];
  P0[1] = par[1];
  Phi = atan2(par[3],par[2]);
  Theta = acos(-sqrt(1-par[2]*par[2]-par[3]*par[3]));
  U[0] = sin(Theta)*cos(Phi);
  U[1] = sin(Theta)*sin(Phi);
  U[2] = cos(Theta);

	TVector3 x0(par[0], par[1], 0.); 
  TVector3 x1(par[0] + par[2], par[1] + par[3], 1.); 
  TVector3 direction = (x1-x0).Unit();

	//fTheta= Theta*180./TMath::Pi();
	//fPhi= Phi*180./TMath::Pi();
	fTheta= direction.Theta()*180./TMath::Pi();
	fPhi= direction.Phi()*180./TMath::Pi();

	fVertexX= x0.X();
	fVertexY= x0.Y();	
	fTx= par[2];
	fTy= par[3];

	fVertexXErr= parErr[0];
	fVertexYErr= parErr[1];	
	fTxErr= parErr[2];
	fTyErr= parErr[3];

	//## theta & phi par error
	//## theta= atan2(sqrt(tx^2+ty^2),1)
	//## phi= atan2(ty,tx)
	//## do the error propagation as FCFt, with F derivative matrix d(theta)/d(par_i), d(phi)/d(par_i) 
	TMatrixD DTheta(1,idim);	
	DTheta(0,0)= 0;
	DTheta(0,1)= 0;	
	DTheta(0,2)= fTx/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
 	DTheta(0,3)= fTy/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
	TMatrixD DThetaT(idim,1);		
	DThetaT.Transpose(DTheta);

	TMatrixD DPhi(1,idim);	
	DPhi(0,0)= 0;
	DPhi(0,1)= 0;		
	DPhi(0,2)= -fTy/(fTx*fTx+fTy*fTy);
 	DPhi(0,3)= fTx/(fTx*fTx+fTy*fTy);
	TMatrixD DPhiT(idim,1);		
	DPhiT.Transpose(DPhi);

	TMatrixD DThetaMultCov(1,idim);
	DThetaMultCov.Mult(DTheta,ParaCovariance);
	TMatrixD DPhiMultCov(1,idim);
	DPhiMultCov.Mult(DPhi,ParaCovariance);
		
	TMatrixD ThetaUncertainty(1,1);
	ThetaUncertainty.Mult(DThetaMultCov,DThetaT);

	TMatrixD PhiUncertainty(1,1);
	PhiUncertainty.Mult(DPhiMultCov,DPhiT);

	fThetaErr= sqrt(ThetaUncertainty(0,0))* 180./TMath::Pi();
	fPhiErr= sqrt(PhiUncertainty(0,0))* 180./TMath::Pi();

	fFitFCNMin= Chi2;

	cout<<"*** FAST TRACK RECO ***"<<endl;
	cout<<"Theta [deg]="<<fTheta<<" +- "<< fThetaErr << "  Phi [deg]="<<fPhi<<" +- "<<fPhiErr<<endl;
	cout<<"X0 [cm]="<< fVertexX <<" +- "<< fVertexXErr << "  Y0[cm]="<<fVertexY<<" +- "<< fVertexYErr<<endl;
	cout<<"Tx="<<fTx<<" +- "<< fTxErr <<"  Ty="<<fTy<<" +- "<< fTyErr <<endl;
	cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 	cout<<"***********************"<<endl;

  return 0;

}//close FastTrackFitter()


int TrackFinder::FastTrackFitterWithMS(Track* aTrack){

	//cout<<"*********************************"<<endl;
	//cout<<"*** FAST TRACK FITTER WITH MS ***"<<endl;
	//cout<<"*********************************"<<endl;
	double x,y,z,time;
	double xErr,yErr,zErr,timeErr;
	nHits= aTrack->GetNpoints();
	
	//## Consistency checks
	if (nHits<3) return -1;	
	if(fSuperPlaneDepth.size()<3) return -1;
	
	//## Reset fit values
	double Chi2 = 0;
	double Ndof = 2*nHits-4;
	double Theta = 0.0;
	double Phi = 0.0;
	double U[3];
	U[0] = 0.0;
	U[1] = 0.0;
	U[2] = 1.0;
	double P0[3];
	P0[0] = 0.0;
	P0[1] = 0.0;
	P0[2] = 0.0;

	//## Get hits and set uncertainties
	double hits[nHits][3];
	double sigma[nHits][3];
	double m[2*nHits][1];
	for(int i=0; i<nHits; i++) {
  	aTrack->GetPoint(i,x,y,z,time);
		hits[i][0] = x;
		hits[i][1] = y;
		hits[i][2] = z;
		
		aTrack->GetPointError(i,xErr,yErr,zErr,timeErr);
		sigma[i][0] = xErr;
		sigma[i][1] = yErr;
		sigma[i][2] = zErr;
		cout<<"Z["<<i<<"]="<< hits[i][2]<<endl;
		cout<<"Sigma["<<i<<"]="<< sigma[i][0]<<"  "<<sigma[i][1]<<endl;

		int ix = i;
		int iy = i+nHits;
		m[ix][0]= x;
		m[iy][0]= y;
	}	

	P0[2] = hits[0][2];

	//## Lenghts (Zi-Z0)
	//double lenz[nHits];
	for(int i=0;i<nHits;i++) lenz[i] = hits[i][2] - P0[2];

	//## Scatterer positions
	const int nS= 3;
	double Ls[nS];
	double Zs[nS];
	double ZsStart[nS];
	
	Ls[0]= fSuperPlaneDepth[0]-fSuperPlaneSizeZ/2.;
	Ls[1]= fSuperPlaneDepth[1]-fSuperPlaneDepth[0]-fSuperPlaneSizeZ;
	Ls[2]= fSuperPlaneDepth[2]-fSuperPlaneDepth[1]-fSuperPlaneSizeZ;
	ZsStart[0]= 0.;
	ZsStart[1]= -(fSuperPlaneDepth[0]+fSuperPlaneSizeZ/2.);
	ZsStart[2]= -(fSuperPlaneDepth[1]+fSuperPlaneSizeZ/2.);
	for(int i=0;i<nS;i++){
		Zs[i]= -(fSuperPlaneDepth[i]-fSuperPlaneSizeZ/2.);
		cout<<"Ls["<<i<<"]="<<Ls[i]<<"  Zs["<<i<<"]="<<Zs[i]<<"  ZsStart["<<i<<"]="<<ZsStart[i]<<endl;
	}

	//## Get measurement matrix M
	TMatrixD M = TMatrixD(2*nHits,1,(double*)m," ");
	cout<<"## Measurement Matrix M ##"<<endl;
	M.Print();
	cout<<endl;


	//## Calculate derivative matrix F & G
	//## Derivative of predictions f with respect to the four track parameter
	const int idim = 4;
	double d[2*nHits][idim];
	for(int i=0;i<nHits;i++) {
		int ix = i;
		int iy = i+nHits;
		for(int j=0;j<idim;j++) { d[ix][j] = 0; d[iy][j] = 0;}
		d[ix][0] = 1.;
		d[iy][1] = 1.;
		d[ix][2] = lenz[i];
		d[iy][3] = lenz[i];
	}
	TMatrixD F = TMatrixD(2*nHits,idim,(double*)d," ");
	cout<<"## Derivative Matrix F ##"<<endl;
	F.Print();
	cout<<endl;

	//## Calculate transpose of derivative matrix Ft
	TMatrixD Ft= TMatrixD(idim,2*nHits);
	cout<<"## Transpose Derivative Matrix Ft ##"<<endl;
	Ft.Transpose(F);
	Ft.Print();
	cout<<endl;

	//## Calculate Covariance Matrix V
	double v[2*nHits][2*nHits];
	double momentum= 1000.;//in MeV/c
	double Lr= 21.997;// radiation length of Malargue soil in g/cm^2 
	double H= pow(13.6/momentum,2)/Lr; //MISSING: PUT TRACK INCLINATION!!!

	for(int i=0;i<2*nHits;i++){
		for(int j=0;j<2*nHits;j++) v[i][j]= 0.;
	}
	
	//## diagonal elements
	for(int i=0;i<nHits;i++) {
		int ix = i;
		int iy = i+nHits;
		v[ix][ix]= pow(sigma[i][0],2);
		v[iy][iy]= pow(sigma[i][1],2);
		for(int s=0;s<nS;s++){
			double thisMSContribution= MSCovariance(hits[i][2],hits[i][2], Zs[s], ZsStart[s], Ls[s], H);
			v[ix][ix]+= thisMSContribution; 
			v[iy][iy]+= thisMSContribution;
			//cout<<"ix="<<ix<<"  iy="<<iy<<"  hits[i][2]="<<hits[i][2]<<"  Zs[="<<s<<"]="<<Zs[s]<<"  ZsStart[="<<s<<"]="<<ZsStart[s]<<"  thisMSContribution="<<thisMSContribution<<"  v="<<v[ix][ix]<<endl;
		}
		//cout<<"v["<<ix<<","<<ix<<"]= "<<v[ix][ix]<<endl;
		//cout<<"v["<<iy<<","<<iy<<"]= "<<v[iy][iy]<<endl;
	}

	//## non-diagonal elements
	v[0][1]= MSCovariance(hits[0][2],hits[1][2], Zs[0], ZsStart[0], Ls[0], H);
	v[1][0]= v[0][1];
	v[0+nHits][1+nHits]= v[0][1];
	v[1+nHits][0+nHits]= v[0+nHits][1+nHits];
	//cout<<"v[0][1]="<<v[0][1]<<endl;
	//cout<<"v[0+nHits][1+nHits]="<<v[0+nHits][1+nHits]<<endl;

	v[0][2]= MSCovariance(hits[0][2],hits[2][2], Zs[0], ZsStart[0], Ls[0], H);
	v[2][0]= v[0][2];
	v[0+nHits][2+nHits]= v[0][2];
	v[2+nHits][0+nHits]= v[0+nHits][2+nHits];
	//cout<<"v[0][2]="<<v[0][2]<<endl;
	//cout<<"v[0+nHits][2+nHits]="<<v[0+nHits][2+nHits]<<endl;

	v[1][2]= MSCovariance(hits[1][2],hits[2][2], Zs[0], ZsStart[0], Ls[0], H) + MSCovariance(hits[1][2],hits[2][2], Zs[1], ZsStart[1], Ls[1], H);
	v[2][1]= v[1][2];
	v[1+nHits][2+nHits]= v[1][2];
	v[2+nHits][1+nHits]= v[1+nHits][2+nHits];
	//cout<<"v[1][2]="<<v[1][2]<<endl;
	//cout<<"v[1+nHits][2+nHits]="<<v[1+nHits][2+nHits]<<endl;

	cout<<endl;
	for(int i=0;i<2*nHits;i++) {
		cout<<"v["<<i<<","<<i<<"]= "<<v[i][i]<<endl;
	}

	TMatrixD V = TMatrixD(2*nHits,2*nHits,(double*)v," ");
	cout<<"## Covariance Matrix V ##"<<endl;
	V.Print();
	cout<<endl;

	//## Invert the Covariance Matrix
  double determ = 0.0;
  TMatrixD VInv = TMatrixD(2*nHits,2*nHits,(double*)v," ");
  VInv = VInv.Invert(&determ);
  if (determ<=0) {
		cerr<<"TrackFinder::FastTrackFitterWithMS(): Covariance Matrix inversion failed!"<<endl;
		return -1;
	}
	cout<<"## Inverted Covariance Matrix V ##"<<endl;
	VInv.Print();
	cout<<endl;

	//## Test inverse
	TMatrixD IdentityMatrix= TMatrixD(2*nHits,2*nHits);
	IdentityMatrix.Mult(V,VInv);
	cout<<"## Test VxVInv=I ##"<<endl;
	IdentityMatrix.Print();
	cout<<endl;


	//## Calculate Matrix (Ft VInv F)	= C covariance
	TMatrix FtVInvProduct= TMatrix(idim,2*nHits); 
	FtVInvProduct.Mult(Ft,VInv);
	cout<<"## Ft x VInv ##"<<endl;
	FtVInvProduct.Print();
	cout<<endl;

	TMatrix C= TMatrixD(idim,idim);
	C.Mult(FtVInvProduct,F);
	cout<<"## Covariance C ##"<<endl;
	C.Print();
	cout<<endl;

	//## Invert the covariance C
	determ = 0.0;
  TMatrixD CInv = C;
  CInv = CInv.Invert(&determ);
  if (determ<=0) {
		cerr<<"TrackFinder::FastTrackFitterWithMS(): Final Covariance Matrix inversion failed!"<<endl;
		return -1;
	}
	cout<<"## Inverted Covariance Matrix C ##"<<endl;
	CInv.Print();
	cout<<endl;

	//## Calculate CInv Ft VInv M
	TMatrixD CInvFtProduct= TMatrixD(idim,2*nHits);
	CInvFtProduct.Mult(CInv,Ft);
	cout<<"## CInv x Ft ##"<<endl;
	CInvFtProduct.Print();
	cout<<endl;

	TMatrixD CInvFtVInvProduct= TMatrixD(idim,2*nHits);
	CInvFtVInvProduct.Mult(CInvFtProduct,VInv);
	cout<<"## CInv x Ft x VInv ##"<<endl;
	CInvFtVInvProduct.Print();
	cout<<endl;


	TMatrixD FitParam= TMatrixD(idim,1);
	FitParam.Mult(CInvFtVInvProduct,M);

	cout<<"## Fit Params ##"<<endl;
	FitParam.Print();
	cout<<endl;

	//## Find track solution
	double par[idim];
	double parErr[idim];
 	for (int k=0;k<idim;k++) {
		par[k] = FitParam(k,0);
		parErr[k]= sqrt( CInv(k,k) );
 	}
 	
	
	//## Calculate (M-F lambda)
	TMatrixD FParProduct= TMatrixD(2*nHits,1);
	FParProduct.Mult(F,FitParam);
	cout<<"## F x Lambda ##"<<endl;
	FParProduct.Print();
	cout<<endl;


	TMatrixD R= TMatrixD(2*nHits,1);
	R= M-FParProduct;
	cout<<"## M- F x Lambda ##"<<endl;
  R.Print();
	cout<<endl;


	TMatrixD Rt= TMatrixD(1,2*nHits);
	Rt.Transpose(R);
	cout<<"## (M- F x Lambda)T ##"<<endl;
  Rt.Print();
	cout<<endl;


	//## Calculate Chi2
	TMatrixD RtVInvProduct= TMatrixD(1,2*nHits);
	RtVInvProduct.Mult(Rt,VInv);
	cout<<"## (M- F x Lambda)T x VInv ##"<<endl;
  RtVInvProduct.Print();
	cout<<endl;
	

	TMatrixD Chi2Matrix= TMatrixD(1,1);
	Chi2Matrix.Mult(RtVInvProduct,R); 
	cout<<"## (M- F x Lambda)T x VInv x (M- F x Lambda) ##"<<endl;
  Chi2Matrix.Print();
	cout<<endl;

	Chi2 = Chi2Matrix(0,0);
	
  //## Final result
  P0[0] = par[0];
  P0[1] = par[1];
  Phi = atan2(par[3],par[2]);
  Theta = acos(-sqrt(1-par[2]*par[2]-par[3]*par[3]));
  U[0] = sin(Theta)*cos(Phi);
  U[1] = sin(Theta)*sin(Phi);
  U[2] = cos(Theta);

	TVector3 x0(par[0], par[1], 0.); 
  TVector3 x1(par[0] + par[2], par[1] + par[3], 1.); 
  TVector3 direction = (x1-x0).Unit();
	
	//fTheta= Theta*180./TMath::Pi();
	//fPhi= Phi*180./TMath::Pi();
	fTheta= direction.Theta()*180./TMath::Pi();
	fPhi= direction.Phi()*180./TMath::Pi();

	fVertexX= x0.X();
	fVertexY= x0.Y();	
	fTx= par[2];
	fTy= par[3];

	//TVector3 oppositedirection = TVector3(-direction.X(),-direction.Y(),-direction.Z());
	//fTheta= oppositedirection.Theta()*180./TMath::Pi();
	//fPhi= oppositedirection.Phi()*180./TMath::Pi();

	fVertexXErr= parErr[0];
	fVertexYErr= parErr[1];	
	fTxErr= parErr[2];
	fTyErr= parErr[3];

	//## theta & phi par error
	//## theta= atan2(sqrt(tx^2+ty^2),1)
	//## phi= atan2(ty,tx)
	//## do the error propagation as FCFt, with F derivative matrix d(theta)/d(par_i), d(phi)/d(par_i) 
	TMatrixD DTheta(1,idim);	
	DTheta(0,0)= 0;
	DTheta(0,1)= 0;	
	DTheta(0,2)= fTx/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
 	DTheta(0,3)= fTy/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
	TMatrixD DThetaT(idim,1);		
	DThetaT.Transpose(DTheta);

	TMatrixD DPhi(1,idim);	
	DPhi(0,0)= 0;
	DPhi(0,1)= 0;		
	DPhi(0,2)= -fTy/(fTx*fTx+fTy*fTy);
 	DPhi(0,3)= fTx/(fTx*fTx+fTy*fTy);
	TMatrixD DPhiT(idim,1);		
	DPhiT.Transpose(DPhi);

	TMatrixD DThetaMultCov(1,idim);
	DThetaMultCov.Mult(DTheta,CInv);
	TMatrixD DPhiMultCov(1,idim);
	DPhiMultCov.Mult(DPhi,CInv);
		
	TMatrixD ThetaUncertainty(1,1);
	ThetaUncertainty.Mult(DThetaMultCov,DThetaT);

	TMatrixD PhiUncertainty(1,1);
	PhiUncertainty.Mult(DPhiMultCov,DPhiT);

	fThetaErr= sqrt(ThetaUncertainty(0,0))* 180./TMath::Pi();
	fPhiErr= sqrt(PhiUncertainty(0,0))* 180./TMath::Pi();
	


	fFitFCNMin= Chi2;

	cout<<"*** FAST TRACK RECO WITH MS ***"<<endl;
	cout<<"Theta [deg]="<<fTheta<<" +- "<< fThetaErr << "  Phi [deg]="<<fPhi<<" +- "<<fPhiErr<<endl;
	cout<<"X0 [cm]="<< fVertexX <<" +- "<< fVertexXErr << "  Y0[cm]="<<fVertexY<<" +- "<< fVertexYErr<<endl;
	cout<<"Tx="<<fTx<<" +- "<< fTxErr <<"  Ty="<<fTy<<" +- "<< fTyErr <<endl;
	cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 	cout<<"***********************"<<endl;

	return 0;

}//close FastTrackFitterWithMS()


void TrackFinder::TrackChi2WithMSFcn(int& nPar, double* const grad,double& value, double* par,const int iFlag){
    
	// the TGraph must be a global variable
  Track* aTrack = dynamic_cast<Track*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
  assert(aTrack != 0);
 
	double x,y,z,time;
	double xErr,yErr,zErr,timeErr;
	int nHits= aTrack->GetNpoints();

	//## Reset fit values
	double Chi2 = 0;
	double Ndof = 2*nHits-4;
	double Theta = 0.0;
	double Phi = 0.0;
	double U[3];
	U[0] = 0.0;
	U[1] = 0.0;
	U[2] = 1.0;
	double P0[3];
	P0[0] = 0.0;
	P0[1] = 0.0;
	P0[2] = 0.0;

	//## Get hits and set uncertainties
	double hits[nHits][3];
	double sigma[nHits][3];
	double xHit[nHits];
	double yHit[nHits];
	double zHit[nHits];
	double xHitSigma[nHits];
	double yHitSigma[nHits];
	double zHitSigma[nHits];
	double m[2*nHits][1];
	for(int i=0; i<nHits; i++) {
  	aTrack->GetPoint(i,x,y,z,time);
		hits[i][0] = x;
		hits[i][1] = y;
		hits[i][2] = z;
		xHit[i]= x;
		yHit[i]= y;
		zHit[i]= z;
		
		aTrack->GetPointError(i,xErr,yErr,zErr,timeErr);
		sigma[i][0] = xErr;
		sigma[i][1] = yErr;
		sigma[i][2] = zErr;
		xHitSigma[i]= xErr;
		yHitSigma[i]= yErr;
		zHitSigma[i]= zErr;
		//cout<<"Z["<<i<<"]="<< hits[i][2]<<endl;
		//cout<<"Sigma["<<i<<"]="<< sigma[i][0]<<"  "<<sigma[i][1]<<endl;

		int ix = i;
		int iy = i+nHits;
		m[ix][0]= x;
		m[iy][0]= y;
	}	

	
	
	//## Consistency check on the number of hits
	
	P0[2] = hits[0][2];

	//## Lenghts (Zi-Z0)
	double lenz[nHits];
	for(int i=0;i<nHits;i++) lenz[i] = hits[i][2] - P0[2];

	//## Scatterer positions
	const int nS= 3;
	double Ls[nS];
	double Zs[nS];
	double ZsStart[nS];
	
	Ls[0]= fSuperPlaneDepth[0]-fSuperPlaneSizeZ/2.;
	Ls[1]= fSuperPlaneDepth[1]-fSuperPlaneDepth[0]-fSuperPlaneSizeZ;
	Ls[2]= fSuperPlaneDepth[2]-fSuperPlaneDepth[1]-fSuperPlaneSizeZ;
	ZsStart[0]= 0.;
	//ZsStart[1]= fSuperPlaneDepth[0]+fSuperPlaneSizeZ/2.;
	//ZsStart[2]= fSuperPlaneDepth[1]+fSuperPlaneSizeZ/2.;
	ZsStart[1]= -(fSuperPlaneDepth[0]+fSuperPlaneSizeZ/2.);
	ZsStart[2]= -(fSuperPlaneDepth[1]+fSuperPlaneSizeZ/2.);
	for(int i=0;i<nS;i++){
		//Zs[i]= fSuperPlaneDepth[i]-fSuperPlaneSizeZ/2.;
		Zs[i]= -(fSuperPlaneDepth[i]-fSuperPlaneSizeZ/2.);
		//cout<<"Ls["<<i<<"]="<<Ls[i]<<"  Zs["<<i<<"]="<<Zs[i]<<"  ZsStart["<<i<<"]="<<ZsStart[i]<<endl;
	}

	//## Get measurement matrix M
	TMatrixD M = TMatrixD(2*nHits,1,(double*)m," ");
	//cout<<"## Measurement Matrix M ##"<<endl;
	//M.Print();
	//cout<<endl;


	//## Calculate derivative matrix F & G
	//## Derivative of predictions f with respect to the four track parameter
	const int idim = 4;
	double d[2*nHits][idim];
	for(int i=0;i<nHits;i++) {
		int ix = i;
		int iy = i+nHits;
		for(int j=0;j<idim;j++) { d[ix][j] = 0; d[iy][j] = 0;}
		d[ix][0] = 1.;
		d[iy][1] = 1.;
		d[ix][2] = lenz[i];
		d[iy][3] = lenz[i];
	}
	TMatrixD F = TMatrixD(2*nHits,idim,(double*)d," ");
	cout<<"## Derivative Matrix F ##"<<endl;
	F.Print();
	cout<<endl;

	//## Calculate transpose of derivative matrix Ft
	TMatrixD Ft= TMatrixD(idim,2*nHits);
	Ft.Transpose(F);
	//cout<<"## Transpose Derivative Matrix Ft ##"<<endl;
	//Ft.Print();
	//cout<<endl;

	//## Calculate Covariance Matrix V
	double v[2*nHits][2*nHits];


	//## Path lenghts
	double Lr= 21.997;// radiation length of Malargue soil in g/cm^2 
	double Density= 1.8;// average density of Malargue soil in g/cm^3
	//double l= Ls[k-1];
	//double RadThickness= l*Density/Lr;
	//double EffectiveRadThickness= RadThickness*sqrt(1+tx*tx+ty*ty);

	//## Energy and Momentum
	double MuonMass= 105.65836668;//in MeV/c^2
	double TrueKinEnergy= fGenEnergy;
	double TrueMomentum= sqrt( TrueKinEnergy*(TrueKinEnergy+2.*MuonMass) );//in MeV/c

	double KinEnergy= 5000.;//in MeV
	double Momentum= sqrt( KinEnergy*(KinEnergy+2.*MuonMass) );//in MeV/c
	double AverageMomentum= sqrt( fAverageEnergyInMS*(fAverageEnergyInMS+2.*MuonMass) );//in MeV/c

	
	//double CMS= H*EffectiveRadThickness*pow( (1+0.038*log(EffectiveRadThickness)), 2);
	//double H= pow(13.6/Momentum,2);
	//double H= pow(13.6/Momentum,2)/Lr; //MISSING: PUT TRACK INCLINATION!!!
	
	double H;
	if(fUseAverageEnergyInMS){
		H= pow(13.6/Momentum,2)/Lr;
	}
	else if(fUseTrueEnergyInMS){
		H= pow(13.6/TrueMomentum,2)/Lr;
	}
	else if(fUseEnergyLossCorrectionInMS){
		cerr<<"TrackFinder::FastTrackFitterWithMS(): not yet implemented...exit"<<endl;
		exit(1);
	}
	else{
		cerr<<"TrackFinder::FastTrackFitterWithMS(): Specified fit with MS on but no choice given for MS treatment...exit"<<endl;
		exit(1);
	}


	for(int i=0;i<2*nHits;i++){
		for(int j=0;j<2*nHits;j++) v[i][j]= 0.;
	}
	
	//## diagonal elements
	for(int i=0;i<nHits;i++) {
		int ix = i;
		int iy = i+nHits;
		v[ix][ix]= pow(sigma[i][0],2);
		v[iy][iy]= pow(sigma[i][1],2);
		for(int s=0;s<nS;s++){
			double thisMSContribution= MSCovariance(hits[i][2],hits[i][2], Zs[s], ZsStart[s], Ls[s], H);
			v[ix][ix]+= thisMSContribution; 
			v[iy][iy]+= thisMSContribution;
			//cout<<"ix="<<ix<<"  iy="<<iy<<"  hits[i][2]="<<hits[i][2]<<"  Zs[="<<s<<"]="<<Zs[s]<<"  ZsStart[="<<s<<"]="<<ZsStart[s]<<"  thisMSContribution="<<thisMSContribution<<"  v="<<v[ix][ix]<<endl;
		}
		//cout<<"v["<<ix<<","<<ix<<"]= "<<v[ix][ix]<<endl;
		//cout<<"v["<<iy<<","<<iy<<"]= "<<v[iy][iy]<<endl;
	}

	//## non-diagonal elements
	v[0][1]= MSCovariance(hits[0][2],hits[1][2], Zs[0], ZsStart[0], Ls[0], H);
	v[1][0]= v[0][1];
	v[0+nHits][1+nHits]= v[0][1];
	v[1+nHits][0+nHits]= v[0+nHits][1+nHits];
	
	v[0][2]= MSCovariance(hits[0][2],hits[2][2], Zs[0], ZsStart[0], Ls[0], H);
	v[2][0]= v[0][2];
	v[0+nHits][2+nHits]= v[0][2];
	v[2+nHits][0+nHits]= v[0+nHits][2+nHits];
	
	v[1][2]= MSCovariance(hits[1][2],hits[2][2], Zs[0], ZsStart[0], Ls[0], H) + MSCovariance(hits[1][2],hits[2][2], Zs[1], ZsStart[1], Ls[1], H);
	v[2][1]= v[1][2];
	v[1+nHits][2+nHits]= v[1][2];
	v[2+nHits][1+nHits]= v[1+nHits][2+nHits];
	
	cout<<endl;
	for(int i=0;i<2*nHits;i++) {
		cout<<"v["<<i<<","<<i<<"]= "<<v[i][i]<<endl;
	}

	TMatrixD V = TMatrixD(2*nHits,2*nHits,(double*)v," ");
	//cout<<"## Covariance Matrix V ##"<<endl;
	//V.Print();
	//cout<<endl;

	//## Invert the Covariance Matrix
  double determ = 0.0;
  TMatrixD VInv = TMatrixD(2*nHits,2*nHits,(double*)v," ");
  VInv = VInv.Invert(&determ);
  if (determ<=0) {
		cerr<<"TrackFinder::FastTrackFitterWithMS(): Covariance Matrix inversion failed!"<<endl;
		//return -1;
	}
	//cout<<"## Inverted Covariance Matrix V ##"<<endl;
	//VInv.Print();
	//cout<<endl;


	//## Calculate Matrix (Ft VInv F)	= C covariance
	TMatrix FtVInvProduct= TMatrix(idim,2*nHits); 
	FtVInvProduct.Mult(Ft,VInv);
	TMatrix C= TMatrixD(idim,idim);
	C.Mult(FtVInvProduct,F);
	
	//## Invert the covariance C
	determ = 0.0;
  TMatrixD CInv = C;
  CInv = CInv.Invert(&determ);
  if (determ<=0) {
		cerr<<"TrackFinder::FastTrackFitterWithMS(): Final Covariance Matrix inversion failed!"<<endl;
		//return -1;
	}
	//cout<<"## Inverted Covariance Matrix C ##"<<endl;
	//CInv.Print();
	//cout<<endl;

	//## Calculate CInv Ft VInv M
	TMatrixD CInvFtProduct= TMatrixD(idim,2*nHits);
	CInvFtProduct.Mult(CInv,Ft);
	TMatrixD CInvFtVInvProduct= TMatrixD(idim,2*nHits);
	CInvFtVInvProduct.Mult(CInvFtProduct,VInv);

	TMatrixD FitParam= TMatrixD(idim,1);
	//FitParam.Mult(CInvFtVInvProduct,M);

	cout<<"## Fit Params ##"<<endl;
	FitParam.Print();
	cout<<endl;

	//## Find track solution
 	for (int k=0;k<nPar;k++) {
		FitParam(k,0)= par[k];
 	}
 	
	
	//## Calculate (M-F lambda)
	TMatrixD FParProduct= TMatrixD(2*nHits,1);
	FParProduct.Mult(F,FitParam);
	TMatrixD R= TMatrixD(2*nHits,1);
	R= M-FParProduct;
	TMatrixD Rt= TMatrixD(1,2*nHits);
	Rt.Transpose(R);


	//## Calculate Chi2
	TMatrixD RtVInvProduct= TMatrixD(1,2*nHits);	
	RtVInvProduct.Mult(Rt,VInv);
	TMatrixD Chi2Matrix= TMatrixD(1,1);
	Chi2Matrix.Mult(RtVInvProduct,R); 
	Chi2 = Chi2Matrix(0,0);
	

  //## Final result
  P0[0] = par[0];
  P0[1] = par[1];
  Phi = atan2(par[3],par[2]);
  Theta = acos(-sqrt(1-par[2]*par[2]-par[3]*par[3]));
  U[0] = sin(Theta)*cos(Phi);
  U[1] = sin(Theta)*sin(Phi);
  U[2] = cos(Theta);

	TVector3 x0(par[0], par[1], 0.); 
  TVector3 x1(par[0] + par[2], par[1] + par[3], 1.); 
  TVector3 direction = (x1-x0).Unit();
	
	//fTheta= Theta*180./TMath::Pi();
	//fPhi= Phi*180./TMath::Pi();
	fTheta= direction.Theta()*180./TMath::Pi();
	fPhi= direction.Phi()*180./TMath::Pi();

	fVertexX= x0.X();
	fVertexY= x0.Y();
	fTx= par[2];
	fTy= par[3];	

	//TVector3 oppositedirection = TVector3(-direction.X(),-direction.Y(),-direction.Z());
	//fTheta= oppositedirection.Theta()*180./TMath::Pi();
	//fPhi= oppositedirection.Phi()*180./TMath::Pi();

	fFitFCNMin= Chi2;

	cout<<"*** FAST TRACK RECO WITH MS ***"<<endl;
	cout<<"Theta [deg]="<<fTheta<<"  Phi [deg]="<<fPhi<<endl;
	cout<<"Theta [deg]="<<Theta<<"  Phi [deg]="<<Phi<<endl;
	cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 	cout<<"*******************************"<<endl;


	value = Chi2;

}//close TrackFinder::TrackChi2WithMSFcn()



bool TrackFinder::Checker(){

	if(fCandidateTrackPointCollection.empty()){
		cout<<"TrackFinder::Checker(): No track points given...exit"<<endl;
		return false;
	}

	// store seed det plane
	int lastSeedPlane=-1;
	std::vector<int> seedPlane;
	seedPlane.clear();
	seedPlane.resize(0);

	fTrackPointCollectionInPlane.clear();
	fTrackPointCollectionInPlane.resize(0);

	int entry_index= -1;

	for(unsigned int i=0;i<fCandidateTrackPointCollection.size();i++){
		int thisDetPlane= fCandidateTrackPointCollection[i].fDetectorPlaneId;
		
		if(thisDetPlane>lastSeedPlane){		
			seedPlane.push_back(thisDetPlane);
			fTrackPointCollectionInPlane.push_back( std::vector<TrackPoint>() );
			lastSeedPlane= thisDetPlane;
			entry_index++;
		}
		fTrackPointCollectionInPlane[entry_index].push_back(fCandidateTrackPointCollection[i]);

	}//close loop i


	if(fTrackPointCollectionInPlane.size()<3){
		cout<<"TrackFinder::Checker(): Cannot identify tracks from only "<<fTrackPointCollectionInPlane.size()<<" XY planes...exit"<<endl;
		return false;
	}


	if(fVerbosity>1)
	{
		cout<<"SEED PLANES"<<endl;
		for(unsigned int k=0;k<seedPlane.size();k++) cout<<"seedPlane["<<k<<"]="<<seedPlane[k]<<endl;

		cout<<"TRACK POINT IN PLANES"<<endl;
		for(unsigned int i=0;i<fTrackPointCollectionInPlane.size();i++) {
			for(unsigned int j=0;j<fTrackPointCollectionInPlane[i].size();j++) {
				cout<<"DET PLANE="<<fTrackPointCollectionInPlane[i][j].fDetectorPlaneId<<"  Pos=("<<(fTrackPointCollectionInPlane[i][j].fPosition).X()<<","<<(fTrackPointCollectionInPlane[i][j].fPosition).Y()<<","<<(fTrackPointCollectionInPlane[i][j].fPosition).Z()<<")"<<"  PosErr=("<<(fTrackPointCollectionInPlane[i][j].fPositionErr).X()<<","<<(fTrackPointCollectionInPlane[i][j].fPositionErr).Y()<<","<<(fTrackPointCollectionInPlane[i][j].fPositionErr).Z()<<")"<<endl;
			}
		}//close for i 
	}//close if verbosity

	return true;

}//close TrackFinder::Checker()


int TrackFinder::GetCombinationNo(){
        
	int n = 1;
  for(unsigned int i=0;i<fTrackPointCollectionInPlane.size();i++) n*= fTrackPointCollectionInPlane[i].size();
        
	return n;

}//close TrackFinder::GetCombinationNo()


void TrackFinder::CalculateCombinationMatrix(){
        
	int combNo = GetCombinationNo();
	int planeNo = (int)(fTrackPointCollectionInPlane.size());
  
	//clear existing combination matrix
	fCombinationMatrix.clear();
	fCombinationMatrix.resize(0);

	//resize new combination matrix
	for(int i=0;i<combNo;i++)	{
		fCombinationMatrix.push_back( std::vector<int>() );
		fCombinationMatrix[i].resize(planeNo);
	}      		

	//fill combination matrix
  int repetition = 1;
  for(int col= planeNo-1; col >= 0; col--) {
  	int max = (int)(fTrackPointCollectionInPlane[col].size());
    int value = 0;
    for(int row = 0; row < combNo;) {
    	for(int i = 0; i < repetition; i++, row++) {
      	fCombinationMatrix[row][col] = value;
      }
      value = value+1 >= max ? 0 : value+1;
    }
    repetition *= max;
  }
	
	//print combination matrix
	if(fVerbosity>1){
		cout<<"CombMatrix"<<endl;
		for(unsigned int i=0;i<fCombinationMatrix.size();i++)	{
			for(unsigned int j=0;j<fCombinationMatrix[i].size()-1;j++)	{
				cout<<fCombinationMatrix[i][j]<<"  ";
			}	
			cout<<fCombinationMatrix[i][fCombinationMatrix[i].size()-1]<<endl;
		}
	}	

}//close TrackFinder::CalculateCombinationMatrix()()

void TrackFinder::RemoveFakeHits(){

	//## Loop over track hits and removing those not due to the muon

	std::vector<TrackPoint>::iterator it= fTrackPointCollection.begin();
	cout<<"KFTrackFinder::RemoveFakeHits(): "<< fTrackPointCollection.size() <<"  total hits are present"<<endl;
	int nSpurious= 0;
	std::vector<int> pointToBeErased;
	pointToBeErased.clear();	
	pointToBeErased.resize(0);

	for(unsigned int i=0;i<fTrackPointCollection.size();i++){	
		int isMuonHit= fTrackPointCollection[i].fIsMuon;
		cout<<"point "<<i<<" ==>" <<isMuonHit<<endl;
		(fTrackPointCollection[i].fPosition).Print();
		if(isMuonHit==0){
			pointToBeErased.push_back(i);	
			nSpurious++;
		}
	}//end loop track points
		
	for(unsigned s=0;s<pointToBeErased.size();s++){
		//cout<<"pointToBeErased "<<pointToBeErased[s]<<endl;
		fTrackPointCollection.erase(it+pointToBeErased[s]);	
		it--;
	}


	cout<<endl;
	cout<<"KFTrackFinder::RemoveFakeHits(): "<< nSpurious <<"  fake hits removed, "<<fTrackPointCollection.size()<<"  left"<<endl;
	for(unsigned int i=0;i<fTrackPointCollection.size();i++){	
		int isMuonHit= fTrackPointCollection[i].fIsMuon;
		cout<<"point "<<i<<" ==>" <<isMuonHit<<endl;
		(fTrackPointCollection[i].fPosition).Print();
	}


}//close KFTrackFinder::RemoveFakeHits()

bool TrackFinder::FindCluster(){
	
	//## Search for clusters in each plane
	//## A cluster is defined by a group of hit points in adhiacent strips
	//## 1) Start from a track point in collection
	//## 2) Loop over all existing clusters (empty at init step) and search if current points
	//##    is adhiacent to any point of the looping cluster. If so, add the current point to the found cluster
	//##    and remove the point from the list
	//## 3) If cannot find any parent cluster, then form a cluster with this point and remove it from the list
	//## 4) Loop over all remaining track points in plane. If any is adhiacent to the current point, add to
	//##    its cluster. Then remove the added points from the list
	//## 5) Iterate from step 1) until all points have been assigned 
	//## 6) Do the same for the other planes
	
	if(fTrackPointCollection.empty()){
		cout<<"KFTrackFinder::ClusterFinder(): No track points...exit"<<endl;
		return false;
	}

	//## Store list of hit points per plane
	std::vector< std::vector<TrackPoint> > fTrackPointListInPlane;	
	fTrackPointListInPlane.clear();
	fTrackPointListInPlane.resize(0);

	int lastPlaneId= -1;
	int plane_index= -1;
		
	//cout<<"TrackPoint collection"<<endl;
	for(unsigned int i=0;i<fTrackPointCollection.size();i++){
		int thisDetPlane= fTrackPointCollection[i].fDetectorPlaneId;
		if(thisDetPlane>lastPlaneId){		
			fTrackPointListInPlane.push_back( std::vector<TrackPoint>() );
			lastPlaneId= thisDetPlane;
			plane_index++;
		}
		fTrackPointListInPlane[plane_index].push_back(fTrackPointCollection[i]);	
		//(fTrackPointCollection[i].fPosition).Print();
	}//end loop track points

	

	//cout<<"HitSeparationInCluster="<<fHitSeparationInCluster<<endl;
	fClusterTrackPointCollection.clear();
	fClusterTrackPointCollection.resize(0);
	

	//## Outer loop over detector planes
	for(unsigned int i=0;i<fTrackPointListInPlane.size();i++){
		//init cluster collection for this plane
		std::vector<Cluster> fClusterCollection;
		fClusterCollection.clear();
		fClusterCollection.resize(0);

		//cout<<"PLANE "<<i+1<<endl;

		bool IsClusterSearchFinished= false;
		while(!IsClusterSearchFinished){

			//## Loop over existing track points (some are going to be removed during this loop, when assigned)
			//for(unsigned int j=0;j<fTrackPointListInPlane[i].size();j++) {
			int j=0;
				TrackPoint thisPoint= fTrackPointListInPlane[i][j];
				double thisX= (thisPoint.fPosition).X();
				double thisY= (thisPoint.fPosition).Y();	
				double thisXAtCenter= (thisPoint.fPositionAtCenter).X();
				double thisYAtCenter= (thisPoint.fPositionAtCenter).Y();
				//cout<<"current point ("<<thisX<<","<<thisY<<")  @center=("<<thisXAtCenter<<","<<thisYAtCenter<<")"<<endl;

				std::vector<int> pointToBeErased;
				pointToBeErased.clear();	
				pointToBeErased.resize(0);
			
				//## Loop over existing clusters
				bool isAssignedPointInCluster= false;
				int cluster_index= -1;

				for(unsigned int k=0;k<fClusterCollection.size();k++) {
					Cluster thisCluster= fClusterCollection[k];
					std::vector<TrackPoint> thisClusterPointList= thisCluster.fClusterPoints;
	
					//## Loop over existing points in this cluster
					for(unsigned int l=0;l<thisClusterPointList.size();l++) {
						double thisClusterPointX= (thisClusterPointList[l].fPosition).X();
						double thisClusterPointY= (thisClusterPointList[l].fPosition).Y();
						double thisClusterPointXAtCenter= (thisClusterPointList[l].fPositionAtCenter).X();
						double thisClusterPointYAtCenter= (thisClusterPointList[l].fPositionAtCenter).Y();

						//Check if current point belongs to this cluster
						//double Dx= fabs(thisX-thisClusterPointX);
						//double Dy= fabs(thisY-thisClusterPointY);
						double Dx= fabs(thisXAtCenter-thisClusterPointXAtCenter);
						double Dy= fabs(thisYAtCenter-thisClusterPointYAtCenter);

						//cout<<"Dx="<<Dx<<"  Dy="<<Dy<<endl; 
						if(Dx<=fHitSeparationInCluster && Dy<=fHitSeparationInCluster){
							//add current point to this cluster
							fClusterCollection[k].AddPoint(thisPoint);
							isAssignedPointInCluster= true;
							cluster_index= k;
							//remove the current point from the list and exit
							//fTrackPointListInPlane[i].erase(fTrackPointListInPlane[i].begin()+j);	
							break;
						}
					}//end loop l points in current cluster

					if(isAssignedPointInCluster) break;
				}//end loop k existing clusters 
			
				//If the current point was not assigned
				//create a cluster and add this point
				//cout<<"isAssignedPointInCluster? "<<isAssignedPointInCluster<<"  cluster_index="<<cluster_index<<endl;
				Cluster newCluster;
				if(!isAssignedPointInCluster){
					newCluster.AddPoint(thisPoint);
					//remove the current point from the list
					//fTrackPointListInPlane[i].erase(fTrackPointListInPlane[i].begin()+j)
				}

				pointToBeErased.push_back(j);

				//loop over all following remaining points
				for(unsigned int s=j+1;s<fTrackPointListInPlane[i].size();s++) {
					TrackPoint nextPoint= fTrackPointListInPlane[i][s];
					double nextX= (nextPoint.fPosition).X();
					double nextY= (nextPoint.fPosition).Y();
					double nextXAtCenter= (nextPoint.fPositionAtCenter).X();
					double nextYAtCenter= (nextPoint.fPositionAtCenter).Y();
					//cout<<"next point ("<<nextX<<","<<nextY<<")  @center=("<<nextXAtCenter<<","<<nextYAtCenter<<")"<<endl;

					//double Dx= fabs(thisX-nextX);
					//double Dy= fabs(thisY-nextY);
					double Dx= fabs(thisXAtCenter-nextXAtCenter);
					double Dy= fabs(thisYAtCenter-nextYAtCenter);
					//cout<<"Dx="<<Dx<<"  Dy="<<Dy<<endl; 
					if(Dx<=fHitSeparationInCluster && Dy<=fHitSeparationInCluster){
						//add current point to the current point cluster
						if(!isAssignedPointInCluster) newCluster.AddPoint(nextPoint);
						else fClusterCollection[cluster_index].AddPoint(nextPoint);
						pointToBeErased.push_back(s);
					}
				}//end loop s remaining points

				
				std::vector<TrackPoint>::iterator it= fTrackPointListInPlane[i].begin();
				for(unsigned s=0;s<pointToBeErased.size();s++){
					//cout<<"pointToBeErased "<<pointToBeErased[s]<<endl;
					fTrackPointListInPlane[i].erase(it+pointToBeErased[s]);	
					it--;
				}
		
				if(!isAssignedPointInCluster) fClusterCollection.push_back(newCluster);
				//cout<<"fClusterCollection.size()="<<fClusterCollection.size()<<"  pointToBeErased.size()="<<pointToBeErased.size()<<endl;
	
				if(fTrackPointListInPlane[i].size()==0) IsClusterSearchFinished= true;
			//}	//end loop hit points per each plane

		}//end while cluster search

		//cout<<"CLUSTER LIST"<<endl;
		for(unsigned int k=0;k<fClusterCollection.size();k++) {
			Cluster thisCluster= fClusterCollection[k];
			std::vector<TrackPoint> thisClusterPointList= thisCluster.fClusterPoints;
			TrackPoint ClusterCentroid(0,0,0);
			//ClusterCentroid.fPosition= TVector3(0,0,0);
			//ClusterCentroid.fPositionErr= TVector3(0,0,0);
			int nPointsInCluster= (int)(thisClusterPointList.size());

			//## Split cluster?
			if(fSplitClustersInTracking && nPointsInCluster>fSplitClusterThreshold){
				//## Split cluster in single hits
				for(unsigned int l=0;l<thisClusterPointList.size();l++) {
					fClusterTrackPointCollection.push_back(thisClusterPointList[l]);
				}//end loop nPoints in current cluster
			}
			else{
				//## Merge hits in cluster, take baricenter				


			//cout<<"*** CLUSTER no "<<k<<endl; 
				for(unsigned int l=0;l<thisClusterPointList.size();l++) {
					double thisClusterPointX= (thisClusterPointList[l].fPosition).X();
					double thisClusterPointY= (thisClusterPointList[l].fPosition).Y();
					double thisClusterPointZ= (thisClusterPointList[l].fPosition).Z();
					double thisClusterPointErrX= (thisClusterPointList[l].fPositionErr).X();
					double thisClusterPointErrY= (thisClusterPointList[l].fPositionErr).Y();
					double thisClusterPointErrZ= (thisClusterPointList[l].fPositionErr).Z();
					double thisClusterPointTimeX= thisClusterPointList[l].fTimeX;
					double thisClusterPointTimeY= thisClusterPointList[l].fTimeY;
					double thisClusterPointEdepX= thisClusterPointList[l].fEdepX;
					double thisClusterPointEdepY= thisClusterPointList[l].fEdepY;
					//cout<<"P=("<<thisClusterPointX<<","<<thisClusterPointY<<","<<thisClusterPointZ<<")"<<endl;
					//cout<<"PErr=("<<thisClusterPointErrX<<","<<thisClusterPointErrY<<","<<thisClusterPointErrZ<<")"<<endl;
					ClusterCentroid.fPosition+= thisClusterPointList[l].fPosition;
					//ClusterCentroid.fPositionErr+= (thisClusterPointList[l].fPositionErr)*(thisClusterPointList[l].fPositionErr);
					ClusterCentroid.fPositionErr+= TVector3(thisClusterPointErrX*thisClusterPointErrX,thisClusterPointErrY*thisClusterPointErrY,thisClusterPointErrZ*thisClusterPointErrZ);
					ClusterCentroid.fEdepX+= thisClusterPointEdepX;
					ClusterCentroid.fEdepY+= thisClusterPointEdepY;
					if(thisClusterPointTimeX<ClusterCentroid.fTimeX) ClusterCentroid.fTimeX= thisClusterPointTimeX;
					if(thisClusterPointTimeY<ClusterCentroid.fTimeY) ClusterCentroid.fTimeY= thisClusterPointTimeY;
				}//end loop nPoints in current cluster

				//cout<<endl;

				if(nPointsInCluster!=0){
					ClusterCentroid.fDetectorPlaneId= thisClusterPointList[0].fDetectorPlaneId;
					ClusterCentroid.fPosition*=1./(double)(nPointsInCluster);
					ClusterCentroid.fPositionErr= TVector3(sqrt(ClusterCentroid.fPositionErr.X()),sqrt(ClusterCentroid.fPositionErr.Y()),sqrt(ClusterCentroid.fPositionErr.Z()));
					ClusterCentroid.fEdepX*=1./(double)(nPointsInCluster);
					ClusterCentroid.fEdepY*=1./(double)(nPointsInCluster);
					fClusterTrackPointCollection.push_back(ClusterCentroid);
				}
			}//end else split cluster
		
		}//end loop clusters
		cout<<endl;
	}//end loop planes

	
	//## Store hits and cluster collection info
	//cout<<"ClusterTrackPointCollection"<<endl;
	nClusterHit= (int)(fClusterTrackPointCollection.size());
	for(int i=0;i<nClusterHit;i++){
		//(fClusterTrackPointCollection[i].fPosition).Print();
		//(fClusterTrackPointCollection[i].fPositionErr).Print();
		ClusterHitX[i]= (fClusterTrackPointCollection[i].fPosition).X();
		ClusterHitY[i]= (fClusterTrackPointCollection[i].fPosition).Y();
		ClusterHitZ[i]= (fClusterTrackPointCollection[i].fPosition).Z();
	}
	
	nHit= (int)(fTrackPointCollection.size());
	for(int i=0;i<nHit;i++){
		HitX[i]= (fTrackPointCollection[i].fPosition).X();
		HitY[i]= (fTrackPointCollection[i].fPosition).Y();
		HitZ[i]= (fTrackPointCollection[i].fPosition).Z();
		HitEdepX[i]= fTrackPointCollection[i].fEdepX;
		HitEdepY[i]= fTrackPointCollection[i].fEdepY;
		HitTimeX[i]= fTrackPointCollection[i].fTimeX;
		HitTimeY[i]= fTrackPointCollection[i].fTimeY;
		IsMuonHit[i]= fTrackPointCollection[i].fIsMuon;
	}


	return true;
}//close TrackFinder::ClusterFinder()


void TrackFinder::SetTrackPointUncertainty(){

	//## Assign point uncertainty to track points
	TVector3 pointErr= TVector3(fHitSigmaX,fHitSigmaY,fHitSigmaZ);
	for(unsigned int i=0;i< fTrackPointCollection.size();i++){
		fTrackPointCollection[i].fPositionErr= pointErr;
		//fTrackPointCollection[i].fPositionErr.Print();
	}

}//close TrackFinder::SetTrackPointUncertainty()


double TrackFinder::MSCovariance(double Zi, double Zj, double Zs, double ZsStart, double L, double H){

	//## Zi, Zj= z hit positions
	//## Zs= z position of end of scatterer layer
	//## L= z length of scatterer layer
	//## H is the multiple scattering RMS angle per unit length
	//## Cms= (13.6 MeV/pc)^2 x t x [1 + 0.038 ln(t)]^2
	//## H= (13.6 MeV/pc)^2 / Lr, Lr radiation length 

	Zi= fabs(Zi);
	Zj= fabs(Zj);
	Zs= fabs(Zs);
	ZsStart= fabs(ZsStart);	
	

	if(Zi>Zj){
		cerr<<"TrackFinder::MSCovariance: Please specify Zi<=Zj!"<<endl;
		exit(1);
	}
	if(ZsStart>Zs){
		cerr<<"TrackFinder::MSCovariance: Please specify ZsStart<Zs!"<<endl;
		exit(1);
	}
	

	//if(Zj>Zs && Zi<ZsStart) return 0.;
	//if(Zi<Zs || Zj<Zs) return 0.;	

	double V= 0;

	if(Zi>=ZsStart && Zj<=Zs){
		V= 0.5*H* pow(Zi-ZsStart,2)* (Zj-ZsStart-(Zi-ZsStart)/3.);
	}
	else if(Zi<=Zs && Zi>=ZsStart && Zj>=Zs){
		V= 0.5*H* pow(Zi-ZsStart,2)* (L-(Zi-ZsStart)/3.+ (Zj-Zs));
	}
	else if(Zj<=Zs && Zj>=ZsStart && Zi>=Zs){
		V= 0.5*H* pow(Zj-ZsStart,2)* (L-(Zj-ZsStart)/3.+ (Zi-Zs));
	}
	else if(Zi>=Zs && Zj>=Zs){
		V= H*L* (pow(L,2)/3.+0.5*L*(Zi-Zs + Zj-Zs)+ (Zj-Zs)*(Zi-Zs));
	}
	else{
		return 0;
	}
	

	return V;

}//close TrackFinder::MSCovariance()


