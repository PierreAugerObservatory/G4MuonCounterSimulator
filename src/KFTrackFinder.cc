/**
* @file KFTrackFinder.cc
* @class KFTrackFinder
* @brief Reconstruct the particle tracks in the detector
*
* Reconstruct the particle tracks in the detector
* @author S. Riggi
* @date 07/09/2010
*/

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
#include <assert.h>

using namespace std;


bool KFTrackFinder::fIncludeMomentumInTracking;
bool KFTrackFinder::fTrackWithMSCorrection;
bool KFTrackFinder::fUseEnergyLossCorrectionInMS;
bool KFTrackFinder::fUseAverageEnergyInMS;
double KFTrackFinder::fAverageEnergyInMS;
bool KFTrackFinder::fUseTrueEnergyInMS;

double KFTrackFinder::fSuperPlaneSizeZ;
std::vector<double> KFTrackFinder::fSuperPlaneDepth;
int KFTrackFinder::fFitStatus;
double KFTrackFinder::fFitFCNMin;

TVector3 KFTrackFinder::fDirection;
double KFTrackFinder::fTheta;
double KFTrackFinder::fPhi;
double KFTrackFinder::fVertexX;
double KFTrackFinder::fVertexY;
double KFTrackFinder::fTx;
double KFTrackFinder::fTy;
double KFTrackFinder::fThetaStart;
double KFTrackFinder::fPhiStart;
double KFTrackFinder::fVertexXStart;
double KFTrackFinder::fVertexYStart;
double KFTrackFinder::fTxStart;
double KFTrackFinder::fTyStart;

double KFTrackFinder::fGenEnergy;
double KFTrackFinder::fGenTheta;
double KFTrackFinder::fGenPhi;
double KFTrackFinder::fGenVertexX;
double KFTrackFinder::fGenVertexY;
double KFTrackFinder::fGenTx;
double KFTrackFinder::fGenTy;

double KFTrackFinder::fThetaErr;
double KFTrackFinder::fPhiErr;
double KFTrackFinder::fVertexXErr;
double KFTrackFinder::fVertexYErr;
double KFTrackFinder::fTxErr;
double KFTrackFinder::fTyErr;
double KFTrackFinder::fThetaStartErr;
double KFTrackFinder::fPhiStartErr;
double KFTrackFinder::fVertexXStartErr;
double KFTrackFinder::fVertexYStartErr;
double KFTrackFinder::fTxStartErr;
double KFTrackFinder::fTyStartErr;

int KFTrackFinder::nHits;
int KFTrackFinder::nPar;
TMatrixD KFTrackFinder::fKFInitState;
TMatrixD KFTrackFinder::fKFInitCovariance;
std::vector<TMatrixD> KFTrackFinder::fKFState;
std::vector<TMatrixD> KFTrackFinder::fKFCovariance;
std::vector<TMatrixD> KFTrackFinder::fKFPredictedState;
std::vector<TMatrixD> KFTrackFinder::fKFPredictedCovariance;
std::vector<TMatrixD> KFTrackFinder::fKFMeas;	
std::vector<TMatrixD> KFTrackFinder::fKFMeasVariance;
std::vector<TMatrixD> KFTrackFinder::fKFPropagator;
std::vector<TMatrixD> KFTrackFinder::fKFProjector;
std::vector<TMatrixD> KFTrackFinder::fKFCovarianceMS;
std::vector<TMatrixD> KFTrackFinder::fKFMeasResidual;
std::vector<TMatrixD> KFTrackFinder::fKFCovarianceResidual;
std::vector<TMatrixD> KFTrackFinder::fKFPredictedMeasResidual;
std::vector<TMatrixD> KFTrackFinder::fKFPredictedCovarianceResidual;
std::vector<TMatrixD> KFTrackFinder::fKFSmoothedState;
std::vector<TMatrixD> KFTrackFinder::fKFSmoothedCovariance;
std::vector<TMatrixD> KFTrackFinder::fKFSmoothedMeasResidual;
std::vector<TMatrixD> KFTrackFinder::fKFSmoothedCovarianceResidual;

std::vector<double> KFTrackFinder::fKFHitChi2;
std::vector<double> KFTrackFinder::fKFHitPredictedChi2;
double KFTrackFinder::fMomentum;
double KFTrackFinder::fKinEnergy;
double KFTrackFinder::lenz[MAXHITNO];
double KFTrackFinder::TrackHitChi2[MAXHITNO];
double KFTrackFinder::TrackHitPredictedChi2[MAXHITNO];
double KFTrackFinder::Ls[MAXSCATNO];
double KFTrackFinder::Zs[MAXSCATNO];
double KFTrackFinder::ZsStart[MAXSCATNO];


KFTrackFinder::KFTrackFinder(){

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
	fDirection.SetXYZ(-999,-999,-999);
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
	fMaxAllowedAngle= 10.;//in degrees
}

KFTrackFinder::~KFTrackFinder(){

}


bool KFTrackFinder::TrackReconstructor(){

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

}//close KFTrackFinder::FindTrack()


bool KFTrackFinder::FindTrack(){

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
			cout<<"KFTrackFinder::FindTrack(): End track search"<<endl;	
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

		std::vector<TVector3> DirectionFit;
		DirectionFit.clear();
		DirectionFit.resize(0);

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

		std::vector<double> ThetaStartFit;
		ThetaStartFit.clear();
		ThetaStartFit.resize(0);

		std::vector<double> PhiStartFit;
		PhiStartFit.clear();
		PhiStartFit.resize(0);

		std::vector<double> VertexXStartFit;
		VertexXStartFit.clear();
		VertexXStartFit.resize(0);

		std::vector<double> VertexYStartFit;
		VertexYStartFit.clear();
		VertexYStartFit.resize(0);

		std::vector<double> TxStartFit;
		TxStartFit.clear();
		TxStartFit.resize(0);

		std::vector<double> TyStartFit;
		TyStartFit.clear();
		TyStartFit.resize(0);

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

		std::vector<double> ThetaStartFitErr;
		ThetaStartFitErr.clear();
		ThetaStartFitErr.resize(0);

		std::vector<double> PhiStartFitErr;
		PhiStartFitErr.clear();
		PhiStartFitErr.resize(0);

		std::vector<double> VertexXStartFitErr;
		VertexXStartFitErr.clear();
		VertexXStartFitErr.resize(0);

		std::vector<double> VertexYStartFitErr;
		VertexYStartFitErr.clear();
		VertexYStartFitErr.resize(0);

		std::vector<double> TxStartFitErr;
		TxStartFitErr.clear();
		TxStartFitErr.resize(0);

		std::vector<double> TyStartFitErr;
		TyStartFitErr.clear();
		TyStartFitErr.resize(0);

  	int trackCounter=0;	

		std::vector< std::vector<double> > TrackFilterChi2;
		std::vector< std::vector<double> > TrackPredictedChi2;
		//double TrackFilterChi2[MAXTRACKNO][MAXHITNO];
		//double TrackPredictedChi2[MAXTRACKNO][MAXHITNO];

		//## Loop and fit all possible combinations
		//std::vector< std::vector<int> >::iterator combMatrixIt= fCombinationMatrix.begin();

		for(int i = 0; i < fCombinationMatrix.size(); i++) {
			int rowIndex = 0;
   		int colIndex = 0;
			std::vector<int> row (fCombinationMatrix[i].begin(),fCombinationMatrix[i].end());
				
			Track* aCandidateTrack= new Track;
			
    	//while(colIndex < fTrackPointCollectionInPlane.size()) {
			for(unsigned int j=0;j<fTrackPointCollectionInPlane.size(); j++){
				//TrackPoint aPoint= fTrackPointCollectionInPlane[rowIndex++][row[colIndex++]];
				TrackPoint aPoint= fTrackPointCollectionInPlane[j][fCombinationMatrix[i][j]];
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

			/*
			if(!AcceptCandidateTrack(aCandidateTrack)){
				delete aCandidateTrack;
				//fCombinationMatrix.erase(combMatrixIt+i);
				//combMatrixIt--;
				continue;
			}
			*/

			CandidateTrackList.push_back(aCandidateTrack);
	
			cout<<"*** Fitting track no. "<< trackCounter << " ***"<<endl;
			//## track fitter
			//fFitStatus= LSTrackFitter(aCandidateTrack);
			fFitStatus= KFTrackFitter(aCandidateTrack);
			//fFitStatus= KFTrackFitterWithMomentum(aCandidateTrack);

			//cout<<"Track "<<trackCounter<<"  FitFCN="<<fFitFCNMin<<"  FitStatus="<<fFitStatus<<endl;
			TrackFitStatus.push_back(fFitStatus);
			TrackFitFCNMin.push_back(fFitFCNMin);
			DirectionFit.push_back(fDirection);	
			ThetaFit.push_back(fTheta);	
			PhiFit.push_back(fPhi);
			VertexXFit.push_back(fVertexX);	
			VertexYFit.push_back(fVertexY);	
			TxFit.push_back(fTx);
			TyFit.push_back(fTy);

			ThetaStartFit.push_back(fThetaStart);	
			PhiStartFit.push_back(fPhiStart);
			VertexXStartFit.push_back(fVertexXStart);	
			VertexYStartFit.push_back(fVertexYStart);	
			TxStartFit.push_back(fTxStart);
			TyStartFit.push_back(fTyStart);

			ThetaFitErr.push_back(fThetaErr);	
			PhiFitErr.push_back(fPhiErr);
			VertexXFitErr.push_back(fVertexXErr);	
			VertexYFitErr.push_back(fVertexYErr);	
			TxFitErr.push_back(fTxErr);
			TyFitErr.push_back(fTyErr);

			ThetaStartFitErr.push_back(fThetaStartErr);	
			PhiStartFitErr.push_back(fPhiStartErr);
			VertexXStartFitErr.push_back(fVertexXStartErr);	
			VertexYStartFitErr.push_back(fVertexYStartErr);	
			TxStartFitErr.push_back(fTxStartErr);
			TyStartFitErr.push_back(fTyStartErr);
	
			TrackFilterChi2.push_back ( std::vector<double>() );
			TrackPredictedChi2.push_back ( std::vector<double>() );
			
			for(unsigned int ll=0;ll<fKFHitChi2.size();ll++){
				//TrackFilterChi2[trackCounter][ll]= fKFHitChi2[ll];
				//TrackPredictedChi2[trackCounter][ll]= fKFHitPredictedChi2[ll];

				TrackFilterChi2[trackCounter].push_back(fKFHitChi2[ll]);
				TrackPredictedChi2[trackCounter].push_back(fKFHitPredictedChi2[ll]);
			}
			
			trackCounter++;
			
    	cout<<endl;
		}//close for combination matrix

		//cout<<"*** search best fit among candidate tracks ***"<<endl;
		//## search best fit among candidate tracks
		double min= 1.e+99;
		int bestTrackIndex= -1;
		for(unsigned int k=0;k<TrackFitFCNMin.size();k++){
			//cout<<"Track "<<k<<"  FitFCN="<<TrackFitFCNMin[k]<<"  FitStatus="<<TrackFitStatus[k]<<endl;
			if(TrackFitFCNMin[k]<min && TrackFitStatus[k]==0) {
				min= TrackFitFCNMin[k];
				bestTrackIndex= k;
			}
		}

		//cout<<"bestTrackIndex="<<bestTrackIndex<<endl;

		if(bestTrackIndex!=-1){
			//## fill the list of reconstructed tracks
			CandidateTrackList[bestTrackIndex]->SetId(rectrackCounter);
			CandidateTrackList[bestTrackIndex]->SetTheta(ThetaFit[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetPhi(PhiFit[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetFitStatus(TrackFitStatus[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetFitChi2(TrackFitFCNMin[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetDirection(DirectionFit[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetVertexPos(TVector3(VertexXFit[bestTrackIndex],VertexYFit[bestTrackIndex],0.));
			CandidateTrackList[bestTrackIndex]->SetTx(TxFit[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetTy(TyFit[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetThetaStart(ThetaStartFit[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetPhiStart(PhiStartFit[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetVertexPosStart(TVector3(VertexXStartFit[bestTrackIndex],VertexYStartFit[bestTrackIndex],0.));
			CandidateTrackList[bestTrackIndex]->SetTxStart(TxStartFit[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetTyStart(TyStartFit[bestTrackIndex]);

			for(unsigned int ll=0;ll<fKFHitChi2.size();ll++) {
				CandidateTrackList[bestTrackIndex]->AddTrackHitChi2(TrackFilterChi2[bestTrackIndex][ll]);
				CandidateTrackList[bestTrackIndex]->AddTrackHitExpChi2(TrackPredictedChi2[bestTrackIndex][ll]);	
			}

			CandidateTrackList[bestTrackIndex]->SetThetaErr(ThetaFitErr[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetPhiErr(PhiFitErr[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetVertexPosErr(TVector3(VertexXFitErr[bestTrackIndex],VertexYFitErr[bestTrackIndex],0.));
			CandidateTrackList[bestTrackIndex]->SetTxErr(TxFitErr[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetTyErr(TyFitErr[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetThetaStartErr(ThetaStartFitErr[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetPhiStartErr(PhiStartFitErr[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetVertexPosStartErr(TVector3(VertexXStartFitErr[bestTrackIndex],VertexYStartFitErr[bestTrackIndex],0.));
			CandidateTrackList[bestTrackIndex]->SetTxStartErr(TxStartFitErr[bestTrackIndex]);
			CandidateTrackList[bestTrackIndex]->SetTyStartErr(TyStartFitErr[bestTrackIndex]);

			fRecTrackCollection.push_back(CandidateTrackList[bestTrackIndex]);

			rectrackCounter++;

			//## remove assigned track points from the list of candidate points
			//get track points of reconstructed track
			std::vector<int> RecTrackPointDetIdList= CandidateTrackList[bestTrackIndex]->GetPointDetId();
			double x,y,z,t;
			double xErr,yErr,zErr,tErr;
			std::vector<TrackPoint> RecTrackPointList;
			//TrackPoint aRecTrackPoint;
			for(int k=0;k<CandidateTrackList[bestTrackIndex]->GetNpoints();k++){
				CandidateTrackList[bestTrackIndex]->GetPoint(k,x,y,z,t);
				CandidateTrackList[bestTrackIndex]->GetPointError(k,xErr,yErr,zErr,tErr);
				//CandidateTrackList[bestTrackIndex]->GetTrackPoint(k,aRecTrackPoint);	
	
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

			//cout<<"Erase points"<<endl;
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
		}//close if bestTrackIndex
		else IsTrackSearchFinished= true;


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

}//close KFTrackFinder::KFTrackFinder()


bool KFTrackFinder::AcceptCandidateTrack(Track* aCandidateTrack){

	//## Accept a candidate track
	//## Compare relative angles in space of tracks formed with each point
	double x,y,z,t;
	int nPoints= aCandidateTrack->GetNpoints();
	TrackPoint thisTrackPoints[nPoints];
	TVector3 thisTrackPointVector[nPoints];
	TVector3 thisTrackVector[nPoints];
	double RelAngleAmongTracks[nPoints];	

	for(int p=0;p<nPoints;p++){
		aCandidateTrack->GetTrackPoint(p,thisTrackPoints[p]);
		thisTrackPointVector[p]= thisTrackPoints[p].fPosition;
	}

	for(int p=0;p<nPoints-1;p++){
		thisTrackVector[p]= thisTrackPointVector[p+1]-thisTrackPointVector[p];
	}
	thisTrackVector[nPoints-1]= thisTrackPointVector[nPoints-1]-thisTrackPointVector[0];

	//## Calculate relative angle in space
	for(int p=0;p<nPoints-1;p++){
		RelAngleAmongTracks[p]= thisTrackVector[p].Angle(thisTrackVector[p+1]) * 180./TMath::Pi();
		//cout<<"RelAngle "<<p<<"-"<<p+1<<" ="<< RelAngleAmongTracks[p]<<endl;
	}
	RelAngleAmongTracks[nPoints-1]= thisTrackVector[nPoints-1].Angle(thisTrackVector[0]) * 180./TMath::Pi();
	//cout<<"RelAngle 0-"<<nPoints-1<<" ="<< RelAngleAmongTracks[nPoints-1]<<endl;
	
	//## Accept candidate track?
	for(int p=0;p<nPoints-1;p++){
		if(RelAngleAmongTracks[p]>fMaxAllowedAngle)
			return false;
	}
	

	return true;

}// close KFTrackFinder::AcceptCandidateTrack()


int KFTrackFinder::KFTrackFitter(Track* aTrack){

	//cout<<"*********************************"<<endl;
	//cout<<"***     KF TRACK FITTER       ***"<<endl;
	//cout<<"*********************************"<<endl;

	//## Consistency checks
	nHits= aTrack->GetNpoints();
	if (nHits<3) return -1;

	int nIter= nHits+1;
	nPar= 4;
	int nMeas= 2;

	//## Init vectors
	fKFState.clear();
	fKFState.resize(0);
	//fKFState.resize(nIter);

	fKFCovariance.clear();
	fKFCovariance.resize(0);
	//fKFCovariance.resize(nIter);

	fKFSmoothedState.clear();
	fKFSmoothedState.resize(0);
	//fKFSmoothedState.resize(nIter);

	fKFSmoothedCovariance.clear();
	fKFSmoothedCovariance.resize(0);
	//fKFSmoothedCovariance.resize(nIter);

	fKFPredictedState.clear();
	fKFPredictedState.resize(0);
	//fKFPredictedState.resize(nHits);
		
	fKFPredictedCovariance.clear();
	fKFPredictedCovariance.resize(0);
	//fKFPredictedCovariance.resize(nHits);

	fKFProjector.clear();
	fKFProjector.resize(0);
	//fKFProjector.resize(nHits);

	fKFMeas.clear();
	//fKFMeas.resize(nHits);
	fKFMeas.resize(0);

	fKFMeasVariance.clear();
	//fKFMeasVariance.resize(nHits);
	fKFMeasVariance.resize(0);

	fKFPropagator.clear();
	//fKFPropagator.resize(nHits);
	fKFPropagator.resize(0);
	
	fKFCovarianceMS.clear();
	fKFCovarianceMS.resize(0);
	//fKFCovarianceMS.resize(nHits);

	fKFCovarianceResidual.clear();
	fKFCovarianceResidual.resize(0);
	//fKFCovarianceResidual.resize(nHits);

	fKFMeasResidual.clear();
	fKFMeasResidual.resize(0);
	//fKFMeasResidual.resize(nHits);

	fKFSmoothedCovarianceResidual.clear();
	fKFSmoothedCovarianceResidual.resize(0);
	//fKFSmoothedCovarianceResidual.resize(nHits);

	fKFSmoothedMeasResidual.clear();
	fKFSmoothedMeasResidual.resize(0);
	//fKFSmoothedMeasResidual.resize(nHits);

	fKFPredictedCovarianceResidual.clear();
	fKFPredictedCovarianceResidual.resize(0);
	//fKFPredictedCovarianceResidual.resize(nHits);

	fKFPredictedMeasResidual.clear();
	fKFPredictedMeasResidual.resize(0);
	//fKFPredictedMeasResidual.resize(nHits);

	fKFHitChi2.clear();
	fKFHitChi2.resize(0);

	fKFHitPredictedChi2.clear();
	fKFHitPredictedChi2.resize(0);
	
	for(int i=0;i<nHits;i++) {
		TrackHitChi2[i]= 0.;
		TrackHitPredictedChi2[i]= 0.;
		fKFHitChi2.push_back(0.);
		fKFHitPredictedChi2.push_back(0.);
	}
	
	//## Choose initial state vector and covariance
	//int LSStatus= LSTrackFitter(aTrack);
	//if(LSStatus!=0){
	//	cerr<<"KFTrackFinder::KalmanFilterTrackFitter(): WARNING: LSFitter (used to init the KFFitter) failed!"<<endl;
	//}
	
	/*
	cout<<"*** Initial KFState x0 ***"<<endl;
	fKFInitState.Print();
	cout<<"*** Initial KFCovariance C0 ***"<<endl;
	fKFInitCovariance.Print(); 
	*/


	TMatrixD x0State(nPar,1);
	TMatrixD C0(nPar,nPar);
	TMatrixD dummyStateMatrix(nPar,1);
	TMatrixD dummyCovarianceMatrix(nPar,nPar);

	for(int i=0;i<nPar;i++){
		//x0State(i,0)= gRandom->Uniform(-100,100);
		x0State(i,0)= 0.;
		dummyStateMatrix(i,0)= 0.;
		for(int j=0;j<nPar;j++){
			C0(i,j)= 0.;	
			//C0(i,j)= gRandom->Uniform(1,100);	
			dummyCovarianceMatrix(i,j)= 0.;
		}
		C0(i,i)= 1000.;		
	}

	//cout<<"## Initial KSStateVector x0 ##"<<endl;
	//x0State.Print();
	//cout<<"## Initial KSCovariance C0 ##"<<endl;
	//C0.Print(); 
	
	
	TMatrixD dummyResMatrix(nMeas,nMeas);
	TMatrixD dummyResMatrix2(nMeas,1);
	
	for(int i=0;i<nMeas;i++){
		dummyResMatrix2(i,0)= 0.;
		for(int j=0;j<nMeas;j++){
			dummyResMatrix(i,j)= 0.;	
		}	
	}

	
	//## Init state and covariance matrix vectors
	for(int i=0;i<nIter;i++){
		if(i==0){
			fKFState.push_back(x0State);
			fKFCovariance.push_back(C0);
			//fKFState.push_back(fKFInitState);
			//fKFCovariance.push_back(fKFInitCovariance);
		}	
		else{
			fKFState.push_back(dummyStateMatrix);
			fKFCovariance.push_back(dummyCovarianceMatrix);
		}
		fKFSmoothedState.push_back(dummyStateMatrix);
		fKFSmoothedCovariance.push_back(dummyCovarianceMatrix);
	}

	for(int i=0;i<nHits;i++){
		fKFPredictedState.push_back(dummyStateMatrix);
		fKFPredictedCovariance.push_back(dummyCovarianceMatrix);
		fKFCovarianceMS.push_back(dummyCovarianceMatrix);

		fKFCovarianceResidual.push_back(dummyResMatrix);
		fKFMeasResidual.push_back(dummyResMatrix2);

		fKFSmoothedCovarianceResidual.push_back(dummyResMatrix);
		fKFSmoothedMeasResidual.push_back(dummyResMatrix2);

		fKFPredictedCovarianceResidual.push_back(dummyResMatrix);
		fKFPredictedMeasResidual.push_back(dummyResMatrix2);
	}


	//## Fill measurement, propagator and projector matrix
	//## Measurement and Variance Matrix
	double hits[nHits][3];
	double sigma[nHits][3];
	//double lenz[nHits];
	double x,y,z,time;
	double xErr,yErr,zErr,timeErr;
	TMatrixD M(nMeas,1);
	TMatrixD V(nMeas,nMeas);
	for(int i=0; i<nHits; i++) {
  	aTrack->GetPoint(i,x,y,z,time);	
		aTrack->GetPointError(i,xErr,yErr,zErr,timeErr);	
		M(0,0)= x;
		M(1,0)= y;

		//V(0,0)= xErr;
		//V(1,1)= yErr;
		V(0,0)= xErr*xErr;
		V(1,1)= yErr*yErr;
		V(0,1)= 0.;
		V(1,0)= 0.;

		//fKFMeas[i]= M;
		//fKFMeasVariance[i]= V;
		fKFMeas.push_back(M);
		fKFMeasVariance.push_back(V);

		hits[i][0] = x;
		hits[i][1] = y;
		hits[i][2] = z;
			
		sigma[i][0] = xErr;
		sigma[i][1] = yErr;
		sigma[i][2] = zErr;

		lenz[i] = hits[i][2];

		/*
		cout<<"X["<<i<<"]="<< hits[i][0]<<endl;
		cout<<"Y["<<i<<"]="<< hits[i][1]<<endl;
		cout<<"Z["<<i<<"]="<< hits[i][2]<<endl;
		cout<<"Sigma["<<i<<"]="<< sigma[i][0]<<"  "<<sigma[i][1]<<endl;
		*/
	}	

	//## Lenghts Z(k)-Z(k-1)
	for(int i=1;i<nHits;i++) lenz[i] = hits[i][2] - hits[i-1][2];

	//## Scatterer positions
	const int nS= 3;
	Ls[0]= fSuperPlaneDepth[0]-fSuperPlaneSizeZ/2.;
	Ls[1]= fSuperPlaneDepth[1]-fSuperPlaneDepth[0]-fSuperPlaneSizeZ;
	Ls[2]= fSuperPlaneDepth[2]-fSuperPlaneDepth[1]-fSuperPlaneSizeZ;
	ZsStart[0]= 0.;
	ZsStart[1]= -(fSuperPlaneDepth[0]+fSuperPlaneSizeZ/2.);
	ZsStart[2]= -(fSuperPlaneDepth[1]+fSuperPlaneSizeZ/2.);
	for(int i=0;i<nS;i++){
		Zs[i]= -(fSuperPlaneDepth[i]-fSuperPlaneSizeZ/2.);
		//cout<<"Ls["<<i<<"]="<<Ls[i]<<"  Zs["<<i<<"]="<<Zs[i]<<"  ZsStart["<<i<<"]="<<ZsStart[i]<<endl;
	}

	//## Propagator Matrix
	TMatrixD F(nPar,nPar);
	for(int i=0;i<nHits;i++){
		for(int j=0;j<nPar;j++) { 
			for(int k=0;k<nPar;k++) {F(j,k)= 0.;}
			F(j,j)= 1.;
		}
		F(0,1)= 0.;
		F(0,2)= lenz[i];
		F(0,3)= 0.;

		F(1,0)= 0.;
		F(1,2)= 0.;
		F(1,3)= lenz[i];
		
		//fKFPropagator[i]= F;
		fKFPropagator.push_back(F);
	}	

	//## Projector Matrix
	TMatrixD H(nMeas,nPar);
	for(int i=0;i<nHits;i++){
		for(int j=0;j<nMeas;j++){
			for(int k=0;k<nPar;k++) {H(j,k)= 0.;}
		}
		H(0,0)= 1.;
		H(1,1)= 1.;
		//fKFProjector[i]= H;
		fKFProjector.push_back(H);
	}


	//## Printing
	//DebugPrint();

	fFitFCNMin= 0;
	
	double par[nPar];
	double parErr[nPar];

	
	//## Iterate doing prediction and filter steps
	int nKFIterations= 20;
	for(int s=0;s<nKFIterations;s++){
		fFitFCNMin= 0;
		fKinEnergy= 5000.;

		if(s>0){
			//## Re-initialize the state and covariance to the smoothed vector of previous iteration
			fKFState[0]= fKFSmoothedState[0];
			fKFCovariance[0]= fKFSmoothedCovariance[0];
		}

		//cout<<"*******************************"<< endl;
		//cout<<"***      KF Iteration --> "<< s << "  ***"<< endl;
		//cout<<"*******************************"<< endl;
		for(int k=1;k<nIter;k++){
			KFPredictionStep(k);
			KFFilterStep(k);
		}
		for(int k=nHits;k>=0;k--){
			KFSmoothingStep(k);
		}

		TMatrixD thisIterationTrackState= fKFSmoothedState[0];
		par[0]= thisIterationTrackState(0,0);
		par[1]= thisIterationTrackState(1,0);
		par[2]= thisIterationTrackState(2,0);
		par[3]= thisIterationTrackState(3,0);

		TMatrixD thisIterationCovariance= fKFSmoothedCovariance[0];
		parErr[0]= sqrt(thisIterationCovariance(0,0));
		parErr[1]= sqrt(thisIterationCovariance(1,1));
		parErr[2]= sqrt(thisIterationCovariance(2,2));
		parErr[3]= sqrt(thisIterationCovariance(3,3));

		TVector3 x0(par[0], par[1], 0.); 
  	TVector3 x1(par[0] + par[2], par[1] + par[3], 1.); 
  	TVector3 direction = (x1-x0).Unit();
		TVector3 oppositedirection = TVector3(-direction.X(),-direction.Y(),-direction.Z());
		TVector3 oppositeZdirection = TVector3(direction.X(),direction.Y(),-direction.Z());
	
		//## track parameters
		fTheta= direction.Theta()*180./TMath::Pi();
		fPhi= direction.Phi()*180./TMath::Pi();

		fVertexX= x0.X();
		fVertexY= x0.Y();	
		fTx= par[2];
		fTy= par[3];

		//## track parameter errors
		fVertexXErr= parErr[0];
		fVertexYErr= parErr[1];
		fTxErr= parErr[2];
		fTyErr= parErr[3];

		//## theta & phi par error
		//## theta= atan2(sqrt(tx^2+ty^2),1)
		//## phi= atan2(ty,tx)
		//## do the error propagation as FCFt, with F derivative matrix d(theta)/d(par_i), d(phi)/d(par_i) 
		TMatrixD DTheta(1,nPar);	
		DTheta(0,0)= 0;
		DTheta(0,1)= 0;	
		DTheta(0,2)= fTx/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
 		DTheta(0,3)= fTy/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
		TMatrixD DThetaT(nPar,1);		
		DThetaT.Transpose(DTheta);

		TMatrixD DPhi(1,nPar);	
		DPhi(0,0)= 0;
		DPhi(0,1)= 0;		
		DPhi(0,2)= -fTy/(fTx*fTx+fTy*fTy);
 		DPhi(0,3)= fTx/(fTx*fTx+fTy*fTy);
		TMatrixD DPhiT(nPar,1);		
		DPhiT.Transpose(DPhi);

		TMatrixD DThetaMultCov(1,nPar);
		DThetaMultCov.Mult(DTheta,thisIterationCovariance);
		TMatrixD DPhiMultCov(1,nPar);
		DPhiMultCov.Mult(DPhi,thisIterationCovariance);
		
		TMatrixD ThetaUncertainty(1,1);
		ThetaUncertainty.Mult(DThetaMultCov,DThetaT);

		TMatrixD PhiUncertainty(1,1);
		PhiUncertainty.Mult(DPhiMultCov,DPhiT);

		fThetaErr= sqrt(ThetaUncertainty(0,0))* 180./TMath::Pi();
		fPhiErr= sqrt(PhiUncertainty(0,0))* 180./TMath::Pi();

		if(s==0){
			//## Store first iter results
			fThetaStart= fTheta;
			fPhiStart= fPhi;
			fTxStart= fTx;
			fTyStart= fTy;
			fVertexXStart= fVertexX;
			fVertexYStart= fVertexY;

			fVertexXStartErr= fVertexXErr;
			fVertexYStartErr= fVertexYErr;
			fTxStartErr= fTxErr;
			fTyStartErr= fTyErr;
			fThetaStartErr= fThetaErr;		
			fPhiStartErr= fPhiErr;
		}

		//cout<<"*** KALMAN FILTER TRACK RECO (ITER "<<s<<") ***"<<endl;
		//cout<<"Theta [deg]="<<fTheta<<" +- "<< fThetaErr << "  Phi [deg]="<<fPhi<<" +- "<<fPhiErr<<endl;
		//cout<<"X0 [cm]="<< fVertexX <<" +- "<< fVertexXErr << "  Y0[cm]="<<fVertexY<<" +- "<< fVertexYErr<<endl;
		//cout<<"Tx="<<fTx<<" +- "<< fTxErr <<"  Ty="<<fTy<<" +- "<< fTyErr <<endl;
		//cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 		//cout<<"*******************************"<<endl;
	}//end iteration
	
	cout<<endl;
	
	
	/*
	KFPredictionStep(1);
	KFFilterStep(1);

	KFPredictionStep(2);
	KFFilterStep(2);

	KFPredictionStep(3);
	KFFilterStep(3);

	KFSmoothingStep(3);
	KFSmoothingStep(2);
	KFSmoothingStep(1);
	KFSmoothingStep(0);
	*/

	//## Printing
	//DebugPrint();

	
	//TMatrixD finalTrackState= fKFState[nIter-1];
	TMatrixD finalTrackState= fKFSmoothedState[0];
	par[0]= finalTrackState(0,0);
	par[1]= finalTrackState(1,0);
	par[2]= finalTrackState(2,0);
	par[3]= finalTrackState(3,0);

	TMatrixD finalTrackCovariance= fKFSmoothedCovariance[0];
	parErr[0]= sqrt(finalTrackCovariance(0,0));
	parErr[1]= sqrt(finalTrackCovariance(1,1));
	parErr[2]= sqrt(finalTrackCovariance(2,2));
	parErr[3]= sqrt(finalTrackCovariance(3,3));

	TVector3 x0(par[0], par[1], 0.); 
  TVector3 x1(par[0] + par[2], par[1] + par[3], 1.); 
  TVector3 direction = (x1-x0).Unit();
	TVector3 oppositedirection = TVector3(-direction.X(),-direction.Y(),-direction.Z());
	TVector3 oppositeZdirection = TVector3(direction.X(),direction.Y(),-direction.Z());
	//fTheta= oppositedirection.Theta()*180./TMath::Pi();
	//fPhi= oppositedirection.Phi()*180./TMath::Pi();

	//fTheta= Theta*180./TMath::Pi();
	//fPhi= Phi*180./TMath::Pi();
	fTheta= direction.Theta()*180./TMath::Pi();
	fPhi= direction.Phi()*180./TMath::Pi();

	fDirection= direction;

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
	TMatrixD DThetaFinal(1,nPar);	
	DThetaFinal(0,0)= 0;
	DThetaFinal(0,1)= 0;		
	DThetaFinal(0,2)= fTx/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
 	DThetaFinal(0,3)= fTy/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
	TMatrixD DThetaFinalT(nPar,1);		
	DThetaFinalT.Transpose(DThetaFinal);
	//cout<<"*** Theta Derivative ***"<<endl;
	//DThetaFinal.Print();
	//cout<<endl;

	TMatrixD DPhiFinal(1,nPar);	
	DPhiFinal(0,0)= 0;
	DPhiFinal(0,1)= 0;		
	DPhiFinal(0,2)= -fTy/(fTx*fTx+fTy*fTy);
 	DPhiFinal(0,3)= fTx/(fTx*fTx+fTy*fTy);
	TMatrixD DPhiFinalT(nPar,1);		
	DPhiFinalT.Transpose(DPhiFinal);
	//cout<<"*** Phi Derivative ***"<<endl;
	//DPhiFinal.Print();
	//cout<<endl;

	TMatrixD DThetaFinalMultCov(1,nPar);
	DThetaFinalMultCov.Mult(DThetaFinal,finalTrackCovariance);
	TMatrixD DPhiFinalMultCov(1,nPar);
	DPhiFinalMultCov.Mult(DPhiFinal,finalTrackCovariance);
		
	//cout<<"*** Final Covariance ***"<<endl;
	//finalTrackCovariance.Print();
	//cout<<endl;

	//cout<<"*** DThetaFinalMultCov ***"<<endl;
	//DThetaFinalMultCov.Print();
	//cout<<endl;

	TMatrixD ThetaFinalUncertainty(1,1);
	ThetaFinalUncertainty.Mult(DThetaFinalMultCov,DThetaFinalT);

	TMatrixD PhiFinalUncertainty(1,1);
	PhiFinalUncertainty.Mult(DPhiFinalMultCov,DPhiFinalT);

	fThetaErr= sqrt(ThetaFinalUncertainty(0,0))* 180./TMath::Pi();
	fPhiErr= sqrt(PhiFinalUncertainty(0,0))* 180./TMath::Pi();

	cout<<"*** KALMAN FILTER TRACK RECO ***"<<endl;
	cout<<"Theta [deg]="<<fTheta<<" +- "<< fThetaErr << "  Phi [deg]="<<fPhi<<" +- "<<fPhiErr<<endl;
	cout<<"X0 [cm]="<< fVertexX <<" +- "<< fVertexXErr << "  Y0[cm]="<<fVertexY<<" +- "<< fVertexYErr<<endl;
	cout<<"Tx="<<fTx<<" +- "<< fTxErr <<"  Ty="<<fTy<<" +- "<< fTyErr <<endl;
	cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 	cout<<"*******************************"<<endl;

	return 0;

}//close KalmanFilterTrackFitter()



void KFTrackFinder::KFPredictionStep(int k){
	
	//cout<<"*******************************"<< endl;
	//cout<<"** KFPredictionStep --> "<< k << "**"<< endl;
	//cout<<"*******************************"<< endl;

	//## Calculating predicted state k from k-1 state: x_pred(k)= F(k)x(k-1)
	TMatrixD F= fKFPropagator[k-1];
	TMatrixD X= fKFState[k-1];
	TMatrixD C= fKFCovariance[k-1];
	int nRow= F.GetNrows();
	int nCol= X.GetNcols();

	TMatrixD XPred(nRow,nCol);
	XPred.Mult(F,X);
	//cout<<"*** Prediction State XPred ***"<<endl;
	//XPred.Print();
	//cout<<endl;

	//## Copying to global predicted state vector
	fKFPredictedState[k-1]= XPred;

	//## Calculate covariance matrix due to Multiple Scattering
	MSCovariance(k);
	TMatrixD Q= fKFCovarianceMS[k-1];
	//cout<<"*** Noise Q ***"<<endl;
	//Q.Print();
	//cout<<endl;

	//## Calculating predicted covariance state k from k-1 state: C_pred(k)= F(k)C(k-1)Ft(k) + Q(k)	
	TMatrixD Ft(nPar,nPar);
	Ft.Transpose(F);

	TMatrixD FMultC(nPar,nPar);
	FMultC.Mult(F,C);

	TMatrixD FCMultFt(nPar,nPar);
	FCMultFt.Mult(FMultC,Ft);
	//cout<<"*** Prediction Det Covariance FxCxFt ***"<<endl;
	//FCMultFt.Print();
	//cout<<endl;

	TMatrixD CPred(nPar,nPar);
	if(fTrackWithMSCorrection)
		CPred= FCMultFt + Q;
	else	
		CPred= FCMultFt;

	//cout<<"*** Prediction Covariance CPred ***"<<endl;
	//CPred.Print();
	//cout<<endl;

	//## Copying to global predicted covariance matrix
	fKFPredictedCovariance[k-1]= CPred;

	//## Calculating meas residuals
	TMatrixD H= fKFProjector[k-1];
	//cout<<"*** Projector H ***"<<endl;
	//H.Print();
	//cout<<endl;

	TMatrixD M= fKFMeas[k-1];
	TMatrixD V= fKFMeasVariance[k-1];	
	TMatrixD Ht(H.GetNcols(),H.GetNrows());
	Ht.Transpose(H);
	//cout<<"*** Transpose Projector Ht ***"<<endl;
	//Ht.Print();
	//cout<<endl;


	int nMeas= M.GetNrows();
	TMatrixD HMultXPred(nMeas,1);
	HMultXPred.Mult(H,XPred);
	//cout<<"*** H x XPred ***"<<endl;
	//HMultXPred.Print();
	//cout<<endl;
	

	TMatrixD rPred(nMeas,1);
	rPred= M-HMultXPred;
	//cout<<"*** Predicted Measurement Residual rPred ***"<<endl;
	//rPred.Print();
	//cout<<endl;
	
	fKFPredictedMeasResidual[k-1]= rPred;

	//## Calculating covariance residuals
  TMatrixD HMultCPred(nMeas,nPar);
	HMultCPred.Mult(H,CPred);

	TMatrixD HCPredMultHt(nMeas,nMeas);
	HCPredMultHt.Mult(HMultCPred,Ht);

	TMatrixD RPred(nMeas,nMeas);
	RPred= V + HCPredMultHt;

	//cout<<"*** Predicted Covariance Measurement Residual RPred ***"<<endl;
	//RPred.Print();
	//cout<<endl;

	fKFPredictedCovarianceResidual[k-1]= RPred;

	//## Calculating predicted Chi2: rtPred x RPred^-1 x rPred 
	TMatrixD rtPred(1,nMeas);
	rtPred.Transpose(rPred);

	double determ = 0.;
  TMatrixD RPredInv = RPred;
  RPredInv = RPredInv.Invert(&determ);
  if (determ<=0) {
		cerr<<"KFTrackFinder::KFPredictionStep(): WARNING: predicted residual inversion failed!"<<endl;
	}

	TMatrixD rtPredMultRPredInv(1,nMeas);
	rtPredMultRPredInv.Mult(rtPred,RPredInv);

	TMatrixD PredChi2Matrix(1,1);
	PredChi2Matrix.Mult(rtPredMultRPredInv,rPred);

	//cout<<"*** Predicted Chi2 ***"<<endl;
	//PredChi2Matrix.Print();
	//cout<<endl;

	double PredChi2= PredChi2Matrix(0,0);
	fKFHitPredictedChi2[k-1]= PredChi2;

}//close KFTrackFinder::KFPredictionStep

void KFTrackFinder::KFFilterStep(int k){
	
	//cout<<"*******************************"<< endl;
	//cout<<"** KFFilterStep --> "<< k << "**"<< endl;
	//cout<<"*******************************"<< endl;

	//## Calculate Gain Matrix: K(k)= CPred(k)Ht x (V(k)+H(k)CPred(k)Ht(k))^-1= CPred(k)Ht R(k)^-1
	TMatrixD XPred= fKFPredictedState[k-1];
	TMatrixD CPred= fKFPredictedCovariance[k-1];
	TMatrixD H= fKFProjector[k-1];
	TMatrixD Ht(H.GetNcols(),H.GetNrows());
	Ht.Transpose(H);
	TMatrixD M= fKFMeas[k-1];
	TMatrixD V= fKFMeasVariance[k-1];	

	int nMeas= M.GetNrows();

	TMatrixD RPred= fKFPredictedCovarianceResidual[k-1];
	TMatrixD rPred= fKFPredictedMeasResidual[k-1];
		

	double determ = 0.;
  TMatrixD RPredInv = RPred;
  RPredInv = RPredInv.Invert(&determ);
  if (determ<=0) {
		cerr<<"KFTrackFinder::KFFilterStep(): WARNING: predicted residual inversion failed!"<<endl;
	}
	//cout<<"*** Inverted Predicted Measurement Residual RPred^-1 ***"<<endl;
	//RPredInv.Print();
	//cout<<endl;

	TMatrixD CPredMultHt(nPar,nMeas);
	CPredMultHt.Mult(CPred,Ht);
	
	TMatrixD K(nPar,nMeas);
	K.Mult(CPredMultHt,RPredInv);
	//cout<<"## Gain K ##"<<endl;
	//K.Print();
	//cout<<endl;
	
	//## Calculating state vector x(k)= xPred(k)+ K(k) x r(k)
	TMatrixD KMultrPred(nPar,1);
	KMultrPred.Mult(K,rPred);
	
	TMatrixD X(nPar,1);
	X= XPred + KMultrPred;
	//cout<<"## Filtered State X ##"<<endl;
	//X.Print();
	//cout<<endl;

	fKFState[k]= X;

	//## Calculating covariance matrix C(k)= (1-K(k)H(k))x CPred(k)
	TMatrixD KMultH(nPar,nPar);
	KMultH.Mult(K,H);

	TMatrixD Identity(nPar,nPar);
	for(int i=0;i<nPar;i++){
		for(int j=0;j<nPar;j++){	
			if(i==j) Identity(i,j)= 1.;
			else Identity(i,j)= 0.;
		}
	}

	TMatrixD IMinusKMultH(nPar,nPar);
	IMinusKMultH= Identity-KMultH;

	//TMatrixD KHMultCPred(nPar,nPar);
	//KHMultCPred.Mult(KMultH,CPred);
	
	TMatrixD C(nPar,nPar);
	//C= CPred-KHMultCPred;
	C.Mult(IMinusKMultH,CPred);
	//cout<<"*** Filtered Covariance C ***"<<endl;
	//C.Print();
	//cout<<endl;
	
	fKFCovariance[k]= C;

	//## Calculating filtered residuals
	TMatrixD HMultK(nMeas,nMeas);
	HMultK.Mult(H,K);
	TMatrixD HKMultrPred(nMeas,1);
	HKMultrPred.Mult(HMultK,rPred);

	TMatrixD r= rPred-HKMultrPred;
	//cout<<"*** Filtered Measurement Residual r ***"<<endl;
	//r.Print();
	//cout<<endl;

	//## Calculating filtered covariance residuals
	TMatrixD HKMultV(nMeas,nMeas);
	HKMultV.Mult(HMultK,V);

	TMatrixD R= V-HKMultV;
	//cout<<"*** Filtered Measurement Covariance Residual R ***"<<endl;
	//R.Print();
	//cout<<endl;

	//## Calculating Chi2
	TMatrixD rt(1,nMeas);
	rt.Transpose(r);
	
	determ = 0.;
  TMatrixD RInv = R;
  RInv = RInv.Invert(&determ);
  if (determ<=0) {
		cerr<<"KFTrackFinder::KFFilterStep(): WARNING: predicted residual inversion failed!"<<endl;
	}
	//cout<<"*** Inverted Measurement Residual R^-1 ***"<<endl;
	//RInv.Print();
	//cout<<endl;

	TMatrixD rtMultRInv(1,nMeas);
	rtMultRInv.Mult(rt,RInv);
	
	TMatrixD Chi2Matrix(1,1);
	Chi2Matrix.Mult(rtMultRInv,r);
	//cout<<"*** Chi2 ***"<<endl;
	//Chi2Matrix.Print();
	//cout<<endl;

	double Chi2= Chi2Matrix(0,0);
	fFitFCNMin+= Chi2;
	//TrackHitChi2[k-1]= Chi2;
	fKFHitChi2[k-1]= Chi2;

}//close KFTrackFinder::KFFilterStep


void KFTrackFinder::KFSmoothingStep(int k){

	//cout<<"*******************************"<< endl;
	//cout<<"** KFSmoothingStep --> "<< k << "**"<< endl;
	//cout<<"*******************************"<< endl;
	if(k>nHits){
		cerr<<"KFTrackFinder::KSSmoothingStep(): ERROR: the index given to smoother must be < "<< nHits <<" ...exit"<<endl;
		exit(1);
	}
	else if(k==nHits){
		fKFSmoothedState[k]= fKFState[k]; 
		fKFSmoothedCovariance[k]= fKFCovariance[k];
	}

	else{
		//## Calculate smoother matrix A(k)= C(k)Ft(k+1) CPred(k+1)^-1 
		TMatrixD C= fKFCovariance[k];
	
		TMatrixD F= fKFPropagator[k];
		int nPar= F.GetNrows();
		TMatrix Ft(nPar,nPar);
		Ft.Transpose(F);

		TMatrix CPred= fKFPredictedCovariance[k];
		double determ = 0.;
  	TMatrixD CPredInv = CPred;
  	CPredInv = CPredInv.Invert(&determ);
  	if (determ<=0) {
			cerr<<"KFTrackFinder::KFSmoothingStep(): WARNING: predicted covariance inversion failed!"<<endl;
		}
	
		TMatrixD CMultFt(nPar,nPar);
		CMultFt.Mult(C,Ft);

		TMatrixD A(nPar,nPar);
		A.Mult(CMultFt,CPredInv);
		//cout<<"*** Smoother Matrix A ***"<<endl;
		//A.Print();
		//cout<<endl;

		TMatrixD At(nPar,nPar);
		At.Transpose(A);

		//## Calculate smoothed state vector: Xsmooth(k,n)= X(k) + A(k) (X(k+1,n)-X(k+1,k)) 
		TMatrixD X= fKFState[k];
		//cout<<"*** State ***"<<endl;
		//X.Print();
		//cout<<endl;
		
		TMatrixD XPred= fKFPredictedState[k];
		TMatrixD XSmoothedNext= fKFSmoothedState[k+1];

		TMatrixD StateDiff(nPar,1);
		StateDiff.Minus(XSmoothedNext,XPred);
		//cout<<"*** State Diff ***"<<endl;
		//StateDiff.Print();
		//cout<<endl;

		TMatrixD AMultStateDiff(nPar,1);
		AMultStateDiff.Mult(A,StateDiff);
	
		TMatrixD XSmoothed(nPar,1);
		XSmoothed= X+AMultStateDiff;
		//cout<<"*** State Smoothed ***"<<endl;
		//XSmoothed.Print();
		//cout<<endl;

		fKFSmoothedState[k]= XSmoothed; 

		//## Calculate smoothed covariance: Csmooth(k,n)= C(k) + A(k) (C(k+1,n)-C(k+1,k)) At(k)
		TMatrixD CSmoothedNext= fKFSmoothedCovariance[k+1];
		//cout<<"*** Covariance Smoothed Next ***"<<endl;
		//CSmoothedNext.Print();
		//cout<<endl;

		TMatrixD CovDiff(nPar,nPar);
		CovDiff.Minus(CSmoothedNext,CPred);

		TMatrixD AMultCovDiff(nPar,nPar);
		AMultCovDiff.Mult(A,CovDiff);

		TMatrixD ACovDiffMultAt(nPar,nPar);
		ACovDiffMultAt.Mult(AMultCovDiff,At);

		TMatrixD CSmoothed(nPar,nPar);
		CSmoothed= C + ACovDiffMultAt;
		//cout<<"*** Covariance Smoothed ***"<<endl;
		//CSmoothed.Print();
		//cout<<endl;

		fKFSmoothedCovariance[k]= CSmoothed;
	}
	

}//close KFTrackFinder::KSSmoothingStep()


void KFTrackFinder::MSCovariance(int k){
	
	//cout<<"### KFTrackFinder::MSCovariance()"<<endl;

	//## Check
	if(k==0){
		cerr<<"KFTrackFinder::MSCovariance(): index must be >0"<<endl;
		exit(1);	
	}

	TMatrixD X= fKFState[k-1]; 
	//cout<<"*** Printing X ***"<<endl;
	//X.Print();

	int nPar= X.GetNrows();
	double tx= X(2,0); 
	double ty= X(3,0);
	//cout<<"tx="<<tx<<"  ty="<<ty<<endl;	
	
	//## Path lenghts
	double Lr= 21.997;// radiation length of Malargue soil in g/cm^2 
	double Density= 1.8;// average density of Malargue soil in g/cm^3
	double l= Ls[k-1];
	double ZPath= fabs(lenz[k-1]);
	double ZGrammagePath= ZPath*Density;
	double GrammagePathLength= ZGrammagePath*sqrt(1+tx*tx+ty*ty);
	double averageEnergyDeposit= 1.8;//MeV/(g/cm^2)
	
	//## Energy and Momentum
	double MuonMass= 105.65836668;//in MeV/c^2
	double TrueKinEnergy= fGenEnergy;
	double TrueMomentum= sqrt( TrueKinEnergy*(TrueKinEnergy+2.*MuonMass) );//in MeV/c
	
	double KinEnergy= 5000.;//in MeV
	double Momentum= sqrt( KinEnergy*(KinEnergy+2.*MuonMass) );//in MeV/c
	double AverageMomentum= sqrt( fAverageEnergyInMS*(fAverageEnergyInMS+2.*MuonMass) );//in MeV/c
	
	double correctedKinEnergy= fKinEnergy-averageEnergyDeposit*GrammagePathLength;
	fKinEnergy= correctedKinEnergy;
	fMomentum= sqrt( fKinEnergy*(fKinEnergy+2.*MuonMass) );//in MeV/c
	//cout<<"ZPath [cm]="<<ZPath<<"  GrammagePathLength[g/cm^2]="<<GrammagePathLength<<endl;
	//cout<<"TrueMomentum [MeV]="<<TrueMomentum<<endl;
	//cout<<"KinEnergy [MeV]="<<fKinEnergy<<"  eLoss [MeV]="<<averageEnergyDeposit*GrammagePathLength<<endl;
	//cout<<"NewKinEnergy [MeV]="<<fKinEnergy<<"  NewMomentum[MeV]="<<fMomentum<<endl;
	
	//double H= pow(13.6/Momentum,2);
	//double H= pow(13.6/TrueMomentum,2);
	//double H= pow(13.6/fMomentum,2);

	double H;
	if(fUseAverageEnergyInMS){
		H= pow(13.6/AverageMomentum,2);
	}
	else if(fUseTrueEnergyInMS){
		H= pow(13.6/TrueMomentum,2);
	}
	else if(fUseEnergyLossCorrectionInMS){
		H= pow(13.6/fMomentum,2);
	}
	else{
		cerr<<"KFTrackFinder::MSCovariance(): Specified fit with MS on but no choice given for MS treatment...exit"<<endl;
		exit(1);
	}

	//## Calculate matrix elements
	double RadThickness= l*Density/Lr;
	double EffectiveRadThickness= RadThickness*sqrt(1+tx*tx+ty*ty);
	double EffectivePathLength= l*sqrt(1+tx*tx+ty*ty);
	double CMS= H*EffectiveRadThickness*pow( (1+0.038*log(EffectiveRadThickness)), 2);
 	//cout<<"l="<<l<<"  RadThickness="<<RadThickness<<"  EffectiveRadThickness="<<EffectiveRadThickness<<"  CMS="<<CMS<<endl;

	double cov_txtx= (1+tx*tx)*(1+tx*tx+ty*ty)*CMS;
	double cov_tyty= (1+ty*ty)*(1+tx*tx+ty*ty)*CMS;
	double cov_txty= tx*ty*(1+tx*tx+ty*ty)*CMS;
	//cout<<"cov_txtx="<<cov_txtx<<"  cov_tyty="<<cov_tyty<<"  cov_txty="<<cov_txty<<endl;

	double D= -1;

	TMatrixD Q(nPar,nPar);
	double leff= EffectivePathLength;
	//double leff= EffectiveRadThickness;
	Q(0,0)= cov_txtx*leff*leff/3.;
	Q(0,1)= cov_txty*leff*leff/3.;
	Q(1,0)= Q(0,1);
	Q(0,2)= cov_txtx*D*leff/2.;
	Q(2,0)= Q(0,2);
	Q(0,3)= cov_txty*D*leff/2.;
	Q(3,0)= Q(0,3);
		
	Q(1,1)= cov_tyty*leff*leff/3.;
	Q(1,2)= cov_txty*D*leff/2.;
	Q(2,1)= Q(1,2);
	Q(1,3)= cov_tyty*D*leff/2.;
	Q(3,1)= Q(1,3);

	Q(2,2)= cov_txtx;
	Q(2,3)= cov_txty;
	Q(3,2)= Q(2,3);
	
	Q(3,3)= cov_tyty;

	fKFCovarianceMS[k-1]= Q;

}//close KFTrackFinder::MSCovariance()



void KFTrackFinder::DebugPrint(){

	cout<<"*** Printing Measurement Matrix M ***"<<endl;
	for(unsigned int i=0;i<fKFMeas.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFMeas[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing Measurement Variance Matrix V ***"<<endl;
	for(unsigned int i=0;i<fKFMeasVariance.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFMeasVariance[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing Propagator Matrix F ***"<<endl;
	for(unsigned int i=0;i<fKFPropagator.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFPropagator[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing Projector Matrix H ***"<<endl;
	for(unsigned int i=0;i<fKFProjector.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFProjector[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing State X ***"<<endl;
	for(unsigned int i=0;i<fKFState.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFState[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing Smoothed State XSmoothed ***"<<endl;
	for(unsigned int i=0;i<fKFSmoothedState.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFSmoothedState[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing Predicted State XPred ***"<<endl;
	for(unsigned int i=0;i<fKFPredictedState.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFPredictedState[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing Covariance C ***"<<endl;
	for(unsigned int i=0;i<fKFCovariance.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFCovariance[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing Smoothed Covariance XSmoothed ***"<<endl;
	for(unsigned int i=0;i<fKFSmoothedCovariance.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFSmoothedCovariance[i].Print();	
	}
	cout<<endl;

	cout<<"*** Printing Predicted Covariance CPred ***"<<endl;
	for(unsigned int i=0;i<fKFPredictedCovariance.size();i++){
		cout<<"--> state= "<< i << endl;
		fKFPredictedCovariance[i].Print();	
	}
	cout<<endl;


}//close KFTrackFinder::DebugPrint()


int KFTrackFinder::LSTrackFitter(Track* aTrack){

	//cout<<"*********************************"<<endl;
	//cout<<"***     LS TRACK FITTER       ***"<<endl;
	//cout<<"*********************************"<<endl;
	double x,y,z,time;
	double xErr,yErr,zErr,timeErr;
	int nHits= aTrack->GetNpoints();
	
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
		//cout<<"Z["<<i<<"]="<< hits[i][2]<<endl;
		//cout<<"Sigma["<<i<<"]="<< sigma[i][0]<<"  "<<sigma[i][1]<<endl;

		int ix = i;
		int iy = i+nHits;
		m[ix][0]= x;
		m[iy][0]= y;
	}	

	P0[2] = hits[0][2];

	//## Lenghts (Zi-Z0)
	double lenz[nHits];
	for(int i=0;i<nHits;i++) lenz[i] = hits[i][2] - P0[2];

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
	//cout<<"## Derivative Matrix F ##"<<endl;
	//F.Print();
	//cout<<endl;

	//## Calculate transpose of derivative matrix Ft
	TMatrixD Ft= TMatrixD(idim,2*nHits);
	Ft.Transpose(F);
	//cout<<"## Transpose Derivative Matrix Ft ##"<<endl;
	//Ft.Print();
	//cout<<endl;

	//## Calculate Covariance Matrix V
	double v[2*nHits][2*nHits];
	
	for(int i=0;i<2*nHits;i++){
		for(int j=0;j<2*nHits;j++) v[i][j]= 0.;
	}
	
	//## diagonal elements
	for(int i=0;i<nHits;i++) {
		int ix = i;
		int iy = i+nHits;
		v[ix][ix]= pow(sigma[i][0],2);
		v[iy][iy]= pow(sigma[i][1],2);
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
		cerr<<"KFTrackFinder::LSTrackFitter(): Covariance Matrix inversion failed!"<<endl;
		return -1;
	}
	//cout<<"## Inverted Covariance Matrix V ##"<<endl;
	//VInv.Print();
	//cout<<endl;

	//## Calculate Matrix (Ft VInv F)	= C covariance
	TMatrix FtVInvProduct= TMatrix(idim,2*nHits); 
	FtVInvProduct.Mult(Ft,VInv);
	//cout<<"## Ft x VInv ##"<<endl;
	//FtVInvProduct.Print();
	//cout<<endl;

	TMatrix C= TMatrixD(idim,idim);
	C.Mult(FtVInvProduct,F);
	//cout<<"## Covariance C ##"<<endl;
	//C.Print();
	//cout<<endl;

	//## Invert the covariance C
	determ = 0.0;
  TMatrixD CInv = C;
  CInv = CInv.Invert(&determ);
  if (determ<=0) {
		cerr<<"KFTrackFinder::LSTrackFitter(): Final Covariance Matrix inversion failed!"<<endl;
		return -1;
	}
	//cout<<"## Inverted Covariance Matrix C ##"<<endl;
	//CInv.Print();
	//cout<<endl;

	//## Calculate CInv Ft VInv M
	TMatrixD CInvFtProduct= TMatrixD(idim,2*nHits);
	CInvFtProduct.Mult(CInv,Ft);
	//cout<<"## CInv x Ft ##"<<endl;
	//CInvFtProduct.Print();
	//cout<<endl;

	TMatrixD CInvFtVInvProduct= TMatrixD(idim,2*nHits);
	CInvFtVInvProduct.Mult(CInvFtProduct,VInv);
	//cout<<"## CInv x Ft x VInv ##"<<endl;
	//CInvFtVInvProduct.Print();
	//cout<<endl;


	TMatrixD FitParam= TMatrixD(idim,1);
	FitParam.Mult(CInvFtVInvProduct,M);

	//cout<<"## Fit Params ##"<<endl;
	//FitParam.Print();
	//cout<<endl;

	//## Find track solution
	double par[idim];
 	for (int k=0;k<idim;k++) {
		par[k] = FitParam(k,0);
 	}
 	
	
	//## Calculate (M-F lambda)
	TMatrixD FParProduct= TMatrixD(2*nHits,1);
	FParProduct.Mult(F,FitParam);
	//cout<<"## F x Lambda ##"<<endl;
	//FParProduct.Print();
	//cout<<endl;


	TMatrixD R= TMatrixD(2*nHits,1);
	R= M-FParProduct;
	//cout<<"## M- F x Lambda ##"<<endl;
  //R.Print();
	//cout<<endl;


	TMatrixD Rt= TMatrixD(1,2*nHits);
	Rt.Transpose(R);
	//cout<<"## (M- F x Lambda)T ##"<<endl;
  //Rt.Print();
	//cout<<endl;


	//## Calculate Chi2
	TMatrixD RtVInvProduct= TMatrixD(1,2*nHits);
	RtVInvProduct.Mult(Rt,VInv);
	//cout<<"## (M- F x Lambda)T x VInv ##"<<endl;
  //RtVInvProduct.Print();
	//cout<<endl;
	

	TMatrixD Chi2Matrix= TMatrixD(1,1);
	Chi2Matrix.Mult(RtVInvProduct,R); 
	//cout<<"## (M- F x Lambda)T x VInv x (M- F x Lambda) ##"<<endl;
  //Chi2Matrix.Print();
	//cout<<endl;

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
	
	
	fTheta= direction.Theta()*180./TMath::Pi();
	fPhi= direction.Phi()*180./TMath::Pi();

	fVertexX= x0.X();
	fVertexY= x0.Y();	
	fTx= par[2];
	fTy= par[3];


	fFitFCNMin= Chi2;

	fKFInitState.ResizeTo(idim,1);
	fKFInitState(0,0)= fVertexX;
	fKFInitState(1,0)= fVertexY;
	fKFInitState(2,0)= fTx;
	fKFInitState(3,0)= fTy;

	fKFInitCovariance.ResizeTo(idim,idim);
	fKFInitCovariance= CInv;

	cout<<"*** LS TRACK RECO ***"<<endl;
	cout<<"Theta [deg]="<<fTheta<<"  Phi [deg]="<<fPhi<<endl;	
	cout<<"X0 [cm]="<<fVertexX<<"  Y0 [cm]="<<fVertexY<<endl;
	cout<<"Tx="<<fTx<<"  Ty="<<fTy<<endl;
	cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 	cout<<"*******************************"<<endl;

	return 0;

}//close LSTrackFitter()


bool KFTrackFinder::Checker(){

	if(fCandidateTrackPointCollection.empty()){
		cout<<"KFTrackFinder::Checker(): No track points given...exit"<<endl;
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
		cout<<"KFTrackFinder::Checker(): Cannot identify tracks from only "<<fTrackPointCollectionInPlane.size()<<" XY planes...exit"<<endl;
		return false;
	}
	else{
		cout<<"*** TRACK POINTS INFO ***"<<endl;
		for(unsigned int i=0;i<fTrackPointCollectionInPlane.size();i++) {
			cout<<"PLANE "<<i+1<<": "<< fTrackPointCollectionInPlane[i].size() << " points"<<endl;
		}		
		cout<<"==> "<<GetCombinationNo()<<" potential candidate tracks to be fitted"<<endl;
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

}//close KFTrackFinder::Checker()


int KFTrackFinder::GetCombinationNo(){
        
	int n = 1;
  for(unsigned int i=0;i<fTrackPointCollectionInPlane.size();i++) n*= fTrackPointCollectionInPlane[i].size();
        
	return n;

}//close KFTrackFinder::GetCombinationNo()

/*
void KFTrackFinder::CalculateCombinationMatrix(){
        
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
	//int nGoodCombinations= 0;
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
	}	//close if verbosity

}//close KFTrackFinder::CalculateCombinationMatrix()()
*/

void KFTrackFinder::CalculateCombinationMatrix(){

	int combNo = GetCombinationNo();
	int planeNo = (int)(fTrackPointCollectionInPlane.size());
	int dummyComb[planeNo];  

	//## clear existing combination matrix
	fCombinationMatrix.clear();
	fCombinationMatrix.resize(0);


	//## compute combination indexes
	std::vector< std::vector<int> > combIndex;
	
	for(int i=0;i<planeNo;i++){
		combIndex.push_back( std::vector<int>() );
		for(int j=0;j<fTrackPointCollectionInPlane[i].size();j++){
			combIndex[i].push_back(j);	
		}
	}

	//## Combine hits from all planes	
	//## Select good hit combination to reduce allocated memory and speed-up the tracking
	int nGoodCombination= 0;

  for (unsigned int i=0;i<combIndex[0].size();++i) {
  	for (unsigned int j=0;j<combIndex[1].size();++j) {
			//cout<<"Comb: ";
			for (unsigned int k=0;k<combIndex[2].size();++k) {
				dummyComb[0]= combIndex[0][i];
				dummyComb[1]= combIndex[1][j];
				dummyComb[2]= combIndex[2][k];

				//## Accept hit combination?
				//cout<<dummyComb[0]<<"  "<<dummyComb[1]<<"  "<<dummyComb[2]<<endl;

				if(! AcceptHitCombination(dummyComb) )
					continue;

				fCombinationMatrix.push_back( std::vector<int>() );
				fCombinationMatrix[nGoodCombination].push_back(combIndex[0][i]);	
				fCombinationMatrix[nGoodCombination].push_back(combIndex[1][j]);
				fCombinationMatrix[nGoodCombination].push_back(combIndex[2][k]);		
    		nGoodCombination++;
			}//end loop k
    }//end loop j
  }//end loop i

	cout<<"Tot/Sel combinations: "<<combNo<<"/"<<nGoodCombination<<endl; 


}//close KFTrackFinder::CalculateCombinationMatrix()


bool KFTrackFinder::AcceptHitCombination(int* aCandidateCombination){

	//## Accept a candidate track
	//## Compare relative angles in space of tracks formed with each point
	double x,y,z,t;
	int nPoints= (int)(fTrackPointCollectionInPlane.size());
	TrackPoint thisTrackPoints[nPoints];
	TVector3 thisTrackPointVector[nPoints];
	TVector3 thisTrackVector[nPoints];
	double RelAngleAmongTracks[nPoints];	

	
	for(int p=0;p<nPoints;p++){
		thisTrackPoints[p]= fTrackPointCollectionInPlane[p][*(aCandidateCombination+p)];
		thisTrackPointVector[p]= thisTrackPoints[p].fPosition;
	}

	for(int p=0;p<nPoints-1;p++){
		thisTrackVector[p]= thisTrackPointVector[p+1]-thisTrackPointVector[p];
	}
	thisTrackVector[nPoints-1]= thisTrackPointVector[nPoints-1]-thisTrackPointVector[0];

	//## Calculate relative angle in space
	for(int p=0;p<nPoints-1;p++){
		RelAngleAmongTracks[p]= thisTrackVector[p].Angle(thisTrackVector[p+1]) * 180./TMath::Pi();
		//cout<<"RelAngle "<<p<<"-"<<p+1<<" ="<< RelAngleAmongTracks[p]<<endl;
	}
	RelAngleAmongTracks[nPoints-1]= thisTrackVector[nPoints-1].Angle(thisTrackVector[0]) * 180./TMath::Pi();
	//cout<<"RelAngle 0-"<<nPoints-1<<" ="<< RelAngleAmongTracks[nPoints-1]<<endl;
	
	//## Accept candidate track?
	for(int p=0;p<nPoints-1;p++){
		if(RelAngleAmongTracks[p]>fMaxAllowedAngle)
			return false;
	}
	

	return true;

}// close KFTrackFinder::AcceptHitCombination()



void KFTrackFinder::RemoveFakeHits(){

	//## Loop over track hits and removing those not due to the muon

	std::vector<TrackPoint>::iterator it= fTrackPointCollection.begin();
	cout<<"KFTrackFinder::RemoveFakeHits(): "<< fTrackPointCollection.size() <<"  total hits are present"<<endl;
	int nSpurious= 0;
	std::vector<int> pointToBeErased;
	pointToBeErased.clear();	
	pointToBeErased.resize(0);

	for(unsigned int i=0;i<fTrackPointCollection.size();i++){	
		int isMuonHit= fTrackPointCollection[i].fIsMuon;
		//cout<<"point "<<i<<" ==>" <<isMuonHit<<endl;
		//(fTrackPointCollection[i].fPosition).Print();
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


	//cout<<endl;
	cout<<"KFTrackFinder::RemoveFakeHits(): "<< nSpurious <<"  fake hits removed, "<<fTrackPointCollection.size()<<"  left"<<endl;
	for(unsigned int i=0;i<fTrackPointCollection.size();i++){	
		int isMuonHit= fTrackPointCollection[i].fIsMuon;
		//cout<<"point "<<i<<" ==>" <<isMuonHit<<endl;
		//(fTrackPointCollection[i].fPosition).Print();
	}

}//close KFTrackFinder::RemoveFakeHits()


bool KFTrackFinder::FindCluster(){
	
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
						//cout<<"Joined to current point"<<endl;
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

		//### MERGE CLUSTER HITS
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

		//cout<<endl;
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
}//close KFTrackFinder::ClusterFinder()


void KFTrackFinder::SetTrackPointUncertainty(){

	//## Assign point uncertainty to track points
	TVector3 pointErr= TVector3(fHitSigmaX,fHitSigmaY,fHitSigmaZ);
	for(unsigned int i=0;i< fTrackPointCollection.size();i++){
		fTrackPointCollection[i].fPositionErr= pointErr;
		//fTrackPointCollection[i].fPositionErr.Print();
	}

}//close KFTrackFinder::SetTrackPointUncertainty()


int KFTrackFinder::KFTrackFitterWithMomentum(Track* aTrack){

	double x,y,z,time;
	double xErr,yErr,zErr,timeErr;
	int Npoints= aTrack->GetNpoints();
	if (Npoints<3) return false;

	//## fit the graph now 
  TVirtualFitter* fit = TVirtualFitter::Fitter(0,1);
  fit->SetObjectFit(aTrack);
  fit->SetFCN(KFTrackFinder::TrackChi2Fcn);
   
  double arglist[10];
  arglist[0] = 3;
  fit->ExecuteCommand("SET PRINT",arglist,1);
  
  double pStart[1] = {0};
  fit->SetParameter(0,"momentum",pStart[0],0.01,-1,3);
   
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
  double parFit[1];
  for(int i=0;i<1;++i) parFit[i] = fit->GetParameter(i);

	cout<<"*** FAST TRACK RECO WITH MS ***"<<endl;
	cout<<"Theta [deg]="<<fTheta<<"  Phi [deg]="<<fPhi<<endl;
	cout<<"Momentum [MeV]="<<fMomentum<<endl;
	cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 	cout<<"*******************************"<<endl;

	return fFitStatus;
  
}//close KFTrackFinder::TrackFitterWithMS()


void KFTrackFinder::TrackChi2Fcn(int& nFitPar, double* const grad,double& value, double* fitpar,const int iFlag){
    
	Track* aTrack = dynamic_cast<Track*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
  assert(aTrack != 0);

	cout<<"fitpar[0]="<<fitpar[0]<<endl;
	fMomentum= pow(10,fitpar[0])*1.e-3;//convert in MeV
	
 
	//## Consistency checks
	nHits= aTrack->GetNpoints();
	
	int nIter= nHits+1;
	nPar= 4;
	int nMeas= 2;

	//## Init vectors
	fKFState.clear();
	fKFState.resize(0);
	//fKFState.resize(nIter);

	fKFCovariance.clear();
	fKFCovariance.resize(0);
	//fKFCovariance.resize(nIter);

	fKFSmoothedState.clear();
	fKFSmoothedState.resize(0);
	//fKFSmoothedState.resize(nIter);

	fKFSmoothedCovariance.clear();
	fKFSmoothedCovariance.resize(0);
	//fKFSmoothedCovariance.resize(nIter);

	fKFPredictedState.clear();
	fKFPredictedState.resize(0);
	//fKFPredictedState.resize(nHits);
		
	fKFPredictedCovariance.clear();
	fKFPredictedCovariance.resize(0);
	//fKFPredictedCovariance.resize(nHits);

	fKFProjector.clear();
	fKFProjector.resize(0);
	//fKFProjector.resize(nHits);

	fKFMeas.clear();
	//fKFMeas.resize(nHits);
	fKFMeas.resize(0);

	fKFMeasVariance.clear();
	//fKFMeasVariance.resize(nHits);
	fKFMeasVariance.resize(0);

	fKFPropagator.clear();
	//fKFPropagator.resize(nHits);
	fKFPropagator.resize(0);
	
	fKFCovarianceMS.clear();
	fKFCovarianceMS.resize(0);
	//fKFCovarianceMS.resize(nHits);

	fKFCovarianceResidual.clear();
	fKFCovarianceResidual.resize(0);
	//fKFCovarianceResidual.resize(nHits);

	fKFMeasResidual.clear();
	fKFMeasResidual.resize(0);
	//fKFMeasResidual.resize(nHits);

	fKFSmoothedCovarianceResidual.clear();
	fKFSmoothedCovarianceResidual.resize(0);
	//fKFSmoothedCovarianceResidual.resize(nHits);

	fKFSmoothedMeasResidual.clear();
	fKFSmoothedMeasResidual.resize(0);
	//fKFSmoothedMeasResidual.resize(nHits);

	fKFPredictedCovarianceResidual.clear();
	fKFPredictedCovarianceResidual.resize(0);
	//fKFPredictedCovarianceResidual.resize(nHits);

	fKFPredictedMeasResidual.clear();
	fKFPredictedMeasResidual.resize(0);
	//fKFPredictedMeasResidual.resize(nHits);

	fKFHitChi2.clear();
	fKFHitChi2.resize(0);

	fKFHitPredictedChi2.clear();
	fKFHitPredictedChi2.resize(0);
	
	for(int i=0;i<nHits;i++) {
		TrackHitChi2[i]= 0.;
		TrackHitPredictedChi2[i]= 0.;
		fKFHitChi2.push_back(0.);
		fKFHitPredictedChi2.push_back(0.);
	}
	
	//## Choose initial state vector and covariance
	//int LSStatus= LSTrackFitter(aTrack);
	//if(LSStatus!=0){
	//	cerr<<"KFTrackFinder::KalmanFilterTrackFitter(): WARNING: LSFitter (used to init the KFFitter) failed!"<<endl;
	//}
	
	/*
	cout<<"*** Initial KFState x0 ***"<<endl;
	fKFInitState.Print();
	cout<<"*** Initial KFCovariance C0 ***"<<endl;
	fKFInitCovariance.Print(); 
	*/


	TMatrixD x0State(nPar,1);
	TMatrixD C0(nPar,nPar);
	TMatrixD dummyStateMatrix(nPar,1);
	TMatrixD dummyCovarianceMatrix(nPar,nPar);

	for(int i=0;i<nPar;i++){
		//x0State(i,0)= gRandom->Uniform(-100,100);
		x0State(i,0)= 0.;
		dummyStateMatrix(i,0)= 0.;
		for(int j=0;j<nPar;j++){
			C0(i,j)= 0.;	
			//C0(i,j)= gRandom->Uniform(1,100);	
			dummyCovarianceMatrix(i,j)= 0.;
		}
		C0(i,i)= 1000.;		
	}

	//cout<<"## Initial KSStateVector x0 ##"<<endl;
	//x0State.Print();
	//cout<<"## Initial KSCovariance C0 ##"<<endl;
	//C0.Print(); 
	
	
	TMatrixD dummyResMatrix(nMeas,nMeas);
	TMatrixD dummyResMatrix2(nMeas,1);
	
	for(int i=0;i<nMeas;i++){
		dummyResMatrix2(i,0)= 0.;
		for(int j=0;j<nMeas;j++){
			dummyResMatrix(i,j)= 0.;	
		}	
	}

	
	//## Init state and covariance matrix vectors
	for(int i=0;i<nIter;i++){
		if(i==0){
			fKFState.push_back(x0State);
			fKFCovariance.push_back(C0);
			//fKFState.push_back(fKFInitState);
			//fKFCovariance.push_back(fKFInitCovariance);
		}	
		else{
			fKFState.push_back(dummyStateMatrix);
			fKFCovariance.push_back(dummyCovarianceMatrix);
		}
		fKFSmoothedState.push_back(dummyStateMatrix);
		fKFSmoothedCovariance.push_back(dummyCovarianceMatrix);
	}

	for(int i=0;i<nHits;i++){
		fKFPredictedState.push_back(dummyStateMatrix);
		fKFPredictedCovariance.push_back(dummyCovarianceMatrix);
		fKFCovarianceMS.push_back(dummyCovarianceMatrix);

		fKFCovarianceResidual.push_back(dummyResMatrix);
		fKFMeasResidual.push_back(dummyResMatrix2);

		fKFSmoothedCovarianceResidual.push_back(dummyResMatrix);
		fKFSmoothedMeasResidual.push_back(dummyResMatrix2);

		fKFPredictedCovarianceResidual.push_back(dummyResMatrix);
		fKFPredictedMeasResidual.push_back(dummyResMatrix2);
	}


	//## Fill measurement, propagator and projector matrix
	//## Measurement and Variance Matrix
	double hits[nHits][3];
	double sigma[nHits][3];
	double lenz[nHits];
	double x,y,z,time;
	double xErr,yErr,zErr,timeErr;
	TMatrixD M(nMeas,1);
	TMatrixD V(nMeas,nMeas);
	for(int i=0; i<nHits; i++) {
  	aTrack->GetPoint(i,x,y,z,time);	
		aTrack->GetPointError(i,xErr,yErr,zErr,timeErr);	
		M(0,0)= x;
		M(1,0)= y;

		V(0,0)= xErr;
		V(1,1)= yErr;
		V(0,1)= 0.;
		V(1,0)= 0.;

		//fKFMeas[i]= M;
		//fKFMeasVariance[i]= V;
		fKFMeas.push_back(M);
		fKFMeasVariance.push_back(V);

		hits[i][0] = x;
		hits[i][1] = y;
		hits[i][2] = z;
			
		sigma[i][0] = xErr;
		sigma[i][1] = yErr;
		sigma[i][2] = zErr;

		lenz[i] = hits[i][2];

		/*
		cout<<"X["<<i<<"]="<< hits[i][0]<<endl;
		cout<<"Y["<<i<<"]="<< hits[i][1]<<endl;
		cout<<"Z["<<i<<"]="<< hits[i][2]<<endl;
		cout<<"Sigma["<<i<<"]="<< sigma[i][0]<<"  "<<sigma[i][1]<<endl;
		*/
	}	

	//## Lenghts Z(k)-Z(k-1)
	for(int i=1;i<nHits;i++) lenz[i] = hits[i][2] - hits[i-1][2];

	//## Scatterer positions
	const int nS= 3;
	Ls[0]= fSuperPlaneDepth[0]-fSuperPlaneSizeZ/2.;
	Ls[1]= fSuperPlaneDepth[1]-fSuperPlaneDepth[0]-fSuperPlaneSizeZ;
	Ls[2]= fSuperPlaneDepth[2]-fSuperPlaneDepth[1]-fSuperPlaneSizeZ;
	ZsStart[0]= 0.;
	ZsStart[1]= -(fSuperPlaneDepth[0]+fSuperPlaneSizeZ/2.);
	ZsStart[2]= -(fSuperPlaneDepth[1]+fSuperPlaneSizeZ/2.);
	for(int i=0;i<nS;i++){
		Zs[i]= -(fSuperPlaneDepth[i]-fSuperPlaneSizeZ/2.);
		//cout<<"Ls["<<i<<"]="<<Ls[i]<<"  Zs["<<i<<"]="<<Zs[i]<<"  ZsStart["<<i<<"]="<<ZsStart[i]<<endl;
	}

	//## Propagator Matrix
	TMatrixD F(nPar,nPar);
	for(int i=0;i<nHits;i++){
		for(int j=0;j<nPar;j++) { 
			for(int k=0;k<nPar;k++) {F(j,k)= 0.;}
			F(j,j)= 1.;
		}
		F(0,1)= 0.;
		F(0,2)= lenz[i];
		F(0,3)= 0.;

		F(1,0)= 0.;
		F(1,2)= 0.;
		F(1,3)= lenz[i];
		
		//fKFPropagator[i]= F;
		fKFPropagator.push_back(F);
	}	

	//## Projector Matrix
	TMatrixD H(nMeas,nPar);
	for(int i=0;i<nHits;i++){
		for(int j=0;j<nMeas;j++){
			for(int k=0;k<nPar;k++) {H(j,k)= 0.;}
		}
		H(0,0)= 1.;
		H(1,1)= 1.;
		//fKFProjector[i]= H;
		fKFProjector.push_back(H);
	}


	//## Printing
	//DebugPrint();

	fFitFCNMin= 0;
	double par[nPar];
	double parErr[nPar];

	
	//## Iterate doing prediction and filter steps
	int nKFIterations= 1;
	for(int s=0;s<nKFIterations;s++){
		if(s>0){
			//## Re-initialize the state and covariance to the smoothed vector of previous iteration
			fKFState[0]= fKFSmoothedState[0];
			fKFCovariance[0]= fKFSmoothedCovariance[0];
		}

		//cout<<"*******************************"<< endl;
		//cout<<"***      KF Iteration --> "<< s << "  ***"<< endl;
		//cout<<"*******************************"<< endl;
		for(int k=1;k<nIter;k++){
			KFPredictionStep(k);
			KFFilterStep(k);
		}
		for(int k=nHits;k>=0;k--){
			KFSmoothingStep(k);
		}

		TMatrixD thisIterationTrackState= fKFSmoothedState[0];
		par[0]= thisIterationTrackState(0,0);
		par[1]= thisIterationTrackState(1,0);
		par[2]= thisIterationTrackState(2,0);
		par[3]= thisIterationTrackState(3,0);

		TMatrixD thisIterationCovariance= fKFSmoothedCovariance[0];
		parErr[0]= sqrt(thisIterationCovariance(0,0));
		parErr[1]= sqrt(thisIterationCovariance(1,1));
		parErr[2]= sqrt(thisIterationCovariance(2,2));
		parErr[3]= sqrt(thisIterationCovariance(3,3));

		TVector3 x0(par[0], par[1], 0.); 
  	TVector3 x1(par[0] + par[2], par[1] + par[3], 1.); 
  	TVector3 direction = (x1-x0).Unit();
		TVector3 oppositedirection = TVector3(-direction.X(),-direction.Y(),-direction.Z());
		TVector3 oppositeZdirection = TVector3(direction.X(),direction.Y(),-direction.Z());
	
		//## track parameters
		fTheta= direction.Theta()*180./TMath::Pi();
		fPhi= direction.Phi()*180./TMath::Pi();

		fVertexX= x0.X();
		fVertexY= x0.Y();	
		fTx= par[2];
		fTy= par[3];

		//## track parameter errors
		fVertexXErr= parErr[0];
		fVertexYErr= parErr[1];
		fTxErr= parErr[2];
		fTyErr= parErr[3];

		//## theta & phi par error
		//## theta= atan2(sqrt(tx^2+ty^2),1)
		//## phi= atan2(ty,tx)
		//## do the error propagation as FCFt, with F derivative matrix d(theta)/d(par_i), d(phi)/d(par_i) 
		TMatrixD DTheta(1,nPar);	
		DTheta(0,0)= 0;
		DTheta(0,1)= 0;	
		DTheta(0,2)= fTx/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
 		DTheta(0,3)= fTy/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
		TMatrixD DThetaT(nPar,1);		
		DThetaT.Transpose(DTheta);

		TMatrixD DPhi(1,nPar);	
		DPhi(0,0)= 0;
		DPhi(0,1)= 0;		
		DPhi(0,2)= -fTy/(fTx*fTx+fTy*fTy);
 		DPhi(0,3)= fTx/(fTx*fTx+fTy*fTy);
		TMatrixD DPhiT(nPar,1);		
		DPhiT.Transpose(DPhi);

		TMatrixD DThetaMultCov(1,nPar);
		DThetaMultCov.Mult(DTheta,thisIterationCovariance);
		TMatrixD DPhiMultCov(1,nPar);
		DPhiMultCov.Mult(DPhi,thisIterationCovariance);
		
		TMatrixD ThetaUncertainty(1,1);
		ThetaUncertainty.Mult(DThetaMultCov,DThetaT);

		TMatrixD PhiUncertainty(1,1);
		PhiUncertainty.Mult(DPhiMultCov,DPhiT);

		fThetaErr= sqrt(ThetaUncertainty(0,0))* 180./TMath::Pi();
		fPhiErr= sqrt(PhiUncertainty(0,0))* 180./TMath::Pi();

		if(s==0){
			//## Store first iter results
			fThetaStart= fTheta;
			fPhiStart= fPhi;
			fTxStart= fTx;
			fTyStart= fTy;
			fVertexXStart= fVertexX;
			fVertexYStart= fVertexY;

			fVertexXStartErr= fVertexXErr;
			fVertexYStartErr= fVertexYErr;
			fTxStartErr= fTxErr;
			fTyStartErr= fTyErr;
			fThetaStartErr= fThetaErr;		
			fPhiStartErr= fPhiErr;
		}

		//cout<<"*** KALMAN FILTER TRACK RECO (ITER "<<s<<") ***"<<endl;
		//cout<<"Theta [deg]="<<fTheta<<" +- "<< fThetaErr << "  Phi [deg]="<<fPhi<<" +- "<<fPhiErr<<endl;
		//cout<<"X0 [cm]="<< fVertexX <<" +- "<< fVertexXErr << "  Y0[cm]="<<fVertexY<<" +- "<< fVertexYErr<<endl;
		//cout<<"Tx="<<fTx<<" +- "<< fTxErr <<"  Ty="<<fTy<<" +- "<< fTyErr <<endl;
		//cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 		//cout<<"*******************************"<<endl;
	}//end iteration
	
	cout<<endl;
	

	//## Printing
	//DebugPrint();

	
	//TMatrixD finalTrackState= fKFState[nIter-1];
	TMatrixD finalTrackState= fKFSmoothedState[0];
	par[0]= finalTrackState(0,0);
	par[1]= finalTrackState(1,0);
	par[2]= finalTrackState(2,0);
	par[3]= finalTrackState(3,0);

	TMatrixD finalTrackCovariance= fKFSmoothedCovariance[0];
	parErr[0]= sqrt(finalTrackCovariance(0,0));
	parErr[1]= sqrt(finalTrackCovariance(1,1));
	parErr[2]= sqrt(finalTrackCovariance(2,2));
	parErr[3]= sqrt(finalTrackCovariance(3,3));

	TVector3 x0(par[0], par[1], 0.); 
  TVector3 x1(par[0] + par[2], par[1] + par[3], 1.); 
  TVector3 direction = (x1-x0).Unit();
	TVector3 oppositedirection = TVector3(-direction.X(),-direction.Y(),-direction.Z());
	TVector3 oppositeZdirection = TVector3(direction.X(),direction.Y(),-direction.Z());
	//fTheta= oppositedirection.Theta()*180./TMath::Pi();
	//fPhi= oppositedirection.Phi()*180./TMath::Pi();

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
	TMatrixD DThetaFinal(1,nPar);	
	DThetaFinal(0,0)= 0;
	DThetaFinal(0,1)= 0;		
	DThetaFinal(0,2)= fTx/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
 	DThetaFinal(0,3)= fTy/(sqrt(fTx*fTx+fTy*fTy)*(1+fTx*fTx+fTy*fTy));
	TMatrixD DThetaFinalT(nPar,1);		
	DThetaFinalT.Transpose(DThetaFinal);
	//cout<<"*** Theta Derivative ***"<<endl;
	//DThetaFinal.Print();
	//cout<<endl;

	TMatrixD DPhiFinal(1,nPar);	
	DPhiFinal(0,0)= 0;
	DPhiFinal(0,1)= 0;		
	DPhiFinal(0,2)= -fTy/(fTx*fTx+fTy*fTy);
 	DPhiFinal(0,3)= fTx/(fTx*fTx+fTy*fTy);
	TMatrixD DPhiFinalT(nPar,1);		
	DPhiFinalT.Transpose(DPhiFinal);
	//cout<<"*** Phi Derivative ***"<<endl;
	//DPhiFinal.Print();
	//cout<<endl;

	TMatrixD DThetaFinalMultCov(1,nPar);
	DThetaFinalMultCov.Mult(DThetaFinal,finalTrackCovariance);
	TMatrixD DPhiFinalMultCov(1,nPar);
	DPhiFinalMultCov.Mult(DPhiFinal,finalTrackCovariance);
		
	//cout<<"*** Final Covariance ***"<<endl;
	//finalTrackCovariance.Print();
	//cout<<endl;

	//cout<<"*** DThetaFinalMultCov ***"<<endl;
	//DThetaFinalMultCov.Print();
	//cout<<endl;

	TMatrixD ThetaFinalUncertainty(1,1);
	ThetaFinalUncertainty.Mult(DThetaFinalMultCov,DThetaFinalT);

	TMatrixD PhiFinalUncertainty(1,1);
	PhiFinalUncertainty.Mult(DPhiFinalMultCov,DPhiFinalT);

	fThetaErr= sqrt(ThetaFinalUncertainty(0,0))* 180./TMath::Pi();
	fPhiErr= sqrt(PhiFinalUncertainty(0,0))* 180./TMath::Pi();

	cout<<"*** KALMAN FILTER TRACK RECO ***"<<endl;
	cout<<"Theta [deg]="<<fTheta<<" +- "<< fThetaErr << "  Phi [deg]="<<fPhi<<" +- "<<fPhiErr<<endl;
	cout<<"X0 [cm]="<< fVertexX <<" +- "<< fVertexXErr << "  Y0[cm]="<<fVertexY<<" +- "<< fVertexYErr<<endl;
	cout<<"Tx="<<fTx<<" +- "<< fTxErr <<"  Ty="<<fTy<<" +- "<< fTyErr <<endl;
	cout<<"fMomentum="<<fMomentum<<endl;
	cout<<"FitFCNMin="<<fFitFCNMin<<endl;
 	cout<<"*******************************"<<endl;


	value = fFitFCNMin;

}//close KFTrackFinder::TrackChi2Fcn()


