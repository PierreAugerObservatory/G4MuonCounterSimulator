/**
* @file Utilities.cc
* @class Utilities
* @brief Utility functions for 
*
* Useful functions to be used in the simulation, or to handle output files
* @author Dr. Simone Riggi
* @date 23/08/2010
*/

#include "Utilities.hh"
#include "TScintHit.hh"
#include "TPMTHit.hh"
#include "AnalysisConsts.hh"

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

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

Utilities::Utilities(){

}

Utilities::~Utilities(){

}


void Utilities::MergeFiles(const char * fileListName){

	ifstream inputfilelist;
  inputfilelist.open(fileListName);
  char filename[400];
	char emptyline[]="";

	TFile* outputFile= new TFile("MergedFile.root","RECREATE");

	//*******************************
	//    OUTPUT MERGED FILE
	//*******************************
	double X0, Y0, Z0, Px, Py, Pz, KinEnergy, Mass, Time, Theta, Phi;
	int PDGCode;
	double Xstrip,Ystrip,Zstrip,Xcoat,Ycoat,Zcoat,XYPlaneDistance;
	int Nstrips, Nplanes;
	const int MAXNOPLANES= 100;
	double PlaneDistance[MAXNOPLANES];
	double PlaneDepth[MAXNOPLANES];
	
	std::vector<TScintHit>* fScintHitCollection;
	std::vector<TPMTHit>* fPMTHitCollection;

	TTree* PMTTree_mergedFile= new TTree("SimTree_pmt","SimTree_pmt");
	TTree* ScintTree_mergedFile= new TTree("SimTree_scint","SimTree_scint");
	TTree* GenTree_mergedFile= new TTree("SimTree_gen","SimTree_gen");
	TTree* DetTree_mergedFile= new TTree("SimTree_det","SimTree_det");

  
	//## DET INFO
	DetTree_mergedFile->Branch("Nstrips",&Nstrips,"Nstrips/I");
	DetTree_mergedFile->Branch("Nplanes",&Nplanes,"Nplanes/I");
	DetTree_mergedFile->Branch("XYPlaneDistance",&XYPlaneDistance,"XYPlaneDistance/D");
	DetTree_mergedFile->Branch("PlaneDistance",PlaneDistance,"PlaneDistance[Nplanes]/D");
	DetTree_mergedFile->Branch("PlaneDepth",PlaneDepth,"PlaneDepth[Nplanes]/D");
	DetTree_mergedFile->Branch("Xstrip",&Xstrip,"Xstrip/D");
	DetTree_mergedFile->Branch("Ystrip",&Ystrip,"Ystrip/D");
	DetTree_mergedFile->Branch("Zstrip",&Zstrip,"Zstrip/D");
	DetTree_mergedFile->Branch("Xcoat",&Xcoat,"Xcoat/D");
	DetTree_mergedFile->Branch("Ycoat",&Ycoat,"Ycoat/D");
	DetTree_mergedFile->Branch("Zcoat",&Zcoat,"Zcoat/D");

	//## GEN INFO
  GenTree_mergedFile->Branch("X0",&X0,"X0/D");
  GenTree_mergedFile->Branch("Y0",&Y0,"Y0/D");
  GenTree_mergedFile->Branch("Z0",&Z0,"Z0/D");
	GenTree_mergedFile->Branch("Px",&Px,"Px/D");
  GenTree_mergedFile->Branch("Py",&Py,"Py/D");
  GenTree_mergedFile->Branch("Pz",&Pz,"Pz/D");
	GenTree_mergedFile->Branch("KinEnergy",&KinEnergy,"KinEnergy/D");
	GenTree_mergedFile->Branch("Mass",&Mass,"Mass/D");
	GenTree_mergedFile->Branch("Time",&Time,"Time/D");
	GenTree_mergedFile->Branch("PDGCode",&PDGCode,"PDGCode/I");
	GenTree_mergedFile->Branch("Theta",&Theta,"Theta/D");
	GenTree_mergedFile->Branch("Phi",&Phi,"Phi/D");
		

	//## PMT
  PMTTree_mergedFile->Branch("PMTHitCollection",&fPMTHitCollection);
 
	
  //## SCINT 
	ScintTree_mergedFile->Branch("ScintHitCollection",&fScintHitCollection);

	TFile* inputFile;
	TTree* SimTree_pmt;
	TTree* SimTree_scint;
	TTree* SimTree_gen;	
	TTree* SimTree_det;

	while(!inputfilelist.eof()){

    if (inputfilelist.good() ) {
      inputfilelist.getline(filename, 400);
      if ( inputfilelist.eof() ) break;
      if(strcmp(filename,emptyline) == 0) break;
      else cout<<"### Reading file "<<filename<<" ###"<<endl;

    }
    else break;

		inputFile=new TFile(filename);

		SimTree_pmt= (TTree*)inputFile->Get("SimTree_pmt");
  	SimTree_scint=(TTree*)inputFile->Get("SimTree_scint");
  	SimTree_gen=(TTree*)inputFile->Get("SimTree_gen");
		SimTree_det=(TTree*)inputFile->Get("SimTree_det");

		//## det tree
		SimTree_det->SetBranchAddress("Nstrips",&Nstrips);
		SimTree_det->SetBranchAddress("Nplanes",&Nplanes);
		SimTree_det->SetBranchAddress("XYPlaneDistance",&XYPlaneDistance);
		SimTree_det->SetBranchAddress("PlaneDistance",PlaneDistance);
		SimTree_det->SetBranchAddress("PlaneDepth",PlaneDepth);
		SimTree_det->SetBranchAddress("Xstrip",&Xstrip);
		SimTree_det->SetBranchAddress("Ystrip",&Ystrip);
		SimTree_det->SetBranchAddress("Zstrip",&Zstrip);
		SimTree_det->SetBranchAddress("Xcoat",&Xcoat);
		SimTree_det->SetBranchAddress("Ycoat",&Ycoat);
		SimTree_det->SetBranchAddress("Zcoat",&Zcoat);
	

		//## gen tree
  	SimTree_gen->SetBranchAddress("X0",&X0);
  	SimTree_gen->SetBranchAddress("Y0",&Y0);
  	SimTree_gen->SetBranchAddress("Z0",&Z0);
		SimTree_gen->SetBranchAddress("Px",&Px);
  	SimTree_gen->SetBranchAddress("Py",&Py);
  	SimTree_gen->SetBranchAddress("Pz",&Pz);
		SimTree_gen->SetBranchAddress("KinEnergy",&KinEnergy);
		SimTree_gen->SetBranchAddress("Mass",&Mass);
		SimTree_gen->SetBranchAddress("Time",&Time);
 	 	SimTree_gen->SetBranchAddress("PDGCode",&PDGCode);
  	SimTree_gen->SetBranchAddress("Theta",&Theta);
		SimTree_gen->SetBranchAddress("Phi",&Phi);	

		//initialize vectors
		fScintHitCollection = 0;
		fPMTHitCollection = 0;  

		//## pmt tree
		SimTree_pmt->SetBranchAddress("PMTHitCollection",&fPMTHitCollection); 
	
		//## scintillator tree
		SimTree_scint->SetBranchAddress("ScintHitCollection",&fScintHitCollection);
	
		for(int i=0;i<SimTree_scint->GetEntries();i++){
			SimTree_scint->GetEntry(i);
			SimTree_pmt->GetEntry(i);
			SimTree_gen->GetEntry(i);
			SimTree_det->GetEntry(i);
	
			ScintTree_mergedFile->Fill();
			PMTTree_mergedFile->Fill();
			GenTree_mergedFile->Fill();
			DetTree_mergedFile->Fill();
		}

	}//close while

	
	outputFile->cd();
	ScintTree_mergedFile->Write();
	PMTTree_mergedFile->Write();
	GenTree_mergedFile->Write();
	DetTree_mergedFile->Write();


}//close Utilities::MergeFiles



std::string Utilities::ExecSystemCommand(const char* cmd) {
	
	FILE* pipe = popen(cmd, "r");
  if (!pipe) {	
		cerr<<"ERROR: Cannot get exec input"<<endl;
		exit(1);
	}

  char buffer[128];
  std::string parsedCmdOutput = "";
  while(!feof(pipe)) {
 		if(fgets(buffer, 128, pipe) != NULL)
    	parsedCmdOutput += buffer;
  }
  pclose(pipe);

	return parsedCmdOutput;

}//close Utilities::ExecSystemCommand


