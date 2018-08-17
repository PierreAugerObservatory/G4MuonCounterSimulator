/**
* @file G4MuonCounterPrimaryGenerator.cc
* @class G4MuonCounterPrimaryGenerator
* @brief Generate the primary particles in the simulation
* @author S. Riggi
* @date 05/04/2010
*/

#include "G4MuonCounterPrimaryGenerator.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4MuonCounterConstruction.hh"
#include "Randomize.hh"
#include "G4MuonCounterSimulator.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpticalPhoton.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"

#include <det/Detector.h>
#include <sevt/StationSimData.h>
#include <sdet/Station.h>

#include <utl/CoordinateSystem.h>
#include <utl/CoordinateSystemPtr.h>
#include <utl/Point.h>
#include <utl/Particle.h>
#include <utl/ReferenceEllipsoid.h>
#include <utl/AugerCoordinateSystem.h>
#include <utl/AugerUnits.h>

#include <TMath.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TDirectory.h>
#include <TRandom.h>
#include <TRandom3.h>

#include <time.h>

#include <cstddef>
#include <iostream>
#include <sstream>

using namespace G4MuonCounterSimulatorUSC;
//using namespace fwk;
using namespace det;
using namespace sevt;
using utl::CoordinateSystem;
using utl::CoordinateSystemPtr;
using utl::Point;
using utl::Particle;
using utl::ReferenceEllipsoid;
using utl::AugerCoordinateSystem;
using namespace std;


G4MuonCounterPrimaryGenerator::G4MuonCounterPrimaryGenerator(G4MuonCounterConstruction* det)
:Detector(det)
{

  int n_particle = 1;
	fNpartPerEvent= 1;
  fGenTheta=0.0*deg;
	fMinGenTheta= 0.*deg;
	fMaxGenTheta= 90.*deg;
	fGenPhi=0.0*deg;
	fMinGenPhi= 0.*deg;
	fMaxGenPhi= 360.*deg;

	fMinGenEnergy= 0.*GeV;
	fMaxGenEnergy= 1000.*GeV;

	fGenTime= 0.*ns;

  particleGun = new G4ParticleGun(n_particle);

	//set default generated particle
	particleTable = G4ParticleTable::GetParticleTable();

/*
  G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
  fGenParticle= particle;

  particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));//particle from above the strip
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticleTime(0.0*ns);
*/
	
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);
  
  
	fUseShowerGenerator= true;	
	fRandomizePosition= false;
	fRandomizePositionAroundVertex= false;
	fRandomizeDirection= false;
	fRandomizeDirectionUniform= false;
	fRandomizeDirectionWithAcceptance= false;
	fRandomizeEnergy= false;
	fUseCosmicMuonGenerator= false;

	fSizeOfRandomAreaX=0.;
	fSizeOfRandomAreaY=0.;
	fSizeOfRandomAreaZ=0.;
}


G4MuonCounterPrimaryGenerator::G4MuonCounterPrimaryGenerator(){

  int n_particle = 1;
	fNpartPerEvent= 1;
  fGenTheta=0.0*deg;
	fMinGenTheta= 0.*deg;
	fMaxGenTheta= 90.*deg;
	fGenPhi=0.0*deg;
	fMinGenPhi= 0.*deg;
	fMaxGenPhi= 360.*deg;
	fMinGenEnergy= 0.*GeV;
	fMaxGenEnergy= 1000.*GeV;
	fGenTime= 0.*ns;

  particleGun = new G4ParticleGun(n_particle);
  
	//set default generated particle
	particleTable = G4ParticleTable::GetParticleTable();

	/*  
	G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
  fGenParticle= particle;

  particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));//particle from above the strip
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticleTime(0.0*ns);
*/

  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

	fUseShowerGenerator= false;	
	fRandomizePosition= false;
	fRandomizePositionAroundVertex= false;
	fRandomizeDirection= false;
	fRandomizeDirectionWithAcceptance= false;
	fRandomizeEnergy= false;
	fUseCosmicMuonGenerator= false;

	fSizeOfRandomAreaX=0.;
	fSizeOfRandomAreaY=0.;
	fSizeOfRandomAreaZ=0.;
  
}

G4MuonCounterPrimaryGenerator::~G4MuonCounterPrimaryGenerator(){
  delete particleGun;
  delete gunMessenger;
}


void G4MuonCounterPrimaryGenerator::ShowerGenerator(G4Event* anEvent){
	
	const sdet::Station* theCurrentDetectorStation = G4MuonCounterSimulator::GetCurrentDetectorStation();
  sevt::SEvent::StationIterator theCurrentEventStationIt = G4MuonCounterSimulator::GetCurrentEventStationIt();
  sevt::StationSimData::ParticleIterator theCurrentParticleIt = G4MuonCounterSimulator::GetCurrentParticleIt();
  

  // Iterate over the particles in the hit station and pass their
  // properties to the G4ParticleGun:
  
  CoordinateSystemPtr siteCS = det::Detector::GetInstance().GetSiteCoordinateSystem();
  
  const ReferenceEllipsoid& e = ReferenceEllipsoid::GetWGS84();
  Point stationPos = theCurrentDetectorStation->GetPosition();
  
  CoordinateSystemPtr localCS = AugerCoordinateSystem::Create(stationPos, e, siteCS);
  CoordinateSystemPtr csDir = theCurrentParticleIt->GetDirection().GetCoordinateSystem();
  
  // set properties and pass to particle gun      	
  G4ParticleDefinition* particleDef = particleTable->FindParticle(theCurrentParticleIt->GetName());

  if (!particleDef) {
    std::ostringstream msg;
    msg << "Undefined particle type: "
        << theCurrentParticleIt->GetName();
    WARNING(msg);
    return;
  }

  double x, y, z;   		
  boost::tie(x, y, z) = (theCurrentParticleIt->GetPosition() - stationPos).GetCoordinates(localCS);    

	
	//cout<<"particle: "<<theCurrentParticleIt->GetName()<<"  E[MeV]="<< theCurrentParticleIt->GetKineticEnergy()/utl::MeV<< endl;

	double xStat, yStat, zStat;
  boost::tie(xStat, yStat, zStat) = stationPos.GetCoordinates(localCS);    

	//cout<<"StatId="<<theCurrentDetectorStation->GetId()<<"  StatPosSiteCS x="<<stationPos.GetX(siteCS)<<"  y="<<stationPos.GetY(siteCS)<<"  z="<<stationPos.GetZ(siteCS)<<endl;
	//cout<<"StatId="<<theCurrentDetectorStation->GetId()<<"  StatPosLocalCS x="<<stationPos.GetX(localCS)<<"  y="<<stationPos.GetY(localCS)<<"  z="<<stationPos.GetZ(localCS)<<endl;
	//cout<<"(x,y,z) [cm]= ("<< x/utl::cm <<","<<y/utl::cm<<","<<z/utl::cm<<")"<<endl;
	//cout<<"PartPos x="<<theCurrentParticleIt->GetPosition().GetX(localCS)<<"  y="<<theCurrentParticleIt->GetPosition().GetY(localCS)<<"  z="<<theCurrentParticleIt->GetPosition().GetZ(localCS)<<endl;
	//cout<<"(theta,phi) ParticleCS= ("<<theCurrentParticleIt->GetDirection().GetTheta(csDir)/utl::degree<<","<<theCurrentParticleIt->GetDirection().GetPhi(csDir)/utl::degree<<")"<<endl;
	//cout<<"(theta,phi) LocalCS= ("<<theCurrentParticleIt->GetDirection().GetTheta(localCS)/utl::degree<<","<<theCurrentParticleIt->GetDirection().GetPhi(localCS)/utl::degree<<")"<<endl;
	//cout<<"PartDirLocalCS x="<<theCurrentParticleIt->GetDirection().GetX(localCS)<<"  y="<<theCurrentParticleIt->GetDirection().GetY(localCS)<<"  z="<<theCurrentParticleIt->GetDirection().GetZ(localCS)<<endl;


  // Convert from Auger units to G4 units
  G4ThreeVector position(x * CLHEP::m / utl::m, y * CLHEP::m / utl::m, z * CLHEP::m / utl::m);
  
  G4ThreeVector direction(theCurrentParticleIt->GetDirection().GetX(csDir), 
                          theCurrentParticleIt->GetDirection().GetY(csDir), 
                          theCurrentParticleIt->GetDirection().GetZ(csDir));
  
  particleGun->SetParticleDefinition(particleDef);
  particleGun->SetParticlePosition(position);
  particleGun->SetParticleMomentumDirection(direction);

  // Convert from Auger units to G4 units
  particleGun->SetParticleEnergy(theCurrentParticleIt->GetKineticEnergy() * CLHEP::eV / utl::eV);
  particleGun->SetParticleTime(theCurrentParticleIt->GetTime().GetInterval());
    
  particleGun->GeneratePrimaryVertex(anEvent);
	

	
	G4ThreeVector StationPosition(xStat * CLHEP::m / utl::m, yStat * CLHEP::m / utl::m, zStat * CLHEP::m / utl::m);
	G4ThreeVector EdgeSamplingArea1(-1.*CLHEP::m, -1.*CLHEP::m, 0*CLHEP::m);
	G4ThreeVector EdgeSamplingArea2(+1.*CLHEP::m, -1.*CLHEP::m, 0*CLHEP::m);
	G4ThreeVector EdgeSamplingArea3(+1.*CLHEP::m, +1.*CLHEP::m, 0*CLHEP::m);
	G4ThreeVector EdgeSamplingArea4(-1.*CLHEP::m, +1.*CLHEP::m, 0*CLHEP::m);
	
	//## draw origin & vertex
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager){
		//vertex
    G4Circle circle(position);
    circle.SetScreenSize(1.5);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);//red
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);

		//station origin 
		G4Circle originCircle(StationPosition);
    originCircle.SetScreenSize(1.5);
    originCircle.SetFillStyle(G4Circle::filled);
    G4Colour colour2(1.,1.,0.);
    G4VisAttributes attribs2(colour2);
    originCircle.SetVisAttributes(attribs2);
    pVVisManager->Draw(originCircle);

		/*
		//sampling area edge
		G4Circle samplingArea1Circle(EdgeSamplingArea1);
    samplingArea1Circle.SetScreenSize(1.5);
    samplingArea1Circle.SetFillStyle(G4Circle::filled);
    G4Colour colour3(1.,0.,1.);
    G4VisAttributes attribs3(colour3);
    samplingArea1Circle.SetVisAttributes(attribs3);
    pVVisManager->Draw(samplingArea1Circle);
	
		G4Circle samplingArea2Circle(EdgeSamplingArea2);
    samplingArea2Circle.SetScreenSize(1.5);
    samplingArea2Circle.SetFillStyle(G4Circle::filled);
    samplingArea2Circle.SetVisAttributes(attribs3);
    pVVisManager->Draw(samplingArea2Circle);

		G4Circle samplingArea3Circle(EdgeSamplingArea3);
    samplingArea3Circle.SetScreenSize(1.5);
    samplingArea3Circle.SetFillStyle(G4Circle::filled);
    samplingArea3Circle.SetVisAttributes(attribs3);
    pVVisManager->Draw(samplingArea3Circle);

		G4Circle samplingArea4Circle(EdgeSamplingArea4);
    samplingArea4Circle.SetScreenSize(1.5);
    samplingArea4Circle.SetFillStyle(G4Circle::filled);
    samplingArea4Circle.SetVisAttributes(attribs3);
    pVVisManager->Draw(samplingArea4Circle);
		*/
  }
	

	//## estimate expected hit strip from geometrical considerations
	//std::vector<int> ExpectedHitStripList= CalculateGeomHits(position,direction);

}//close ShowerGenerator()


void G4MuonCounterPrimaryGenerator::ShowerGeneratorList(G4Event* anEvent){
	
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	
	const sdet::Station* theCurrentDetectorStation = G4MuonCounterSimulator::GetCurrentDetectorStation();
  sevt::SEvent::StationIterator theCurrentEventStationIt = G4MuonCounterSimulator::GetCurrentEventStationIt();
  std::vector<sevt::StationSimData::ParticleIterator> theCurrentParticleItList = G4MuonCounterSimulator::GetCurrentParticleItList();
  

  // Iterate over the particles in the hit station and pass their
  // properties to the G4ParticleGun
  CoordinateSystemPtr siteCS = det::Detector::GetInstance().GetSiteCoordinateSystem();
  
  const ReferenceEllipsoid& e = ReferenceEllipsoid::GetWGS84();
  Point stationPos = theCurrentDetectorStation->GetPosition();
  
  CoordinateSystemPtr localCS = AugerCoordinateSystem::Create(stationPos, e, siteCS);
  
	
  // set properties and pass to particle gun      	
	for(unsigned int k=0;k<theCurrentParticleItList.size();k++){  
		CoordinateSystemPtr csDir = theCurrentParticleItList[k]->GetDirection().GetCoordinateSystem();
  
		G4ParticleDefinition* particleDef = particleTable->FindParticle(theCurrentParticleItList[k]->GetName());

  	if (!particleDef) {
    	std::ostringstream msg;
    	msg << "Undefined particle type: "
          << theCurrentParticleItList[k]->GetName();
   		 WARNING(msg);
    	return;
  	}

  	double x, y, z;   		
  	boost::tie(x, y, z) = (theCurrentParticleItList[k]->GetPosition() - stationPos).GetCoordinates(localCS);    

		cout<<"particle: "<<theCurrentParticleItList[k]->GetName()<<"  E[MeV]="<< theCurrentParticleItList[k]->GetKineticEnergy()/utl::MeV<< endl;

		double xStat, yStat, zStat;
  	boost::tie(xStat, yStat, zStat) = stationPos.GetCoordinates(localCS);    


  	// Convert from Auger units to G4 units
  	G4ThreeVector position(x * CLHEP::m / utl::m, y * CLHEP::m / utl::m, z * CLHEP::m / utl::m);
  
  	G4ThreeVector direction(theCurrentParticleItList[k]->GetDirection().GetX(csDir), 
                          	theCurrentParticleItList[k]->GetDirection().GetY(csDir), 
                          	theCurrentParticleItList[k]->GetDirection().GetZ(csDir));
  
  	particleGun->SetParticleDefinition(particleDef);
  	particleGun->SetParticlePosition(position);
  	particleGun->SetParticleMomentumDirection(direction);

 	 	// Convert from Auger units to G4 units
  	particleGun->SetParticleEnergy(theCurrentParticleItList[k]->GetKineticEnergy() * CLHEP::eV / utl::eV);
  	particleGun->SetParticleTime(theCurrentParticleItList[k]->GetTime().GetInterval());
    
  	particleGun->GeneratePrimaryVertex(anEvent);
	
		
		//## draw vertex		
  	if(pVVisManager){
			//vertex
    	G4Circle circle(position);
    	circle.SetScreenSize(1.5);
    	circle.SetFillStyle(G4Circle::filled);
    	G4Colour colour(1.,0.,0.);//red
    	G4VisAttributes attribs(colour);
    	circle.SetVisAttributes(attribs);
    	pVVisManager->Draw(circle);		
  	}

	}//end loop list of particles

}//close ShowerGeneratorList()



void G4MuonCounterPrimaryGenerator::DebugGenerator(G4Event* anEvent){
	
	sevt::StationSimData::ParticleIterator theCurrentParticleIt = G4MuonCounterSimulator::GetCurrentParticleIt();

	/*
	G4ParticleDefinition* particleDef = particleTable->FindParticle(theCurrentParticleIt->GetName());
  if (!particleDef) {
    std::ostringstream msg;
    msg << "Undefined particle type: "
        << theCurrentParticleIt->GetName();
    WARNING(msg);
    return;
  }
	*/

	G4ParticleDefinition* particleDef = particleTable->FindParticle("mu-");
	
	//particle: mu-  E[MeV]=485.364
	//(x,y,z) [cm]= (-37.9112,18.7625,5.53882e-08)
	//(dx,dy,dz)= (-0.938992,-0.184035,-0.290561)
	double x= -37.9112 * utl::cm;
	double y= 18.7625 * utl::cm;	
	double z= 5.53882e-08* utl::cm;

  //G4ThreeVector position(x * CLHEP::m, 0 * CLHEP::m, 0 * CLHEP::m);
	G4ThreeVector position(x * CLHEP::m / utl::m, y * CLHEP::m / utl::m, z * CLHEP::m / utl::m);
  G4ThreeVector direction(-0.938992, -0.184035, -0.290561);
	double energy= 485.364*utl::MeV;
	double time= 0.*CLHEP::ns;
  //cout<<"particle E[MeV]="<< energy/CLHEP::MeV<< endl;
	//cout<<"particle: "<<theCurrentParticleIt->GetName()<<"  E[MeV]="<< theCurrentParticleIt->GetKineticEnergy() /utl::MeV<< endl;

  particleGun->SetParticleDefinition(particleDef);
  particleGun->SetParticlePosition(position);
  particleGun->SetParticleMomentumDirection(direction);

  // Convert from Auger units to G4 units
  //particleGun->SetParticleEnergy(theCurrentParticleIt->GetKineticEnergy() * CLHEP::eV / utl::eV);
  particleGun->SetParticleTime(theCurrentParticleIt->GetTime().GetInterval());
	particleGun->SetParticleEnergy(energy);
  //particleGun->SetParticleTime(time);
    
  particleGun->GeneratePrimaryVertex(anEvent);

	//## draw origin & vertex
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager){
		//vertex
    G4Circle circle(position);
    circle.SetScreenSize(1.5);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);//red
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
	}


}//close DebugGenerator()



void G4MuonCounterPrimaryGenerator::GeneratePrimaries(G4Event* anEvent){

	//DebugGenerator(anEvent);
	
	(G4MuonCounterSimulator::fOneParticlePerStation) ?
		ShowerGenerator(anEvent) :
		ShowerGeneratorList(anEvent);

	
}//close G4MuonCounterPrimaryGenerator::GeneratePrimaries()


void G4MuonCounterPrimaryGenerator::SetOptPhotonPolar(){
 double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}


void G4MuonCounterPrimaryGenerator::SetOptPhotonPolar(double angle)
{
 if (particleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton")
   {
     cout << "--> warning from G4MuonCounterPrimaryGenerator::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << endl;
     return;
   }
     	       
 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = particleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton); 
 double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product; 
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 particleGun->SetParticlePolarization(polar);
}


double G4MuonCounterPrimaryGenerator::GetEnergyFromSr90Spectrum(){

	double GeneratedEnergy;

	const int n= 77;
	//sampled energy in keV
	double SampledEnergy[]={141*keV,163.8*keV,175.2*keV,179*keV,190.5*keV,201.9*keV,209.5*keV,217.1*keV,228.6*keV,236.2*keV,255.2*keV,285.7*keV,297.1*keV,323.8*keV,342.9*keV,373.3*keV,384.8*keV,403.8*keV,411.4*keV,430.5*keV,441.9*keV,468.6*keV,487.6*keV,506.7*keV,525.7*keV,560*keV,590.5*keV,613.3*keV,621*keV,628.6*keV,636.2*keV,659*keV,678.1*keV,704.8*keV,731.4*keV,769.5*keV,796.2*keV,857.1*keV,880*keV,918.1*keV,929.5*keV,960*keV,975.2*keV,998.1*keV,1017*keV,1048*keV,1074*keV,1105*keV,1120*keV,1143*keV,1162*keV,1177*keV,1200*keV,1280*keV,1314*keV,1341*keV,1371*keV,1394*keV,1455*keV,1482*keV,1512*keV,1589*keV,1600*keV,1627*keV,1642*keV,1691*keV,1703*keV,1722*keV,1741*keV,1825*keV,1840*keV,1943*keV,1973*keV,1981*keV,2107*keV,2168*keV,2206*keV,2213*keV};
	
  double SampledSpectrum[]={39.09,187.2,351.9,491.8,738.7,1158,1422,1825,2130,2467,2862,3093,3224,3265,3298,3241,3126,2994,2846,2747,2599,2352,2121,1973,1751,1611,1488,1340,1323,1348,1381,1389,1381,1471,1356,1471,1405,1447,1447,1422,1389,1463,1397,1356,1414,1455,1372,1381,1356,1438,1331,1381,1307,1290,1233,1290,1126,1200,1043,1142,1010,928,977.4,862.1,895.1,829.2,746.9,755.1,672.8,549.4,565.8,376.5,310.7,343.6,96.71,22.63,14.4,0.};	
	
	double Energy,Spectrum;

	TTree* SpectrumTree=new TTree();
	SpectrumTree->Branch("Energy",&Energy,"Energy/D");
	SpectrumTree->Branch("Spectrum",&Spectrum,"Spectrum/D");

	for(int i=0;i<n;i++){
		Energy= SampledEnergy[i];//express in KeV
		Spectrum= SampledSpectrum[i];
		SpectrumTree->Fill();
	}//close for
	
	SpectrumTree->Draw("Spectrum:Energy >> SpectrumHist2D", "", "goff");
	TH2D* SpectrumHist2D= (TH2D*)gDirectory->Get("SpectrumHist2D"); 
	TProfile* SpectrumHistProfile= (TProfile*)SpectrumHist2D->ProfileX();
  TH1D* SpectrumHist=(TH1D*)SpectrumHistProfile->ProjectionX();

	GeneratedEnergy= SpectrumHist->GetRandom();

	cout<<"GeneratedEnergy="<<GeneratedEnergy/MeV<<endl;
	
	//clean memory
	SpectrumTree->Delete();
	SpectrumHist2D->Delete();
	SpectrumHistProfile->Delete();
	SpectrumHist->Delete();

	return GeneratedEnergy;

}//close function



