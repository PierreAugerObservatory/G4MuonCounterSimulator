/**
* @file PrimaryGeneratorMessenger.cc
* @class PrimaryGeneratorMessenger
* @brief Create the UI commands to handle the primary generator parameters, i.e. from a configuration file 
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "PrimaryGeneratorMessenger.hh"
#include "G4MuonCounterPrimaryGenerator.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"


using namespace G4MuonCounterSimulatorUSC;

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(G4MuonCounterPrimaryGenerator* Gun)
:Generator(Gun)
{
  genDir = new G4UIdirectory("/StripSimulation/gen/");
  genDir->SetGuidance("PrimaryGenerator control");
   
  polarCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/optPhotonPolar",this);
  polarCmd->SetGuidance("Set linear polarization");
  polarCmd->SetGuidance("  angle w.r.t. (k,n) plane");
  polarCmd->SetParameterName("angle",true);
  polarCmd->SetUnitCategory("Angle");  
  polarCmd->SetDefaultValue(-360.0);
  polarCmd->SetDefaultUnit("deg");
  polarCmd->AvailableForStates(G4State_Idle);

	genPositionCmd = new G4UIcmdWith3VectorAndUnit("/StripSimulation/gen/position",this);
  genPositionCmd->SetGuidance("Set the position of the generated particles.");
  genPositionCmd->SetParameterName("X","Y","Z",true,true);
  genPositionCmd->SetDefaultUnit("cm");
	genPositionCmd->SetUnitCandidates("micron mm cm m km");
	genPositionCmd->SetDefaultValue(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));

	genDirectionCmd = new G4UIcmdWith3Vector("/StripSimulation/gen/direction",this);
  genDirectionCmd->SetGuidance("Set the direction of the generated particles.");
  genDirectionCmd->SetParameterName("Px","Py","Pz",true,true);
	genDirectionCmd->SetDefaultValue(G4ThreeVector(0.,0.,-1.));
	genDirectionCmd->SetRange("Px != 0 || Py != 0 || Pz != 0");

  genThetaCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/theta",this);
  genThetaCmd->SetGuidance("Set theta of the generated particles");
  genThetaCmd->SetParameterName("theta",true,true);
	genThetaCmd->SetDefaultUnit("deg");
  genThetaCmd->SetUnitCandidates("deg rad");
	genThetaCmd->SetDefaultValue(0.0);

	genMinThetaCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/minTheta",this);
  genMinThetaCmd->SetGuidance("Set minimum theta of the generated particles");
  genMinThetaCmd->SetParameterName("minTheta",true,true);
	genMinThetaCmd->SetDefaultUnit("deg");
  genMinThetaCmd->SetUnitCandidates("deg rad");
	genMinThetaCmd->SetDefaultValue(0.0*deg);

	genMaxThetaCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/maxTheta",this);
  genMaxThetaCmd->SetGuidance("Set maximum theta of the generated particles");
  genMaxThetaCmd->SetParameterName("maxTheta",true,true);
	genMaxThetaCmd->SetDefaultUnit("deg");
  genMaxThetaCmd->SetUnitCandidates("deg rad");
	genMaxThetaCmd->SetDefaultValue(90.0*deg);


	genPhiCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/phi",this);
  genPhiCmd->SetGuidance("Set phi of the generated particles");
  genPhiCmd->SetParameterName("phi",true,true);
	genPhiCmd->SetDefaultUnit("deg");
  genPhiCmd->SetUnitCandidates("deg rad");
	genPhiCmd->SetDefaultValue(0.0);

	genMinPhiCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/minPhi",this);
  genMinPhiCmd->SetGuidance("Set minimum phi of the generated particles");
  genMinPhiCmd->SetParameterName("minPhi",true,true);
	genMinPhiCmd->SetDefaultUnit("deg");
  genMinPhiCmd->SetUnitCandidates("deg rad");
	genMinPhiCmd->SetDefaultValue(0.0*deg);

	genMaxPhiCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/maxPhi",this);
  genMaxPhiCmd->SetGuidance("Set maximum phi of the generated particles");
  genMaxPhiCmd->SetParameterName("maxPhi",true,true);
	genMaxPhiCmd->SetDefaultUnit("deg");
  genMaxPhiCmd->SetUnitCandidates("deg rad");
	genMaxPhiCmd->SetDefaultValue(0.0*deg);


	genEnergyCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/energy",this);
  genEnergyCmd->SetGuidance("Set energy of the generated particles");
  genEnergyCmd->SetParameterName("energy",true,true);
	genEnergyCmd->SetDefaultUnit("GeV");
  genEnergyCmd->SetUnitCandidates("eV keV MeV GeV TeV");
	genEnergyCmd->SetDefaultValue(1.0);

	genMinEnergyCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/minEnergy",this);
  genMinEnergyCmd->SetGuidance("Set the minimum energy of the generated particles");
  genMinEnergyCmd->SetParameterName("minEnergy",true,true);
	genMinEnergyCmd->SetDefaultUnit("GeV");
  genMinEnergyCmd->SetUnitCandidates("eV keV MeV GeV TeV");
	genMinEnergyCmd->SetDefaultValue(0.1);

	genMaxEnergyCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/maxEnergy",this);
  genMaxEnergyCmd->SetGuidance("Set the maximum energy of the generated particles");
  genMaxEnergyCmd->SetParameterName("maxEnergy",true,true);
	genMaxEnergyCmd->SetDefaultUnit("GeV");
  genMaxEnergyCmd->SetUnitCandidates("eV keV MeV GeV TeV");
	genMaxEnergyCmd->SetDefaultValue(1000.);
  
	genTimeCmd = new G4UIcmdWithADoubleAndUnit("/StripSimulation/gen/time",this);
  genTimeCmd->SetGuidance("Set time of the generated particles");
  genTimeCmd->SetParameterName("time",true,true);
	genTimeCmd->SetDefaultUnit("ns");
  genTimeCmd->SetUnitCandidates("ns ms s");
	genTimeCmd->SetDefaultValue(0.0);

	genParticleCmd = new G4UIcmdWithAString("/StripSimulation/gen/particle",this);
  genParticleCmd->SetGuidance("Set particle to be generated.");
  genParticleCmd->SetGuidance(" (geantino is default)");
  genParticleCmd->SetParameterName("particleName",true);
  genParticleCmd->SetDefaultValue("geantino");
  G4String candidateList; 
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4int nPtcl = particleTable->entries();
  for(G4int i=0;i<nPtcl;i++){
    candidateList += particleTable->GetParticleName(i);
    candidateList += " ";
  }
  genParticleCmd->SetCandidates(candidateList);


	genNumberOfParticlesCmd = new G4UIcmdWithAnInteger("/StripSimulation/gen/nparticle",this);
  genNumberOfParticlesCmd->SetGuidance("Set the number of particles injected per event");
  genNumberOfParticlesCmd->SetParameterName("nparticle",true); 
  genNumberOfParticlesCmd->SetDefaultValue(1);
  genNumberOfParticlesCmd->AvailableForStates(G4State_Idle);

	//genParticlesFromListCmd = new G4UIcmdWithABool("/StripSimulation/gen/particleFromList",this);
  //genParticlesFromListCmd->SetGuidance("Enable/Disable the injection of particles from a list (e.g. as set from external generator)");
  //genParticlesFromListCmd->SetParameterName("particleFromList",true,true);
	//genParticlesFromListCmd->SetDefaultValue(false);

	genParticlesFromShowerCmd = new G4UIcmdWithABool("/StripSimulation/gen/particleFromShower",this);
  genParticlesFromShowerCmd->SetGuidance("Enable/Disable the injection of particles from resampled showers");
  genParticlesFromShowerCmd->SetParameterName("particleFromShower",true,true);
	genParticlesFromShowerCmd->SetDefaultValue(false);


  genRandomDirectionCmd = new G4UIcmdWithABool("/StripSimulation/gen/randomDirection",this);
  genRandomDirectionCmd->SetGuidance("Enable/Disable the random generation of particle directions");
  genRandomDirectionCmd->SetParameterName("randomDirection",true,true);
  genRandomDirectionCmd->SetDefaultValue(false);

	genRandomDirectionUniformCmd = new G4UIcmdWithABool("/StripSimulation/gen/randomDirectionUniform",this);
  genRandomDirectionUniformCmd->SetGuidance("Enable/Disable the random generation of particle directions (uniform in space)");
  genRandomDirectionUniformCmd->SetParameterName("randomDirectionUniform",true,true);
  genRandomDirectionUniformCmd->SetDefaultValue(false);


	genRandomDirectionWithAcceptanceCmd = new G4UIcmdWithABool("/StripSimulation/gen/randomDirectionWithAcceptance",this);
  genRandomDirectionWithAcceptanceCmd->SetGuidance("Enable/Disable the random generation of particle directions (according to acceptance)");
  genRandomDirectionWithAcceptanceCmd->SetParameterName("randomDirectionWithAcceptance",true,true);
	genRandomDirectionWithAcceptanceCmd->SetDefaultValue(false);

	genRandomPositionCmd = new G4UIcmdWithABool("/StripSimulation/gen/randomPosition",this);
  genRandomPositionCmd->SetGuidance("Enable/Disable the random generation of particle position (in a 3D grid around the detector)");
  genRandomPositionCmd->SetParameterName("randomPosition",true,true);
	genRandomPositionCmd->SetDefaultValue(false);

	genRandomPositionAroundVertexCmd = new G4UIcmdWithABool("/StripSimulation/gen/randomPositionAroundVertex",this);
  genRandomPositionAroundVertexCmd->SetGuidance("Enable/Disable the random generation of particle position (in a 3D grid around the chosen vertex)");
  genRandomPositionAroundVertexCmd->SetParameterName("randomPositionAroundVertex",true,true);
	genRandomPositionAroundVertexCmd->SetDefaultValue(false);

	genRandomPositionWithAcceptanceCmd = new G4UIcmdWithABool("/StripSimulation/gen/randomPositionWithAcceptance",this);
  genRandomPositionWithAcceptanceCmd->SetGuidance("Enable/Disable the random generation of particle position (according to acceptance)");
  genRandomPositionWithAcceptanceCmd->SetParameterName("randomPositionWithAcceptance",true,true);
	genRandomPositionWithAcceptanceCmd->SetDefaultValue(false);

	genCosmicMuonCmd = new G4UIcmdWithABool("/StripSimulation/gen/genCosmicMuon",this);
  genCosmicMuonCmd->SetGuidance("Enable/Disable the random generation of cosmic muons");
  genCosmicMuonCmd->SetParameterName("genCosmicMuon",true,true);
	genCosmicMuonCmd->SetDefaultValue(false);

	genRandomEnergyCmd = new G4UIcmdWithABool("/StripSimulation/gen/randomEnergy",this);
  genRandomEnergyCmd->SetGuidance("Enable/Disable the random generation of energy");
  genRandomEnergyCmd->SetParameterName("randomEnergy",true,true);
	genRandomEnergyCmd->SetDefaultValue(false);

	genRandomGridSizeCmd = new G4UIcmdWith3VectorAndUnit("/StripSimulation/gen/randomGridSize",this);
  genRandomGridSizeCmd->SetGuidance("Set the size of the 3D grid used to randomize particle vertexes");
  genRandomGridSizeCmd->SetParameterName("Xrand","Yrand","Zrand",true,true);
  genRandomGridSizeCmd->SetDefaultUnit("cm");
	genRandomGridSizeCmd->SetUnitCandidates("micron mm cm m km");
	genRandomGridSizeCmd->SetDefaultValue(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));

}


PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
	delete genDir;
  delete polarCmd;
	delete genPositionCmd;
	delete genDirectionCmd;
	delete genEnergyCmd;
	delete genMinEnergyCmd;
	delete genMaxEnergyCmd;
  delete genThetaCmd;
	delete genMinThetaCmd;
	delete genMaxThetaCmd;
	delete genPhiCmd;
	delete genMinPhiCmd;
	delete genMaxPhiCmd;
	delete genTimeCmd;
	delete genParticleCmd;
	delete genNumberOfParticlesCmd;
  delete genCosmicMuonCmd;
	delete genParticlesFromShowerCmd;
	//delete genParticlesFromListCmd;
	delete genRandomEnergyCmd;
	delete genRandomDirectionCmd;
	delete genRandomDirectionUniformCmd;
	delete genRandomDirectionWithAcceptanceCmd;
	delete genRandomPositionCmd;
	delete genRandomPositionAroundVertexCmd;
	delete genRandomPositionWithAcceptanceCmd;
	delete genRandomGridSizeCmd;
	
}


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
 
  if( command == polarCmd ) {
  	G4double angle = polarCmd->GetNewDoubleValue(newValue);
    if(angle == -360.0*deg ) {
    	Generator->SetOptPhotonPolar();
    } 
		else Generator->SetOptPhotonPolar(angle);      
  }
	
	else if(command==genParticleCmd){
		G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pd = particleTable->FindParticle(newValue);
    if(pd!=NULL) Generator->SetGenParticle(pd);
  }
  else if(command==genDirectionCmd){ 
		Generator->SetGenDirection(genDirectionCmd->GetNew3VectorValue(newValue)); 
	}
  else if(command==genThetaCmd){ 
    Generator->SetThetaOfGenParticles(genThetaCmd->GetNewDoubleValue(newValue));
	}
	else if(command==genMinThetaCmd){ 
    Generator->SetMinThetaOfGenParticles(genMinThetaCmd->GetNewDoubleValue(newValue));
	}
	else if(command==genMaxThetaCmd){ 
    Generator->SetMaxThetaOfGenParticles(genMaxThetaCmd->GetNewDoubleValue(newValue));
	}
	else if(command==genPhiCmd){ 
    Generator->SetPhiOfGenParticles(genPhiCmd->GetNewDoubleValue(newValue));
	}
	else if(command==genMinPhiCmd){ 
    Generator->SetMinPhiOfGenParticles(genMinPhiCmd->GetNewDoubleValue(newValue));
	}
	else if(command==genMaxPhiCmd){ 
    Generator->SetMaxPhiOfGenParticles(genMaxPhiCmd->GetNewDoubleValue(newValue));
	}
  else if(command==genEnergyCmd){ 
		Generator->SetGenEnergy(genEnergyCmd->GetNewDoubleValue(newValue)); 
	}
	else if(command==genMinEnergyCmd){ 
		Generator->SetMinGenEnergy(genMinEnergyCmd->GetNewDoubleValue(newValue)); 
	}
	else if(command==genMaxEnergyCmd){ 
		Generator->SetMaxGenEnergy(genMaxEnergyCmd->GetNewDoubleValue(newValue)); 
	}
  else if(command==genPositionCmd){ 
		Generator->SetGenPosition(genPositionCmd->GetNew3VectorValue(newValue)); 
	}
  else if(command==genTimeCmd){ 
		Generator->SetGenTime(genTimeCmd->GetNewDoubleValue(newValue)); 
	}
	else if(command==genNumberOfParticlesCmd) {
		G4int NpartPerEvent = genNumberOfParticlesCmd->GetNewIntValue(newValue);
		Generator->SetNumberOfGenParticles(NpartPerEvent);
	}
	/*
	else if(command==genParticlesFromListCmd) {
		Generator->SetParticleInjectionFromList(genParticlesFromListCmd->GetNewBoolValue(newValue));
	}
	*/
	else if(command==genParticlesFromShowerCmd) {
		Generator->SetParticleInjectionFromResampledShower(genParticlesFromShowerCmd->GetNewBoolValue(newValue));
	}
	else if(command==genRandomDirectionCmd) {
		Generator->SetRandomizingDirection(genRandomDirectionCmd->GetNewBoolValue(newValue));
	}
	else if(command==genRandomDirectionUniformCmd) {
		Generator->SetRandomizingDirectionUniform(genRandomDirectionUniformCmd->GetNewBoolValue(newValue));
	}
	else if(command==genRandomDirectionWithAcceptanceCmd) {
		Generator->SetRandomizingDirectionWithAcceptance(genRandomDirectionWithAcceptanceCmd->GetNewBoolValue(newValue));
	}
	else if(command==genRandomPositionCmd) {
		Generator->SetRandomizingPosition(genRandomPositionCmd->GetNewBoolValue(newValue));
	}	
	else if(command==genRandomPositionAroundVertexCmd) {
		Generator->SetRandomizingPositionAroundVertex(genRandomPositionAroundVertexCmd->GetNewBoolValue(newValue));
	}
	else if(command==genRandomPositionWithAcceptanceCmd) {
		Generator->SetRandomizingPositionWithAcceptance(genRandomPositionWithAcceptanceCmd->GetNewBoolValue(newValue));
	}
	else if(command==genRandomGridSizeCmd){
		Generator->SetSizeOfRandomInjectionArea(genRandomGridSizeCmd->GetNew3VectorValue(newValue)); 
	}
	else if(command==genCosmicMuonCmd){
		Generator->UseCosmicMuonGenerator(genCosmicMuonCmd->GetNewBoolValue(newValue)); 
	}
	else if(command==genRandomEnergyCmd){
		Generator->SetRandomizingEnergy(genRandomEnergyCmd->GetNewBoolValue(newValue)); 
	}


}


