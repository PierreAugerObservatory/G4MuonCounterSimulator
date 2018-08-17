/**
* @file UserEventInformation.hh
* @class UserEventInformation
* @brief User-defined container of simulation information produced in the event
* @author Dr. Simone Riggi
* @date 05/04/2010
*/
#ifndef _G4MuonCounterSimulatorUSC_UserEventInformation_h
#define _G4MuonCounterSimulatorUSC_UserEventInformation_h 1


#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "vector"

namespace G4MuonCounterSimulatorUSC {

class UserEventInformation : public G4VUserEventInformation
{
public:
	/** 
	\brief Class constructor: initialize structures.
 	*/
  UserEventInformation();
	/**
	* \brief Class destructor: free allocated memory
	*/
  ~UserEventInformation();
  /**
	* \brief Print user information
	*/
  inline void Print()const{};

	//############################
	//### GENERATION INFO
	//############################
	/**
	* \brief Get list of vertexes of primary particles
	*/
	std::vector<G4ThreeVector> GetPrimaryParticlePosition(){return PrimaryParticleVertex;}
	/**
	* \brief Set list of vertexes of primary particles
	*/
	void SetPrimaryParticlePosition(std::vector<G4ThreeVector> vect) {PrimaryParticleVertex= vect;}
	/**
	* \brief Add a vertex entry to the list of primary particle vertexes
	*/
	void InsertPrimaryParticlePosition(G4ThreeVector vect) {PrimaryParticleVertex.push_back(vect);}
	/**
	* \brief Get list of directions of primary particles
	*/
	std::vector<G4ThreeVector> GetPrimaryParticleDirection(){return PrimaryParticleDirection;}	
	/**
	* \brief Set list of directions of primary particles
	*/
	void SetPrimaryParticleDirection(std::vector<G4ThreeVector> vect) {PrimaryParticleDirection= vect;}
	/**
	* \brief Add a direction entry to the list of primary particle directions
	*/
	void InsertPrimaryParticleDirection(G4ThreeVector vect) {PrimaryParticleDirection.push_back(vect);}
	/**
	* \brief Get list of energies of primary particles
	*/
	std::vector<double> GetPrimaryParticleEnergy(){return PrimaryParticleEnergy;}
	/**
	* \brief Set list of energies of primary particles
	*/
	void SetPrimaryParticleEnergy(std::vector<double> vect) {PrimaryParticleEnergy= vect;}	
	/**
	* \brief Add an energy entry to the list of primary particle energies
	*/
	void InsertPrimaryParticleEnergy(double value) {PrimaryParticleEnergy.push_back(value);}



  /**
  * \brief Get list of theta of primary particles
	*/
	std::vector<double> GetPrimaryParticleTheta(){return PrimaryParticleTheta;}
	/**
	* \brief Set list of theta of primary particles
	*/
	void SetPrimaryParticleTheta(std::vector<double> vect) {PrimaryParticleTheta= vect;}	
	/**
	* \brief Add an theta entry to the list of primary particle theta
	*/
	void InsertPrimaryParticleTheta(double value) {PrimaryParticleTheta.push_back(value);}




	/**
	* \brief Get list of times of primary particles
	*/
  std::vector<double> GetPrimaryParticleTime(){return PrimaryParticleTime;}
	/**
	* \brief Set list of times of primary particles
	*/
	void SetPrimaryParticleTime(std::vector<double> vect) {PrimaryParticleTime= vect;}
	/**
	* \brief Add a time entry to the list of primary particle times
	*/
	void InsertPrimaryParticleTime(double value) {PrimaryParticleTime.push_back(value);}
	/**
	* \brief Get list of masses of primary particles
	*/
	std::vector<double> GetPrimaryParticleMass(){return PrimaryParticleMass;}
	/**
	* \brief Set list of masses of primary particles
	*/
	void SetPrimaryParticleMass(std::vector<double> vect) {PrimaryParticleMass= vect;}
	/**
	* \brief Add a mass entry to the list of primary particle masses
	*/
	void InsertPrimaryParticleMass(double value) {PrimaryParticleMass.push_back(value);}

	/**
	* \brief Get list of PDG codes of primary particles
	*/
	std::vector<int> GetPrimaryParticlePDGCode(){return PrimaryParticlePDGCode;}
	/**
	* \brief Set list of PDG codes of primary particles
	*/
	void SetPrimaryParticlePDGCode(std::vector<int> vect) {PrimaryParticlePDGCode= vect;}
	/**
	* \brief Add a PDG code entry to the list of primary particle PDG codes
	*/
	void InsertPrimaryParticlePDGCode(int value) {PrimaryParticlePDGCode.push_back(value);}

	/**
	* \brief Get list of momenta of primary particles
	*/
  std::vector<G4ThreeVector> GetPrimaryParticleMomentum(){return PrimaryParticleMomentum;}
	/**
	* \brief Set list of momenta of primary particles
	*/
	void SetPrimaryParticleMomentum(std::vector<G4ThreeVector> vect) {PrimaryParticleMomentum= vect;}
	/**
	* \brief Add a momentum entry to the list of primary particle momenta
	*/
	void InsertPrimaryParticleMomentum(G4ThreeVector vect) {PrimaryParticleMomentum.push_back(vect);}


	//############################
	//### SCINTILLATION INFO
	//############################
	/**
	* \brief Get list of emission angles of emitted photons
	*/
	std::vector<double> GetPhotonEmissionAngle(){return scint_EmissionAngles;}
	/**
	* \brief Set list of emission angles of emitted photons
	*/
	void SetPhotonEmissionAngle(std::vector<double> vect) {scint_EmissionAngles= vect;}
	/**
	* \brief Add an emission angle entry to the list of photon emission angles
	*/
	void InsertPhotonEmissionAngle(double value) {scint_EmissionAngles.push_back(value);}
	
	/**
	* \brief Get list of emission times of emitted photons
	*/
	std::vector<double> GetPhotonEmissionTime(){return scint_EmissionTime;}
	/**
	* \brief Set list of emission times of emitted photons
	*/
	void SetPhotonEmissionTime(std::vector<double> vect) {scint_EmissionTime= vect;}
	/**
	* \brief Add an emission time entry to the list of photon emission times
	*/
	void InsertPhotonEmissionTime(double value) {scint_EmissionTime.push_back(value);}

	/**
	* \brief Get list of emission wavelength of emitted photons
	*/
	std::vector<double> GetPhotonEmissionWavelength(){return scint_EmissionWavelength;}
	/**
	* \brief Set list of emission wavelength of emitted photons
	*/
	void SetPhotonEmissionWavelength(std::vector<double> vect) {scint_EmissionWavelength= vect;}
	/**
	* \brief Add an emission wavelength entry to the list of photon emission times
	*/
	void InsertPhotonEmissionWavelength(double value) {scint_EmissionWavelength.push_back(value);}

	/**
	* \brief Get list of process type of emitted photons
	*/
  std::vector<int> GetProcessType(){return scint_ProcessType;}
	/**
	* \brief Set list of process type of emitted photons
	*/
	void SetProcessType(std::vector<int> vect) {scint_ProcessType= vect;}
	/**
	* \brief Add a process type entry to the list of photon emission process types
	*/
	void InsertProcessType(int value) {scint_ProcessType.push_back(value);}	

	/**
	* \brief Get list of strip IDs of emitted photons
	*/
	std::vector<int> GetStripNumber(){return scint_StripNumber;}
	/**
	* \brief Set list of strip IDs of emitted photons
	*/
	void SetStripNumber(std::vector<int> vect) {scint_StripNumber= vect;}
	/**
	* \brief Add a strip ID entry to the list of photon emission strip IDs
	*/
	void InsertStripNumber(int value) {scint_StripNumber.push_back(value);}	

	/**
	* \brief Get list of plane IDs of emitted photons
	*/
	std::vector<int> GetPlaneNumber(){return scint_PlaneNumber;}	
	/**
	* \brief Set list of plane IDs of emitted photons
	*/
	void SetPlaneNumber(std::vector<int> vect) {scint_PlaneNumber= vect;}
	/**
	* \brief Add a plane ID entry to the list of photon emission plane IDs
	*/
	void InsertPlaneNumber(int value) {scint_PlaneNumber.push_back(value);}	

	/**
	* \brief Get list of super plane IDs of emitted photons
	*/
	std::vector<int> GetSuperPlaneNumber(){return scint_SuperPlaneNumber;}	
	/**
	* \brief Set list of super plane IDs of emitted photons
	*/
	void SetSuperPlaneNumber(std::vector<int> vect) {scint_SuperPlaneNumber= vect;}
	/**
	* \brief Add a super plane ID entry to the list of photon emission super plane IDs
	*/
	void InsertSuperPlaneNumber(int value) {scint_SuperPlaneNumber.push_back(value);}	
	


	/**
	* \brief Get list of tracklengths of absorbed photons
	*/
	std::vector<double> GetPhotonTrackLength(){return scint_TrackLength;}
	/**
	* \brief Set list of tracklengths of absorbed photons
	*/
	void SetPhotonTrackLength(std::vector<double> vect) {scint_TrackLength= vect;}
	/**
	* \brief Add a tracklength entry to the list of absorbed photon tracklengths
	*/
	void InsertPhotonTrackLength(double value) {scint_TrackLength.push_back(value);}	

	/**
	* \brief Get list of horizontal tracklengths of absorbed photons
	*/
	std::vector<double> GetPhotonHorizontalTrackLength(){return scint_HorizontalTrackLength;}
	/**
	* \brief Set list of horizontal tracklengths of absorbed photons
	*/
	void SetPhotonHorizontalTrackLength(std::vector<double> vect) {scint_HorizontalTrackLength= vect;}
	/**
	* \brief Add a horizontal tracklength entry to the list of absorbed photon tracklengths
	*/
	void InsertPhotonHorizontalTrackLength(double value) {scint_HorizontalTrackLength.push_back(value);}

	/**
	* \brief Get list of strip IDs of absorbed photons
	*/
  std::vector<int> GetStripNumber2(){return scint_StripNumber2;}
	/**
	* \brief Set list of plane IDs of absorbed photons
	*/
	void SetStripNumber2(std::vector<int> vect) {scint_StripNumber2= vect;}
	/**
	* \brief Add a strip ID entry to the list of absorbed photon strip IDs
	*/
	void InsertStripNumber2(int value) {scint_StripNumber2.push_back(value);}	

	/**
	* \brief Get list of plane IDs of absorbed photons
	*/
	std::vector<int> GetPlaneNumber2(){return scint_PlaneNumber2;}
	/**
	* \brief Set list of plane IDs of absorbed photons
	*/
	void SetPlaneNumber2(std::vector<int> vect) {scint_PlaneNumber2= vect;}
	/**
	* \brief Add a plane ID entry to the list of absorbed photon plane IDs
	*/
	void InsertPlaneNumber2(int value) {scint_PlaneNumber2.push_back(value);}		
	/**
	* \brief Get list of super plane IDs of absorbed photons
	*/
	std::vector<int> GetSuperPlaneNumber2(){return scint_SuperPlaneNumber2;}
	/**
	* \brief Set list of super plane IDs of absorbed photons
	*/
	void SetSuperPlaneNumber2(std::vector<int> vect) {scint_SuperPlaneNumber2= vect;}
	/**
	* \brief Add a super plane ID entry to the list of absorbed photon super plane IDs
	*/
	void InsertSuperPlaneNumber2(int value) {scint_SuperPlaneNumber2.push_back(value);}	

	/**
	* \brief Get list of absorption flags of absorbed photons
	*/
	std::vector<int> GetAbsorptionFlag(){return scint_AbsorptionFlag;}
	/**
	* \brief Set list of absorption flags of absorbed photons
	*/
	void SetAbsorptionFlag(std::vector<int> vect) {scint_AbsorptionFlag= vect;}
	/**
	* \brief Add an absorption entry to the list of absorbed photon flags
	*/
	void InsertAbsorptionFlag(int value) {scint_AbsorptionFlag.push_back(value);}		

	
	/**
	* \brief Get strip IDs of all hits
	*/
	std::vector<int> GetStripId(){return scint_StripId;}
	/**
	* \brief Set strip IDs of all hits
	*/
	void SetStripId(std::vector<int> vect) {scint_StripId= vect;}
	/**
	* \brief Add a strip ID entry to the list of all hit strip IDs 
	*/
	void InsertStripId(int value) {scint_StripId.push_back(value);}	
	/**
	* \brief Get plane IDs of all hits
	*/
	std::vector<int> GetPlaneId(){return scint_PlaneId;}
	/**
	* \brief Get plane IDs of all hits
	*/
	void SetPlaneId(std::vector<int> vect) {scint_PlaneId= vect;}
	/**
	* \brief Add a plane ID entry to the list of all hit plane IDs
	*/
	void InsertPlaneId(int value) {scint_PlaneId.push_back(value);}	
	/**
	* \brief Get super plane IDs of all hits
	*/
	std::vector<int> GetSuperPlaneId(){return scint_SuperPlaneId;}
	/**
	* \brief Get super plane IDs of all hits
	*/
	void SetSuperPlaneId(std::vector<int> vect) {scint_SuperPlaneId= vect;}
	/**
	* \brief Add a super plane ID entry to the list of all hit super plane IDs
	*/
	void InsertSuperPlaneId(int value) {scint_SuperPlaneId.push_back(value);}
	
	
	/**
	* \brief Get scintillation photon counts of all hit strips
	*/
	std::vector<int> GetScintillationCounts(){return scint_ScintillationCounts;}
	/**
	* \brief Set scintillation photon counts of all hit strips
	*/
	void SetScintillationCounts(std::vector<int> vect) {scint_ScintillationCounts= vect;}
	/**
	* \brief Add a scintillation photon count entry to the list of all hit strips
	*/
	void InsertScintillationCounts(int value) {scint_ScintillationCounts.push_back(value);}	

	/**
	* \brief Get cerenkov photon counts of all hit strips
	*/
	std::vector<int> GetCerenkovCounts(){return scint_CerenkovCounts;}
	/**
	* \brief Set cerenkov photon counts of all hit strips
	*/
	void SetCerenkovCounts(std::vector<int> vect) {scint_CerenkovCounts= vect;}
	/**
	* \brief Add a cerenkov photon count entry to the list of all hit strips
	*/
	void InsertCerenkovCounts(int value) {scint_CerenkovCounts.push_back(value);}	

	/**
	* \brief Get WLS photon counts of all hit strips
	*/
	std::vector<int> GetWLSCounts(){return scint_WLSCounts;}
	/**
	* \brief Set WLS photon counts of all hit strips
	*/
	void SetWLSCounts(std::vector<int> vect) {scint_WLSCounts= vect;}
	/**
	* \brief Add a WLS photon count entry to the list of all hit strips
	*/
	void InsertWLSCounts(int value) {scint_WLSCounts.push_back(value);}	

	/**
	* \brief Get absorption photon counts of all hit strips
	*/
	std::vector<int> GetAbsorptionCounts(){return scint_AbsorptionCounts;}
	/**
	* \brief Set absorption photon counts of all hit strips
	*/
	void SetAbsorptionCounts(std::vector<int> vect) {scint_AbsorptionCounts= vect;}
	/**
	* \brief Add an absorption photon count entry to the list of all hit strips
	*/
	void InsertAbsorptionCounts(int value) {scint_AbsorptionCounts.push_back(value);}	
	
	/**
	* \brief Get boundary absorption photon counts of all hit strips
	*/
	std::vector<int> GetBoundaryAbsorptionCounts(){return scint_BoundaryAbsorptionCounts;}
	/**
	* \brief Set boundary absorption photon counts of all hit strips
	*/
	void SetBoundaryAbsorptionCounts(std::vector<int> vect) {scint_BoundaryAbsorptionCounts= vect;}
	/**
	* \brief Add a boundary absorption photon count entry to the list of all hit strips
	*/
	void InsertBoundaryAbsorptionCounts(int value) {scint_BoundaryAbsorptionCounts.push_back(value);}	
	
	/**
	* \brief Get boundary absorption step too small photon counts of all hit strips
	*/
	std::vector<int> GetStepTooSmallBoundaryCounts(){return scint_StepTooSmallBoundaryCounts;}
	/**
	* \brief Set boundary absorption step too small photon counts of all hit strips
	*/
	void SetStepTooSmallBoundaryCounts(std::vector<int> vect) {scint_StepTooSmallBoundaryCounts= vect;}
	/**
	* \brief Add a boundary absorption step too small photon count entry to the list of all hit strips
	*/
	void InsertStepTooSmallBoundaryCounts(int value) {scint_StepTooSmallBoundaryCounts.push_back(value);}	
	
	//cumulative
	/**
	* \brief Get total number of photons produced in the simulation
	*/
  int GetPhotonCount(){return photonCount_Scint+photonCount_Ceren+photonCount_WLS;}
	/**
	* \brief Get total number of scintillation photons produced in the simulation
	*/
	int GetPhotonCount_Scint()const {return photonCount_Scint;}
	/**
	* \brief Increment the total number of scintillation photons produced in the simulation (+1)
	*/
  void IncPhotonCount_Scint(){photonCount_Scint++;}
  /**
	* \brief Get total number of cerenkov photons produced in the simulation
	*/
	int GetPhotonCount_Ceren()const {return photonCount_Ceren;}
	/**
	* \brief Increment the total number of cerenkov photons produced in the simulation (+1)
	*/
	void IncPhotonCount_Ceren(){photonCount_Ceren++;}
	/**
	* \brief Get total number of WLS photons produced in the simulation
	*/
	int GetPhotonCount_WLS()const {return photonCount_WLS;}
	/**
	* \brief Increment the total number of WLS photons produced in the simulation (+1)
	*/
	void IncPhotonCount_WLS(){photonCount_WLS++;}
	/**
	* \brief Get total number of lost photons in the simulation
	*/
  int GetPhotonCount_lost()const {return photonCount_lost;}
	/**
	* \brief Increment the total number of lost photons in the simulation (+1)
	*/
	void IncPhotonCount_lost(){photonCount_lost++;}
	/**
	* \brief Get total number of killed photons in the simulation
	*/
  int GetPhotonCount_killed()const {return photonCount_killed;}
	/**
	* \brief Increment the total number of killed photons in the simulation (+1)
	*/
	void IncPhotonCount_killed(){photonCount_killed++;}
	/**
	* \brief Get the total energy deposit in the simulation
	*/
	double GetEDep()const {return totE;}
	/**
	* \brief Increment the total energy deposit in the simulation with current energy loss
	*/
  void IncEDep(double dep){totE+=dep;}
	/**
	* \brief Get total number of absorbed photons in the simulation
	*/
	int GetAbsorptionCount()const {return absorptionCount;}
	/**
	* \brief Increment the total number of absorbed photons in the simulation (+1)
	*/
  void IncAbsorption(){absorptionCount++;}
	/**
	* \brief Get total number of absorbed photons at boundaries in the simulation
	*/
	int GetBoundaryAbsorptionCount() const {return boundaryAbsorptionCount;}
	/**
	* \brief Increment the total number of absorbed photons at boundaries in the simulation (+1)
	*/
  void IncBoundaryAbsorption(){boundaryAbsorptionCount++;}
	/**
	* \brief Get total number of absorbed photons at steps too small in the simulation
	*/
	int GetBoundaryStepTooSmallCount() const {return boundaryStepTooSmallCount;}
	/**
	* \brief Increment the total number of absorbed photons at steps too small in the simulation (+1)
	*/
	void IncBoundaryStepTooSmall(){boundaryStepTooSmallCount++;}

	/**
	* \brief Set energy-weighted hit position
	*/
  void SetEWeightPos(const G4ThreeVector& p){eWeightPos=p;}
	/**
	* \brief Get energy-weighted hit position
	*/
	G4ThreeVector GetEWeightPos(){return eWeightPos;}

	/**
	* \brief ### TO BE DOCUMENTED ###
	*/
	void SetReconPos(const G4ThreeVector& p){reconPos=p;}
	/**
	* \brief ### TO BE DOCUMENTED ###
	*/
  G4ThreeVector GetReconPos(){return reconPos;}
	/**
	* \brief ### TO BE DOCUMENTED ###
	*/
	void SetConvPos(const G4ThreeVector& p){convPos=p;convPosSet=true;}
	/**
	* \brief ### TO BE DOCUMENTED ###
	*/
  G4ThreeVector GetConvPos(){return convPos;}
	/**
	* \brief ### TO BE DOCUMENTED ###
	*/
	void SetPosMax(const G4ThreeVector& p,double edep){posMax=p;edepMax=edep;}
	/**
	* \brief ### TO BE DOCUMENTED ###
	*/
  G4ThreeVector GetPosMax(){return posMax;}
	/**
	* \brief ### TO BE DOCUMENTED ###
	*/
  double GetEDepMax(){return edepMax;}
	/**
	* \brief ### TO BE DOCUMENTED ###
	*/
  double IsConvPosSet(){return convPosSet;}
  
	
	//############################
	//### PMT INFO
	//############################
	/**
	* \brief Get arrival times of photons at the PMT
	*/
  std::vector<double> GetPMTPhotonTime() {return pmt_Time;}
	/**
	* \brief Set arrival times of photons at the PMT
	*/
	void SetPMTPhotonTime(std::vector<double> vect) {pmt_Time= vect;}
	/**
	* \brief Add an arrival time entry to the list of arrival times of photons at the PMT
	*/	
  void InsertPMTPhotonTime(double value) {pmt_Time.push_back(value);}	
	/**
	* \brief Get energy of photons at the PMT
	*/
  std::vector<double> GetPMTPhotonEnergy() {return pmt_Energy;}
	/**
	* \brief Set energy of photons at the PMT
	*/
	void SetPMTPhotonEnergy(std::vector<double> vect) {pmt_Energy= vect;}	
	/**
	* \brief Add an energy entry to the list of energies of photons at the PMT
	*/
  void InsertPMTPhotonEnergy(double value) {pmt_Energy.push_back(value);}	
	/**
	* \brief Get wavelength of photons at the PMT
	*/
  std::vector<double> GetPMTPhotonWavelength() {return pmt_Wavelength;}
	/**
	* \brief Set wavelength of photons at the PMT
	*/
	void SetPMTPhotonWavelength(std::vector<double> vect) {pmt_Wavelength= vect;}	
	/**
	* \brief Add a wavelength entry to the list of wavelengths of photons at the PMT
	*/
  void InsertPMTPhotonWavelength(double value) {pmt_Wavelength.push_back(value);}	
	/**
	* \brief Get arrival position of photons at the PMT
	*/
  std::vector<G4ThreeVector> GetPMTPhotonPosition() {return pmt_Position;}
	/**
	* \brief Set arrival position of photons at the PMT
	*/
	void SetPMTPhotonPosition(std::vector<G4ThreeVector> vect) {pmt_Position= vect;}
  /**
	* \brief Add an arrival position entry to the list of arrival position of photons at the PMT
	*/	
  void InsertPMTPhotonPosition(G4ThreeVector value) {pmt_Position.push_back(value);}	

	/**
	* \brief Get strip number of hit PMT
 	*/
  std::vector<int> GetPMTStripNumber() {return pmt_StripNumber;}
	/**
	* \brief Set strip number of hit PMT
 	*/
	void SetPMTStripNumber(std::vector<int> vect) {pmt_StripNumber= vect;}	
	/**
	* \brief Add strip number of hit PMT
 	*/
  void InsertPMTStripNumber(int value) {pmt_StripNumber.push_back(value);}
	/**
	* \brief Get plane number of hit PMT
 	*/
	std::vector<int> GetPMTPlaneNumber() {return pmt_PlaneNumber;}
	/**
	* \brief Set plane number of hit PMT
 	*/
	void SetPMTPlaneNumber(std::vector<int> vect) {pmt_PlaneNumber= vect;}
	/**
	* \brief Add plane number of hit PMT
 	*/	
  void InsertPMTPlaneNumber(int value) {pmt_PlaneNumber.push_back(value);}
	
	/**
	* \brief Get super plane number of hit PMT
 	*/
	std::vector<int> GetPMTSuperPlaneNumber() {return pmt_SuperPlaneNumber;}
	/**
	* \brief Set super plane number of hit PMT
 	*/
	void SetPMTSuperPlaneNumber(std::vector<int> vect) {pmt_SuperPlaneNumber= vect;}
	/**
	* \brief Add super plane number of hit PMT
 	*/	
  void InsertPMTSuperPlaneNumber(int value) {pmt_SuperPlaneNumber.push_back(value);}


	/**
	* \brief Get PMT number of hit PMT
 	*/
	std::vector<int> GetPMTNumber() {return pmt_PMTNumber;}
	/**
	* \brief Set PMT number of hit PMT
 	*/
	void SetPMTNumber(std::vector<int> vect) {pmt_PMTNumber= vect;}	
	/**
	* \brief Add PMT number of hit PMT
 	*/
  void InsertPMTNumber(int value) {pmt_PMTNumber.push_back(value);}	
	
	/**
	* \brief Get strip Id of hit PMT
 	*/
  std::vector<int> GetPMTStripId() {return pmt_StripId;}
	/**
	* \brief Set strip Id of hit PMT
 	*/
	void SetPMTStripId(std::vector<int> vect) {pmt_StripId= vect;}		
	/**
	* \brief Add strip Id of hit PMT
 	*/
  void InsertPMTStripId(int value) {pmt_StripId.push_back(value);}
	/**
	* \brief Get plane Id of hit PMT
 	*/
	std::vector<int> GetPMTPlaneId() {return pmt_PlaneId;}
	/**
	* \brief Set plane Id of hit PMT
 	*/
	void SetPMTPlaneId(std::vector<int> vect) {pmt_PlaneId= vect;}
	/**
	* \brief Add plane Id of hit PMT
 	*/	
  void InsertPMTPlaneId(int value) {pmt_PlaneId.push_back(value);}	
	
	/**
	* \brief Get super plane Id of hit PMT
 	*/
	std::vector<int> GetPMTSuperPlaneId() {return pmt_SuperPlaneId;}
	/**
	* \brief Set super plane Id of hit PMT
 	*/
	void SetPMTSuperPlaneId(std::vector<int> vect) {pmt_SuperPlaneId= vect;}
	/**
	* \brief Add super plane Id of hit PMT
 	*/	
  void InsertPMTSuperPlaneId(int value) {pmt_SuperPlaneId.push_back(value);}	

	/**
	* \brief Get PMT Id of hit PMT
 	*/
	std::vector<int> GetPMTId() {return pmt_Id;}
	/**
	* \brief Set PMT Id of hit PMT
 	*/
	void SetPMTId(std::vector<int> vect) {pmt_Id= vect;}
	/**
	* \brief Add PMT Id of hit PMT
 	*/	
  void InsertPMTId(int value) {pmt_Id.push_back(value);}
	/**
	* \brief Get photon counts of hit PMT
 	*/
	std::vector<int> GetPMTCounts() {return pmt_Counts;}
	/**
	* \brief Set photon counts of hit PMT
 	*/
	void SetPMTCounts(std::vector<int> vect) {pmt_Counts= vect;}	
	/**
	* \brief Add photon counts of hit PMT
 	*/
  void InsertPMTCounts(int value) {pmt_Counts.push_back(value);}
	/**
	* \brief Get scintillation photon counts of hit PMT
 	*/
	std::vector<int> GetPMTScintillationCounts() {return pmt_ScintillationCounts;}
	/**
	* \brief Set scintillation photon counts of hit PMT
 	*/
	void SetPMTScintillationCounts(std::vector<int> vect) {pmt_ScintillationCounts= vect;}	
	/**
	* \brief Add scintillation photon counts of hit PMT
 	*/
  void InsertPMTScintillationCounts(int value) {pmt_ScintillationCounts.push_back(value);}
	/**
	* \brief Get cerenkov photon counts of hit PMT
 	*/
  std::vector<int> GetPMTCerenkovCounts() {return pmt_CerenkovCounts;}
	/**
	* \brief Set cerenkov photon counts of hit PMT
 	*/
	void SetPMTCerenkovCounts(std::vector<int> vect) {pmt_CerenkovCounts= vect;}	
	/**
	* \brief Add cerenkov photon counts of hit PMT
 	*/
  void InsertPMTCerenkovCounts(int value) {pmt_CerenkovCounts.push_back(value);}
	/**
	* \brief Get WLS photon counts of hit PMT
 	*/
	std::vector<int> GetPMTWLSCounts() {return pmt_WLSCounts;}
	/**
	* \brief Set WLS photon counts of hit PMT
 	*/
	void SetPMTWLSCounts(std::vector<int> vect) {pmt_WLSCounts= vect;}
	/**
	* \brief Add WLS photon counts of hit PMT
 	*/	
  void InsertPMTWLSCounts(int value) {pmt_WLSCounts.push_back(value);}

	//cumulative
	/**
	* \brief Get total photon counts
 	*/
	int GetHitCount()const {return hitCount;}
	/**
	* \brief Add total photon counts
 	*/
	void IncHitCount(int i=1){hitCount+=i;}
	/**
	* \brief Get total scintillation photon counts
 	*/
	int GetHitScintCount()const {return hitScintCount;}
	/**
	* \brief Increment the total scintillation photon counts (+1)
 	*/
	void IncHitScintCount(int i=1){hitScintCount+=i;}
	/**
	* \brief Get total cerenkov photon counts
 	*/
	int GetHitCerenkCount()const {return hitCerenkCount;}
	/**
	* \brief Increment the total cerenkov photon counts (+1)
 	*/
	void IncHitCerenkCount(int i=1){hitCerenkCount+=i;}
	/**
	* \brief Get total WLS photon counts
 	*/
	int GetHitWLSCount()const {return hitWLSCount;}
	/**
	* \brief Increment the total WLS photon counts (+1)
 	*/
	void IncHitWLSCount(int i=1){hitWLSCount+=i;}
	/**
	* \brief Increment the number of PMT above threshold (+1)
 	*/
	void IncPMTSAboveThreshold(){pmtsAboveThreshold++;}
	/**
	* \brief Get number of PMT above threshold
 	*/
  int GetPMTSAboveThreshold(){return pmtsAboveThreshold;}


	/**
	* \brief Get nMuons incoming at the plane surface
 	*/
	std::vector<double> GetInMuonsAtPlaneSurface(){return muonIncomingAtPlaneSurface;}
	/**
	* \brief Add nMuons incoming at the plane surface
 	*/	
  void IncInMuonsAtPlaneSurface(int planeId,double value) {muonIncomingAtPlaneSurface[planeId]+= value;}
	/**
	* \brief Get nMuons outcoming at the plane surface
 	*/
	std::vector<double> GetOutMuonsAtPlaneSurface(){return muonOutcomingAtPlaneSurface;}
	/**
	* \brief Add nMuons outcoming at the plane surface
 	*/	
  void IncOutMuonsAtPlaneSurface(int planeId,double value) {muonOutcomingAtPlaneSurface[planeId]+= value;}
	/**
	* \brief Get nEm incoming at the plane surface
 	*/
	std::vector<double> GetInEmAtPlaneSurface(){return emIncomingAtPlaneSurface;}
	/**
	* \brief Add nEm incoming at the plane surface
 	*/	
  void IncInEmAtPlaneSurface(int planeId,double value) {emIncomingAtPlaneSurface[planeId]+= value;}	
	/**
	* \brief Get nEm outcoming at the plane surface
 	*/
	std::vector<double> GetOutEmAtPlaneSurface(){return emOutcomingAtPlaneSurface;}
	/**
	* \brief Add nEm outcoming at the plane surface
 	*/	
  void IncOutEmAtPlaneSurface(int planeId,double value) {emOutcomingAtPlaneSurface[planeId]+= value;}
	/**
	* \brief Get nHadrons incoming at the plane surface
 	*/
	std::vector<double> GetInHadronsAtPlaneSurface(){return hadronIncomingAtPlaneSurface;}
	/**
	* \brief Add nHadrons incoming at the plane surface
 	*/	
  void IncInHadronsAtPlaneSurface(int planeId,double value) {hadronIncomingAtPlaneSurface[planeId]+= value;}
	/**
	* \brief Get nHadrons outcoming at the plane surface
 	*/
	std::vector<double> GetOutHadronsAtPlaneSurface(){return hadronOutcomingAtPlaneSurface;}
	/**
	* \brief Add nHadrons outcoming at the plane surface
 	*/	
  void IncOutHadronsAtPlaneSurface(int planeId,double value) {hadronOutcomingAtPlaneSurface[planeId]+= value;}

private:

	//## GENERATED INFO
	//G4ThreeVector PrimaryParticleDirection;
	/**
	* \brief Vector of primary particle directions
 	*/
	std::vector<G4ThreeVector> PrimaryParticleDirection;
	/**
	* \brief Vector of primary particle vertexes
 	*/
	std::vector<G4ThreeVector> PrimaryParticleVertex;
	/**
	* \brief Vector of primary particle energies
 	*/
	std::vector<double> PrimaryParticleEnergy;

  /**
	* \brief Vector of primary particle theta
 	*/
	std::vector<double> PrimaryParticleTheta;




	/**
	* \brief Vector of primary particle masses
 	*/
	std::vector<double> PrimaryParticleMass;
	/**
	* \brief Vector of primary particle generation times
 	*/
	std::vector<double> PrimaryParticleTime;
	/**
	* \brief Vector of primary particle PDG codes
 	*/
	std::vector<int> PrimaryParticlePDGCode;
	/**
	* \brief Vector of primary particle momenta
 	*/
	std::vector<G4ThreeVector> PrimaryParticleMomentum;

	//## SCINTILLATOR vars
	//vectors looping over all photons circulating in the strips
	/**
	* \brief Vector of photon emission angles in the strip (loop over all photons circulating in the strips)
 	*/
	std::vector<double> scint_EmissionAngles;
	/**
	* \brief Vector of photon emission times in the strip (loop over all photons circulating in the strips)
 	*/
	std::vector<double> scint_EmissionTime;
	/**
	* \brief Vector of photon emission wavelength in the strip (loop over all photons circulating in the strips)
 	*/
	std::vector<double> scint_EmissionWavelength;
	/**
	* \brief Vector of photon generation process in the strip (loop over all photons circulating in the strips)
 	*/
	std::vector<int> scint_ProcessType;
	/**
	* \brief Vector of strip IDs of emitted photons (loop over all photons circulating in the strips)
 	*/
  std::vector<int> scint_StripNumber;
	/**
	* \brief Vector of plane IDs of emitted photons (loop over all photons circulating in the strips)
 	*/
	std::vector<int> scint_PlaneNumber;
	/**
	* \brief Vector of super plane IDs of emitted photons (loop over all photons circulating in the strips)
 	*/
	std::vector<int> scint_SuperPlaneNumber;

  //vectors looping over all photons absorbed in the strips
	/**
	* \brief Vector of tracklengths of absorbed photons (loop over all photons absorbed in the strips)
 	*/
  std::vector<double> scint_TrackLength;
	/**
	* \brief Vector of horizontal tracklengths of absorbed photons (loop over all photons absorbed in the strips)
 	*/
	std::vector<double> scint_HorizontalTrackLength;
	/**
	* \brief Vector of strip IDs of absorbed photons (loop over all photons absorbed in the strips)
 	*/
	std::vector<int> scint_StripNumber2;
	/**
	* \brief Vector of plane IDs of absorbed photons (loop over all photons absorbed in the strips)
 	*/
	std::vector<int> scint_PlaneNumber2;
	/**
	* \brief Vector of super plane IDs of absorbed photons (loop over all photons absorbed in the strips)
 	*/
	std::vector<int> scint_SuperPlaneNumber2;
	/**
	* \brief Vector of absorption flag of absorbed photons (loop over all photons absorbed in the strips)
 	*/
	std::vector<int> scint_AbsorptionFlag;

	//vectors looping over all strips
	/**
	* \brief Vector of strip IDs (loop over all strips)
 	*/
	std::vector<int> scint_StripId;
	/**
	* \brief Vector of plane IDs (loop over all strips)
 	*/
	std::vector<int> scint_PlaneId;
	/**
	* \brief Vector of super plane IDs (loop over all strips)
 	*/
	std::vector<int> scint_SuperPlaneId;
	/**
	* \brief Vector of scintillation counts (loop over all strips)
 	*/
	std::vector<int> scint_ScintillationCounts;
	/**
	* \brief Vector of cerenkov counts (loop over all strips)
 	*/
	std::vector<int> scint_CerenkovCounts;	
	/**
	* \brief Vector of WLS counts (loop over all strips)
 	*/
	std::vector<int> scint_WLSCounts;
	/**
	* \brief Vector of absorption counts (loop over all strips)
 	*/
	std::vector<int> scint_AbsorptionCounts;
	/**
	* \brief Vector of boundary absorption counts (loop over all strips)
 	*/
	std::vector<int> scint_BoundaryAbsorptionCounts;
	/**
	* \brief Vector of boundary step-too-small counts (loop over all strips)
 	*/
	std::vector<int> scint_StepTooSmallBoundaryCounts;

	//cumulative
	/**
	* \brief Total number of scintillation counts
 	*/
	int photonCount_Scint;
	/**
	* \brief Total number of cerenkov counts
 	*/
  int photonCount_Ceren;
	/**
	* \brief Total number of WLS counts
 	*/
	int photonCount_WLS;
	/**
	* \brief Total number of lost photons
 	*/
  int photonCount_lost;
	/**
	* \brief Total number of killed photons
 	*/
	int photonCount_killed;
	/**
	* \brief Total number of absorbed photons
 	*/	
  int absorptionCount;
	/**
	* \brief Total number of photons absorbed at boundaries
 	*/
  int boundaryAbsorptionCount;
	/**
	* \brief Total number of photons absorbed at boundary with a step-too-small
 	*/
	int boundaryStepTooSmallCount;
	/**
	* \brief Total energy deposit
 	*/
	double totE; 
	/**
	* \brief Strip ID ### TO BE REMOVED ###
 	*/
  int StripNo; 

	//These only have meaning if totE > 0
  //If totE = 0 then these wont be set by EndOfEventAction
	/**
	* \brief Energy-weighted hit position
 	*/
  G4ThreeVector eWeightPos;
	/**
	* \brief ### TO BE DOCUMENTED/REMOVED ###
 	*/
  G4ThreeVector reconPos; //Also relies on hitCount>0
	/**
	* \brief ### TO BE DOCUMENTED/REMOVED ###
 	*/
  G4ThreeVector convPos;//true (initial) converstion position
	/**
	* \brief ### TO BE DOCUMENTED/REMOVED ###
 	*/
  bool convPosSet;
	/**
	* \brief ### TO BE DOCUMENTED/REMOVED ###
 	*/
  G4ThreeVector posMax;
	/**
	* \brief ### TO BE DOCUMENTED/REMOVED ###
 	*/
  double edepMax;
	

	//##PMT vars
	//vectors looping over all photons hitting the PMTs
	/**
	* \brief Vector of photon arrival times at the PMT (over all photons hitting the PMTs)
 	*/
  std::vector<double> pmt_Time; 
	/**
	* \brief Vector of photon energies at the PMT (over all photons hitting the PMTs)
 	*/
	std::vector<double> pmt_Energy;
	/**
	* \brief Vector of photon wavelengths at the PMT (over all photons hitting the PMTs)
 	*/
  std::vector<double> pmt_Wavelength;	
	/**
	* \brief Vector of photon arrival position at the PMT (over all photons hitting the PMTs)
 	*/
  std::vector<G4ThreeVector> pmt_Position; 
	/**
	* \brief Vector of strip ID of photons at the PMT (over all photons hitting the PMTs)
 	*/
  std::vector<int> pmt_StripNumber;
	/**
	* \brief Vector of plane ID of photons at the PMT (over all photons hitting the PMTs)
 	*/
	std::vector<int> pmt_PlaneNumber;
	/**
	* \brief Vector of super plane ID of photons at the PMT (over all photons hitting the PMTs)
 	*/
	std::vector<int> pmt_SuperPlaneNumber;
	/**
	* \brief Vector of PMT ID of photons at the PMT (over all photons hitting the PMTs)
 	*/
	std::vector<int> pmt_PMTNumber;

	//vectors looping over strips
	/**
	* \brief Vector of strip ID of photons at the PMT (over all hit strips)
 	*/
	std::vector<int> pmt_StripId;
	/**
	* \brief Vector of plane ID of photons at the PMT (over all hit strips)
 	*/
	std::vector<int> pmt_PlaneId;
	/**
	* \brief Vector of super plane ID of photons at the PMT (over all hit strips)
 	*/
	std::vector<int> pmt_SuperPlaneId;
	/**
	* \brief Vector of PMT ID of photons at the PMT (over all hit strips)
 	*/
	std::vector<int> pmt_Id;
	/**
	* \brief Vector of PMT counts (over all hit strips)
 	*/
	std::vector<int> pmt_Counts;
	/**
	* \brief Vector of PMT scintillation counts (over all hit strips)
 	*/
	std::vector<int> pmt_ScintillationCounts;
	/**
	* \brief Vector of PMT cerenkov counts (over all hit strips)
 	*/
  std::vector<int> pmt_CerenkovCounts;
	/**
	* \brief Vector of PMT WLS counts (over all hit strips)
 	*/
	std::vector<int> pmt_WLSCounts;

	//cumulative info
	/**
	* \brief Total PMT counts
 	*/
	int hitCount;
	/**
	* \brief Total PMT scintillation ounts
 	*/
	int hitScintCount;
	/**
	* \brief Total PMT cerenkov counts
 	*/
	int hitCerenkCount;
	/**
	* \brief Total PMT WLS counts
 	*/
	int hitWLSCount;
	/**
	* \brief Number of PMTs above threshold
 	*/
	int pmtsAboveThreshold;


	/**
	* \brief Number of muons incoming the detector plane surface
 	*/
	std::vector<double> muonIncomingAtPlaneSurface;
	/**
	* \brief Number of muons outcoming the detector plane surface
 	*/
	std::vector<double> muonOutcomingAtPlaneSurface;
	/**
	* \brief Number of em particles incoming the detector plane surface
 	*/
	std::vector<double> emIncomingAtPlaneSurface;
	/**
	* \brief Number of em particles outcoming the detector plane surface
 	*/
	std::vector<double> emOutcomingAtPlaneSurface;
	/**
	* \brief Number of proton/neutron particles incoming the detector plane surface
 	*/
	std::vector<double> hadronIncomingAtPlaneSurface;
	/**
	* \brief Number of proton/neutron particles outcoming the detector plane surface
 	*/
	std::vector<double> hadronOutcomingAtPlaneSurface;

	

};

}//close namespace

#endif





