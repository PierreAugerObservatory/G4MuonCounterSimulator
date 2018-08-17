/**
* @file UserEventInformation.cc
* @class UserEventInformation
* @brief User-defined container of simulation information produced in the event
* @author Dr. Simone Riggi
* @date 05/04/2010
*/

#include "UserEventInformation.hh"

using namespace G4MuonCounterSimulatorUSC;

UserEventInformation::UserEventInformation()
  :hitCount(0),hitScintCount(0),hitCerenkCount(0),hitWLSCount(0),photonCount_Scint(0),photonCount_Ceren(0),photonCount_WLS(0),photonCount_lost(0),photonCount_killed(0),absorptionCount(0),boundaryAbsorptionCount(0),boundaryStepTooSmallCount(0),totE(0.),StripNo(-1),eWeightPos(0.),reconPos(0.),convPos(0.),convPosSet(false),posMax(0.),pmtsAboveThreshold(0),PrimaryParticleDirection(0.)
{

	//GENERATION VECTORS
	PrimaryParticleDirection.clear();
	PrimaryParticleDirection.resize(0);

	PrimaryParticleMomentum.clear();
	PrimaryParticleMomentum.resize(0);

	PrimaryParticleVertex.clear();
	PrimaryParticleVertex.resize(0);

	PrimaryParticleEnergy.clear();
	PrimaryParticleEnergy.resize(0);


  PrimaryParticleTheta.clear();
	PrimaryParticleTheta.resize(0);


	PrimaryParticleMass.clear();
	PrimaryParticleMass.resize(0);

	PrimaryParticleTime.clear();
	PrimaryParticleTime.resize(0);

	PrimaryParticlePDGCode.clear();
	PrimaryParticlePDGCode.resize(0);

	//SCINTILLATOR VECTORS
	scint_EmissionAngles.clear();
	scint_EmissionAngles.resize(0);

	scint_EmissionTime.clear();
	scint_EmissionTime.resize(0);

	scint_EmissionWavelength.clear();
	scint_EmissionWavelength.resize(0);

	scint_TrackLength.clear();
	scint_TrackLength.resize(0);

	scint_HorizontalTrackLength.clear();
	scint_HorizontalTrackLength.resize(0);

	scint_ProcessType.clear();
	scint_ProcessType.resize(0);

	scint_StripNumber.clear();
	scint_StripNumber.resize(0);

	scint_PlaneNumber.clear();
	scint_PlaneNumber.resize(0);

	scint_SuperPlaneNumber.clear();
	scint_SuperPlaneNumber.resize(0);

	scint_StripNumber2.clear();
	scint_StripNumber2.resize(0);

	scint_PlaneNumber2.clear();
	scint_PlaneNumber2.resize(0);

	scint_SuperPlaneNumber2.clear();
	scint_SuperPlaneNumber2.resize(0);

	scint_StripId.clear();
	scint_StripId.resize(0);

	scint_PlaneId.clear();
	scint_PlaneId.resize(0);

	scint_SuperPlaneId.clear();
	scint_SuperPlaneId.resize(0);

	scint_ScintillationCounts.clear();
	scint_ScintillationCounts.resize(0);

	scint_CerenkovCounts.clear();
	scint_CerenkovCounts.resize(0);

	scint_WLSCounts.clear();
	scint_WLSCounts.resize(0);

	scint_AbsorptionCounts.clear();
	scint_AbsorptionCounts.resize(0);

	scint_BoundaryAbsorptionCounts.clear();
	scint_BoundaryAbsorptionCounts.resize(0);

	scint_StepTooSmallBoundaryCounts.clear();
	scint_StepTooSmallBoundaryCounts.resize(0);

	//PMT VECTORS
	pmt_Time.clear();
	pmt_Time.resize(0);

	pmt_Energy.clear();
	pmt_Energy.resize(0);

	pmt_Position.clear();
	pmt_Position.resize(0);

	pmt_Wavelength.clear();
	pmt_Wavelength.resize(0);

	pmt_StripNumber.clear();
	pmt_StripNumber.resize(0);

	pmt_PlaneNumber.clear();
	pmt_PlaneNumber.resize(0);

	pmt_SuperPlaneNumber.clear();
	pmt_SuperPlaneNumber.resize(0);

	pmt_PMTNumber.clear();
	pmt_PMTNumber.resize(0);

	pmt_StripId.clear();
	pmt_StripId.resize(0);

	pmt_PlaneId.clear();
	pmt_PlaneId.resize(0);

	pmt_SuperPlaneId.clear();
	pmt_SuperPlaneId.resize(0);
	
	pmt_Id.clear();
	pmt_Id.resize(0);

	pmt_Counts.clear();
	pmt_Counts.resize(0);

	pmt_ScintillationCounts.clear();
	pmt_ScintillationCounts.resize(0);

	pmt_CerenkovCounts.clear();
	pmt_CerenkovCounts.resize(0);

	pmt_WLSCounts.clear();
	pmt_WLSCounts.resize(0);

	//plane vectors
	muonIncomingAtPlaneSurface.clear();
	muonIncomingAtPlaneSurface.assign(10,0);
	muonOutcomingAtPlaneSurface.clear();
	muonOutcomingAtPlaneSurface.assign(10,0);
	emIncomingAtPlaneSurface.clear();
	emIncomingAtPlaneSurface.assign(10,0);
	emOutcomingAtPlaneSurface.clear();
	emOutcomingAtPlaneSurface.assign(10,0);
	hadronIncomingAtPlaneSurface.clear();
	hadronIncomingAtPlaneSurface.assign(10,0);
	hadronOutcomingAtPlaneSurface.clear();
	hadronOutcomingAtPlaneSurface.assign(10,0);


}

UserEventInformation::~UserEventInformation()
{
	
}


