/**
   \file
   implementation of class ShowerRegenerator, adapted from CachedShowerRegenerator
   \version $Id$
*/

static const char CVSId[] =
  "$Id$";


#include "ShowerRegenerator.h"
#include "LogGaussSmearing.h"

#include <utl/config.h>

#include <cmath>
#include <sstream>
#include <vector>
#include <limits>

#include <fwk/CentralConfig.h>
#include <fwk/CoordinateSystemRegistry.h>
#include <fwk/RandomEngineRegistry.h>

#include <evt/Event.h>
#include <evt/ShowerSimData.h>

#include <sevt/SEvent.h>
#include <sevt/Station.h>
#include <sevt/StationSimData.h>

//#include "MyDetector.h"
//#include <MDetector.h>
//#include <MStation.h>

#include <det/Detector.h>
#include <sdet/SDetector.h>
#include <sdet/Station.h>

#include <utl/ErrorLogger.h>
#include <utl/GeometryUtilities.h>
#include <utl/MathConstants.h>
#include <utl/Particle.h>
#include <utl/PhysicalConstants.h>
#include <utl/Point.h>
#include <utl/UTMPoint.h>
#include <utl/RandomEngine.h>
#include <utl/TimeStamp.h>
#include <utl/Math.h>
#include <utl/Reader.h>
#include <utl/TabularStream.h>
#include <utl/TabulatedFunction.h>

#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Random/RandPoisson.h>
#include <CLHEP/Random/RandGauss.h>

#include <TFile.h>
#include <TTree.h>

using namespace G4MuonCounterShowerRegenerator;
using namespace fwk;
using namespace evt;
using namespace det;
//using namespace mdet;
using namespace sevt;
using namespace sevt;
using namespace utl;
using namespace std;

using CLHEP::RandFlat;
using CLHEP::RandPoisson;

// defs for particle type switch
#define PHOTONS utl::Particle::ePhoton
#define ELECTRONS utl::Particle::eElectron: \
             case utl::Particle::ePositron
#define MUONS utl::Particle::eMuon: \
         case utl::Particle::eAntiMuon
#define BARYONS utl::Particle::eProton: \
           case utl::Particle::eAntiProton: \
           case utl::Particle::eNeutron: \
           case utl::Particle::eAntiNeutron: \
           case utl::Particle::eLambda: \
           case utl::Particle::eAntiLambda
#define MESONS utl::Particle::ePiZero: \
          case utl::Particle::ePiPlus: \
          case utl::Particle::ePiMinus: \
          case utl::Particle::eEta: \
          case utl::Particle::eKaon0L: \
          case utl::Particle::eKaon0S: \
          case utl::Particle::eKaonPlus: \
          case utl::Particle::eKaonMinus

#ifndef OFFLINE_MINOR_VERSION
#define OFFLINE_MINOR_VERSION 0 //whatever value smaller than release 6 Asterix-Obelix
#endif


namespace G4MuonCounterShowerRegenerator {

  template<typename Map, typename T>
  inline
  void
  InsertValue(Map& map, const int sId, const T& value)
  {
    const typename Map::iterator sIt = map.find(sId);
    if (sIt != map.end())
      sIt->second(value);
    else
      map.insert(make_pair(sId, value));
  }

  inline
  double
  Round(const double div, const double val)
  {
    return round(div*val)/div;
  }

  /// Calculate time of arrival of the plan front at point x
  /**
      t_fp = t_c - (x - x_c)_z / c
      x_c = core position
      t_c = core time
  */
  inline
  double
  PlaneFrontTime(const CoordinateSystemPtr showerCS, const Point& corePosition,
                 const Point& position)
  {
    return (corePosition.GetZ(showerCS) - position.GetZ(showerCS)) / kSpeedOfLight;
  }

} // namespace ShowerRegenerator


ShowerRegenerator::ShowerRegenerator() :
  fMaxParticles(numeric_limits<unsigned int>::max()),
  fInnerRadiusCut(0.),
  fOuterRadiusCut(1.e6*km),
  fElectronEnergyCut(0),
  fMuonEnergyCut(0),
  fPhotonEnergyCut(0),
  fHadronEnergyCut(0),
  fMesonEnergyCut(0),
	fIsElectronSkipped(0),
  fIsMuonSkipped(0),
  fIsPhotonSkipped(0),
  fIsHadronSkipped(0),
  fIsMesonSkipped(0),
  fDeltaROverR(0),
  fDeltaPhi(0),
  fHorizontalParticleCut(0),
  fUseStationPositionMatrix(true),
  fPhiGranularity(0),
  fRGranularity(0),
	fSamplingDetAreaSizeX(0),
	fSamplingDetAreaSizeY(0),
  fMuToVem(0),
  fEToVem(0),
  fUseWeightLimiting(false),
  fWeightThreshold(0),
  fAccumulatedWeightLimit(0),
	fSaveInjParticlesInfo(0),
	fOutputFile(0),
	fInjParticleInfo(0),
  fRandomEngine(0),
  fTimeSlotSize(1.*ns),
  fNOutOfRangeWeight(0),
  fWeightLimit(0)
{
}


VModule::ResultFlag 
ShowerRegenerator::Init()
{

  unityWeightCtr = 0; // temp
  largeWeightCtr = 0; // temp

  Branch topB= CentralConfig::GetInstance()->GetTopBranch("ShowerRegenerator");

  //AttributeMap useAtt;
	map<string, string> useAtt;
  useAtt["use"] = string("yes");

  Branch limitB = topB.GetChild("LimitParticlesPerCycle", useAtt);
  if (limitB) limitB.GetData(fMaxParticles);

  Branch distanceCutB = topB.GetChild("DistanceCuts");

  Branch inB = distanceCutB.GetChild("InnerRadiusCut", useAtt);
  if (inB) inB.GetData(fInnerRadiusCut);

  Branch outB = distanceCutB.GetChild("OuterRadiusCut", useAtt);
  if (outB) outB.GetData(fOuterRadiusCut);

  Branch energyCutB = topB.GetChild("EnergyCuts");
  energyCutB.GetChild("ElectronEnergyCut").GetData(fElectronEnergyCut);
  energyCutB.GetChild("MuonEnergyCut").GetData(fMuonEnergyCut);
  energyCutB.GetChild("PhotonEnergyCut").GetData(fPhotonEnergyCut);
  energyCutB.GetChild("HadronEnergyCut").GetData(fHadronEnergyCut);
  energyCutB.GetChild("MesonEnergyCut").GetData(fMesonEnergyCut);

	Branch particleSkipB = topB.GetChild("ParticleSkip");
  particleSkipB.GetChild("SkipElectrons").GetData(fIsElectronSkipped);
  particleSkipB.GetChild("SkipMuons").GetData(fIsMuonSkipped);
  particleSkipB.GetChild("SkipPhotons").GetData(fIsPhotonSkipped);
  particleSkipB.GetChild("SkipHadrons").GetData(fIsHadronSkipped);
  particleSkipB.GetChild("SkipMesons").GetData(fIsMesonSkipped);

  Branch algoB = topB.GetChild("AlgorithmParameters");
  algoB.GetChild("DeltaROverR").GetData(fDeltaROverR);
  algoB.GetChild("DeltaPhi").GetData(fDeltaPhi);
  algoB.GetChild("HorizontalParticleCut").GetData(fHorizontalParticleCut);
  algoB.GetChild("UseStationPositionMatrix").GetData(fUseStationPositionMatrix);
  algoB.GetChild("PhiGranularity").GetData(fPhiGranularity);
  algoB.GetChild("RGranularity").GetData(fRGranularity);
  algoB.GetChild("LogGaussSmearingWidth").GetData(fLogGaussSmearingWidth);

	algoB.GetChild("SamplingDetAreaSizeX").GetData(fSamplingDetAreaSizeX);
	algoB.GetChild("SamplingDetAreaSizeY").GetData(fSamplingDetAreaSizeY);
	
  Branch limitSwitchB = topB.GetChild("ResamplingWeightLimiting", useAtt);
  if (limitSwitchB) {
    fUseWeightLimiting = true;
    limitSwitchB.GetChild("WeightThreshold").GetData(fWeightThreshold);
    limitSwitchB.GetChild("AccumulatedWeightLimit").GetData(fAccumulatedWeightLimit);
  }

	Branch saveInjParticleInfoB = topB.GetChild("SaveInjParticlesInfo", useAtt);
  if (saveInjParticleInfoB) {
    fSaveInjParticlesInfo = true;
    saveInjParticleInfoB.GetChild("OutFileName").GetData(fOutputFileName);
  }

	if(fSaveInjParticlesInfo) {
		fOutputFile= new TFile(fOutputFileName.c_str(),"RECREATE");
		fInjParticleInfo= new TTree("InjParticleInfo","InjParticleInfo");
		fInjParticleInfo->Branch("x",&x,"x/D");
		fInjParticleInfo->Branch("y",&y,"y/D");	
		fInjParticleInfo->Branch("xstat",&xstat,"xstat/D");
		fInjParticleInfo->Branch("ystat",&ystat,"ystat/D");	
		fInjParticleInfo->Branch("E",&E,"E/D");
		fInjParticleInfo->Branch("t",&t,"t/D");
		fInjParticleInfo->Branch("dx",&dx,"dx/D");		
		fInjParticleInfo->Branch("dy",&dy,"dy/D");
		fInjParticleInfo->Branch("dz",&dz,"dz/D");
		fInjParticleInfo->Branch("id",&StatId,"id/I");
	
		fStationInfo= new TTree("StationInfo","StationInfo");
		fStationInfo->Branch("x",&StatX,"x/D");
		fStationInfo->Branch("y",&StatY,"y/D");
		fStationInfo->Branch("id",&StatId,"id/I");

		fShowerInfo= new TTree("ShowerInfo","ShowerInfo");
		fShowerInfo->Branch("CoreX",&CoreX,"CoreX/D");
		fShowerInfo->Branch("CoreY",&CoreY,"CoreY/D");
		fShowerInfo->Branch("minRadiusCut",&minRadiusCut,"minRadiusCut/I");	
		fShowerInfo->Branch("maxRadiusCut",&maxRadiusCut,"maxRadiusCut/I");	
	}

  fRandomEngine = &RandomEngineRegistry::GetInstance().Get(RandomEngineRegistry::eDetector).GetEngine();

  if (!fUseStationPositionMatrix)
    WARNING("Not using StationPositionMatrix. "
            "All stations will be searched for injection.");

  return eSuccess;
}


bool
ShowerRegenerator::IsParticleEnergyLow(const int pType,
                                             const double pEnergy)
  const
{
  switch (pType) {
  case PHOTONS:
    return pEnergy < fPhotonEnergyCut;
  case ELECTRONS:
    return pEnergy < fElectronEnergyCut;
  case MUONS:
    return pEnergy < fMuonEnergyCut;
  case BARYONS:
    return pEnergy < fHadronEnergyCut;
  case MESONS:
    return pEnergy < fMesonEnergyCut;
  default:
    return true;
  }
}


bool
ShowerRegenerator::IsParticleSkipped(const int pType)
	const
{
  switch (pType) {
  case PHOTONS:
    return fIsPhotonSkipped;
  case ELECTRONS:
    return fIsElectronSkipped;
  case MUONS:
    return fIsMuonSkipped;
  case BARYONS:
    return fIsHadronSkipped;
  case MESONS:
    return fIsMesonSkipped;
  default:
    return true;
  }
}


void
ShowerRegenerator::InitNewShower(Event& event)
{

	const CoordinateSystemPtr siteCS = det::Detector::GetInstance().GetSiteCoordinateSystem();

  fShowerData = new ShowerData(StationPositionMatrix(fPhiGranularity, fRGranularity));
  fShowerData->fWeightCounterMap.clear();

  if (fLogGaussSmearingWidth) fShowerData->fLogGauss= new LogGaussSmearing(fLogGaussSmearingWidth, fRandomEngine);

  const ShowerSimData& simShower = event.GetSimShower();

  fShowerData->fParticleIt = simShower.GroundParticlesBegin();
  fShowerData->fParticlesEnd = simShower.GroundParticlesEnd();

  const double showerZenith = simShower.GetZenith();
  const double samplingAreaFactor= 4 * fDeltaPhi * fDeltaROverR / cos(showerZenith);

	CoreX= simShower.GetPosition().GetX(siteCS);
	CoreY= simShower.GetPosition().GetY(siteCS);
	minRadiusCut= simShower.GetMinRadiusCut(); 	
	maxRadiusCut= simShower.GetMaxRadiusCut();

	if(fSaveInjParticlesInfo)
		fShowerInfo->Fill();

  if (!event.HasSEvent()) event.MakeSEvent();
  SEvent& sEvent = event.GetSEvent();

  const CoordinateSystemPtr showerCS = simShower.GetShowerCoordinateSystem();
  const Point& corePosition = simShower.GetPosition();
  const TimeStamp& coreTime = simShower.GetTimeStamp();
  //const mdet::MDetector& mDetector = MyDetector::GetInstance().GetMDetector();
	const sdet::SDetector& sDetector = det::Detector::GetInstance().GetSDetector();
	

  int nStations = 0;
  int nInner = 0;
  int nOuter = 0;

  //for (mdet::MDetector::StationIterator sdIt = mDetector.StationsBegin();
  //     sdIt != mDetector.StationsEnd(); ++sdIt) {

	for (sdet::SDetector::StationIterator sdIt = sDetector.StationsBegin();
       sdIt != sDetector.StationsEnd(); ++sdIt) {

    ++nStations;

    const int sId = sdIt->GetId();

    if (!sEvent.HasStation(sId))
      sEvent.MakeStation(sId);

    //MStation& station = sEvent.GetStation(sId);
		Station& station = sEvent.GetStation(sId);

		//save station info
		StatId= sId;
		StatX= sdIt->GetPosition().GetX(siteCS);
		StatY= sdIt->GetPosition().GetY(siteCS);
		if(fSaveInjParticlesInfo)
			fStationInfo->Fill();


    // clear particle lists for the station. The stuff above
    // presumably only needs to be called once per shower, but
    // particle clearing has to happen each time through.
    if (!station.HasSimData()) station.MakeSimData();
    StationSimData& sSim = station.GetSimData();

    const Point& sPosition = sdIt->GetPosition();

    const TimeStamp planeTime = coreTime + TimeInterval(PlaneFrontTime(showerCS, corePosition, sPosition));
    
		#if (OFFLINE_MINOR_VERSION>5)
			sSim.SetPlaneFrontTime(planeTime);
		#else
			sSim.SetTime(planeTime);
		#endif
		
    const double sR = sPosition.GetRho(showerCS);

    if (sR < fInnerRadiusCut) {
      ++nInner;
      continue;
    }

    if (sR > fOuterRadiusCut) {
      ++nOuter;
      continue;
    }

    // DV: Particle weights are calculated by the shower simulation according
    // to the particle densities as measured in the ground-particle plane.
    // Thus all reweighting area ratios calculated in the resampler have to
    // have area values expressed in the corresponding ground-particle plane
    // projections, ie if there is a small missalignement of the tank vertical
    // and the ground-particle plane vertical, it has to be taken into account
    // with appropriate cosine of this small angle.
    // First calculate the ground-particle area of the resampling pie
    // section of a station:
    const double samplingArea = samplingAreaFactor * Sqr(sR);

    // it is best to create station matrix based on log(r^2)
    const double r1 = 2 * log(sR * (1 - fDeltaROverR));
    const double r2 = 2 * log(sR * (1 + fDeltaROverR));

    fShowerData->fStationMatrix.PushBack(*sdIt, station,
                                         sPosition.GetPhi(showerCS), fDeltaPhi,
                                         sR, r1, r2,
                                         samplingArea);

  }//end loop over stations

  fShowerData->fStationMatrix.CreateMatrix(fUseStationPositionMatrix);
  fShowerData->fMinR2 = exp(fShowerData->fStationMatrix.GetMinR());

  ostringstream info;
  info << "Out of " << nStations << " stations: inner cut " << nInner << ", outer " << nOuter;
  INFO(info);
}


VModule::ResultFlag
ShowerRegenerator::Run(Event& event)
{
	
  if (!event.HasSimShower()) {
    ERROR("Current event does not have a simulated shower.");
    return eFailure;
  }

  const ShowerSimData& simShower = event.GetSimShower();

  // sim event finished ?
  if (fShowerData && fShowerData->fParticleIt == fShowerData->fParticlesEnd) {
    fShowerData = 0;
    // break inner loop
    return eBreakLoop;
  }

  if (!fShowerData)
    InitNewShower(event);

  {
    ostringstream info;
    info << "Regenerating batch of "
         << fMaxParticles << " particles.";
    INFO(info);
  }

  const sdet::SDetector& sDetector = det::Detector::GetInstance().GetSDetector();
  SEvent& sEvent = event.GetSEvent();

  // clear particle list for the stations
  for (sdet::SDetector::StationIterator sdIt = sDetector.StationsBegin();
       sdIt != sDetector.StationsEnd(); ++sdIt)
    sEvent.GetStation(sdIt->GetId()).GetSimData().ClearParticleList();

  const CoordinateSystemPtr showerCS = simShower.GetShowerCoordinateSystem();
  const Point& corePosition = simShower.GetPosition();
  const CoordinateSystemPtr grParticleCS = simShower.GetLocalCoordinateSystem();
	const CoordinateSystemPtr siteCS = det::Detector::GetInstance().GetSiteCoordinateSystem();

	
  // Based on the cuts read in for the outer and inner radius in Init(), we
  // build the list of stations that lie within the region of interest in
  // the shower plane. We then loop only over these stations, since anything
  // outside the limits defined by the outer and inner radius cuts will never
  // receive any regenerated particles anyway.

  // keep count to avoid excessive memory consumption
  unsigned int nParticlesInjected = 0;

  ShowerParticleIterator& particleIt = fShowerData->fParticleIt;

  for (; particleIt != fShowerData->fParticlesEnd; ++particleIt) {

    // over the cache limit?
    if (nParticlesInjected >= fMaxParticles) {
      ostringstream info;
      info << nParticlesInjected << " particles processed.";
      INFO(info);
      return eSuccess;
    }

    const Point& pPosition = particleIt->GetPosition();
    const double pR2ShowerFrame = pPosition.GetRho2(showerCS);

    if (pR2ShowerFrame <= fShowerData->fMinR2) continue;

    const double pLogR2ShowerFrame = log(pR2ShowerFrame);

    const double pPhiShowerFrame = pPosition.GetPhi(showerCS);

    const StationPositionMatrix::StationInfoPtrList& sList =
      fShowerData->fStationMatrix.GetStationList(pPhiShowerFrame, pLogR2ShowerFrame);

    if (sList.empty()) continue;

    if (IsParticleEnergyLow(particleIt->GetType(), particleIt->GetKineticEnergy())) continue;
		if (IsParticleSkipped(particleIt->GetType())) continue;

		//loop over involved stations
    for (StationPositionMatrix::StationInfoPtrList::const_iterator sIt = sList.begin();
         sIt != sList.end(); ++sIt) {

      const StationInfo& sInfo = **sIt;

      if (!sInfo.IsIn(pPhiShowerFrame, pLogR2ShowerFrame)) continue;

			//## Check if the given sampling det area at ground (used for muon counters)
			//## is larger than the resampling area. 
			//## It should be greater or at least equal to the muon counters real area and smaller than the resampling area.
			//## If not throw an error!
			//## Better do this check at the very beginning of the module
			double detAreaResampling= fSamplingDetAreaSizeX*fSamplingDetAreaSizeY;
			//cout<<"sInfo.GetResamplingArea()="<<sInfo.GetResamplingArea()/m2<<"  detAreaResampling="<<detAreaResampling/m2<<endl;
			if(sInfo.GetResamplingArea() < detAreaResampling){
				ERROR("Given sampling detector area at ground is larger than resampling area.");
    		return eFailure;   
			}

      const double pWeight = particleIt->GetWeight();

      Particle newParticle(*particleIt);

      // This decision should be taken elsewhere, after
      // calculating AvgN
      //
      //original      newParticle.SetWeight(1);

      // all particle positions within the resampling area are equivalent, define
      // "weight density", ie the probability to find a particle with such properties
      // within some area
      const double pWeightDensity = pWeight / sInfo.GetResamplingArea();

			//get SD station info
      const sdet::Station& dStation = sInfo.GetDetStation();
      const double sThickness = dStation.GetThickness();
      const double sRadius = dStation.GetRadius() + sThickness;
      const double sHeight = dStation.GetHeight() + sThickness;
			const CoordinateSystemPtr dStationCS = dStation.GetLocalCoordinateSystem();

			
			/*
			//## get corresponding muon counter info
			//## create muon counter coordinate system
			//## center at ground at the basis of the corresponding tank
			const ReferenceEllipsoid e = ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84);
   		CoordinateSystemPtr ecef = e.GetECEF();

			const UTMPoint MuonCounterStationUTM(6075445*m, 446954*m, 1444*m, 19, 'H', e);
   	 	const Point MuonCounterStation(MuonCounterStationUTM.GetPoint());

			const CoordinateSystemPtr MuonCounterCS =
      AugerCoordinateSystem::Create(MuonCounterStationUTM, siteCS);

    	const CoordinateSystemPtr MuonCounterLocalCS =
      LocalCoordinateSystem::Create(MuonCounterStation);	
			
			const Point dStationPositionInSiteCS = dStation.GetPosition();
			*/


      StationSimData& sSim = sInfo.GetEvtStation().GetSimData();


      const unsigned int sId = dStation.GetId();

      const Vector& pDirection = newParticle.GetDirection();
      const double pStationCosine = pDirection.GetZ(dStationCS);
      const double pPlaneCosine = pDirection.GetZ(grParticleCS);

      // particle time relative to the core time
      const double pTime = particleIt->GetTime().GetInterval();
			

      // top entry

      // tank top surface seen by the particle?
      if (pStationCosine < 0 && pPlaneCosine < fHorizontalParticleCut) {

        // projected top area of a tank
        //const double effArea = kPi*Sqr(sRadius) * pStationCosine/pPlaneCosine;

				//## effective area for muon counters
				const double effArea = detAreaResampling * pStationCosine/pPlaneCosine;
				
        const double avgN = pWeightDensity * effArea;

        InsertValue(fShowerData->fWeightStat, sId, avgN);

        unsigned int n;
        double weight;

        if (fUseWeightLimiting) {
          boost::tie(n, weight) = ParticleNumberAndWeight(avgN, sId);

          if (int(weight) == 1)
            ++unityWeightCtr;
          else
            largeWeightCtr += weight;

        } 
				else {
          n = RandPoisson::shoot(fRandomEngine, avgN);
          weight = 1;
        }

        newParticle.SetWeight(weight);

        for (unsigned int i=0;i<n;++i) {

					// random generation in tank circular surface
          //const double r = sRadius * sqrt(RandFlat::shoot(fRandomEngine, 0, 1));
          //const double phi = RandFlat::shoot(fRandomEngine, 0, kTwoPi);
          //const Point pShiftedPosition_tank(r, phi, sHeight, dStationCS, Point::kCylindrical);
					//const Point pShiftedPosition_tank(r, phi, 0., dStationCS, Point::kCylindrical);

					// random generation in ground area
					const double xrand= -fSamplingDetAreaSizeX/2. + fSamplingDetAreaSizeX*RandFlat::shoot(fRandomEngine, 0, 1);
					const double yrand= -fSamplingDetAreaSizeY/2. + fSamplingDetAreaSizeY*RandFlat::shoot(fRandomEngine, 0, 1);
					const Point pShiftedPosition(xrand, yrand, 0., dStationCS, Point::kCartesian);
					//const Point pShiftedPosition(xrand, yrand, sHeight, dStationCS);
					const UTMPoint pShiftedPositionUTM(pShiftedPosition,ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));

					const Point pShiftedPositionInSiteCS= pShiftedPosition;
					pShiftedPositionInSiteCS.TransformTo(siteCS);
					const UTMPoint pShiftedPositionInSiteCSUTM(pShiftedPositionInSiteCS,ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));

					const Point dStationPos= dStation.GetPosition();
					const UTMPoint dStationPosUTM(dStationPos,ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));
					
					/*
					cout<<"StationId="<<dStation.GetId()<<"  PosStatCS=("<<dStation.GetPosition().GetX(dStationCS)/m<<","<<dStation.GetPosition().GetY(dStationCS)/m<<","<<dStation.GetPosition().GetZ(dStationCS)/m<<")  PosPartStatCS=("<<pShiftedPosition.GetX(dStationCS)/m<<","<<pShiftedPosition.GetY(dStationCS)/m<<","<<pShiftedPosition.GetZ(dStationCS)/m<<")"<<endl;
					cout<<"StationId="<<dStation.GetId()<<"  PosSiteCS=("<<dStation.GetPosition().GetX(siteCS)/m<<","<<dStation.GetPosition().GetY(siteCS)/m<<","<<dStation.GetPosition().GetZ(siteCS)/m<<")  PosPartSiteCS=("<<pShiftedPosition.GetX(siteCS)/m<<","<<pShiftedPosition.GetY(siteCS)/m<<","<<pShiftedPosition.GetZ(siteCS)/m<<")"<<endl;
					cout<<"StationId="<<dStation.GetId()<<"  UTMPos=("<<dStationPosUTM.GetNorthing()<<","<<dStationPosUTM.GetEasting()<<","<<dStationPosUTM.GetHeight()<<")  PosPart=("<<pShiftedPositionUTM.GetNorthing()<<","<<pShiftedPositionUTM.GetEasting()<<","<<pShiftedPositionUTM.GetHeight()<<")"<<endl;
					cout<<"StationId="<<dStation.GetId()<<"  UTMPos=("<<dStationPosUTM.GetNorthing()<<","<<dStationPosUTM.GetEasting()<<","<<dStationPosUTM.GetHeight()<<")  PosPart=("<<pShiftedPositionInSiteCSUTM.GetNorthing()<<","<<pShiftedPositionInSiteCSUTM.GetEasting()<<","<<pShiftedPositionInSiteCSUTM.GetHeight()<<")"<<endl;
					*/
          newParticle.SetPosition(pShiftedPosition);

          // Shift the time to account for differences in position between the
          // hit tank and the position the particle was placed at the ground
          // by the shower simulation. Also, include a shift for multiple
          // particle entries from regeneration of one weighted particle. So far,
          // no one has provided any real justification for what distribution
          // should be used for the energy shifting.

					//double pShiftedTime_tank= pTime + PlaneFrontTime(showerCS, pPosition, pShiftedPosition_tank);

					double pShiftedTime= pTime + PlaneFrontTime(showerCS, pPosition, pShiftedPosition);
					//cout<<"pTime="<<pTime<<"  pShiftedTime_tank="<<pShiftedTime_tank<<"  pShiftedTime="<<pShiftedTime<<endl;         

          if (i && fShowerData->fLogGauss) {
            const double planeFrontTime = PlaneFrontTime(showerCS, corePosition, pShiftedPosition);
						//const double planeFrontTime_tank = PlaneFrontTime(showerCS, corePosition, pShiftedPosition_tank);
						//cout<<"planeFrontTime_tank="<<planeFrontTime_tank<<"  planeFrontTime="<<planeFrontTime<<endl;
            pShiftedTime = fShowerData->fLogGauss->GetSmearedTime(planeFrontTime, pShiftedTime);
          }
          newParticle.SetTime(pShiftedTime);

          InsertValue(fShowerData->fTimeStat, sId, pShiftedTime);

          sSim.AddParticle(newParticle);


					//summary inj particle info
					x= pShiftedPosition.GetX(siteCS)/m;
					y= pShiftedPosition.GetY(siteCS)/m;
					xstat= pShiftedPosition.GetX(dStationCS)/m;
					ystat= pShiftedPosition.GetY(dStationCS)/m;
					t= pShiftedTime/ns;
					dx= newParticle.GetDirection().GetX(dStationCS);
					dy= newParticle.GetDirection().GetY(dStationCS);
					dz= newParticle.GetDirection().GetZ(dStationCS);
					E= newParticle.GetKineticEnergy()/GeV;
					StatId= dStation.GetId();

					if(fSaveInjParticlesInfo && fInjParticleInfo) 
						fInjParticleInfo->Fill();

        }//end loop n particles injection

        fShowerData->fTopParticles[sId] += n;
        nParticlesInjected += n;

      } // if top entries
	
			/*
			//# comment here as for muon counter only the horizontal sampling area is present
      // side entry
      const double pStationSine = pDirection.GetRho(dStationCS);
      const double pPlaneAbsCosine = abs(pPlaneCosine);

      if (pPlaneAbsCosine < fHorizontalParticleCut) 
				++fShowerData->fNHorizontalParticles;

      else if (pStationSine) {

        // projected side area of a tank
        const double effArea = 2*sRadius*sHeight * pStationSine/pPlaneAbsCosine;
        const double avgN = pWeightDensity * effArea;
        //original      const unsigned int n = RandPoisson::shoot(fRandomEngine, avgN);
        const double pPhi = pDirection.GetPhi(dStationCS) + kPiOnTwo;

        unsigned int n;
        double weight;

        if (fUseWeightLimiting) {
          boost::tie(n, weight) = ParticleNumberAndWeight(avgN, sId);
        } 
				else {
          n = RandPoisson::shoot(fRandomEngine, avgN);
          weight = 1;
        }

        newParticle.SetWeight(weight);

        for (unsigned int i = 0; i < n; ++i) {

          // For the X and Y component we need to find the angle made between
          // the (tank) X-Y axis and the incoming particle direction, then add
          // an angle sampled from a cosine distribution.

          const double phi = pPhi + acos(RandFlat::shoot(fRandomEngine, -1, 1));
          const double z = RandFlat::shoot(fRandomEngine, -sThickness, sHeight);
          const Point pShiftedPosition(sRadius, phi, z, dStationCS, Point::kCylindrical);
          newParticle.SetPosition(pShiftedPosition);

          // Shift the time (read comment above).
          double pShiftedTime =
            pTime + PlaneFrontTime(showerCS, pPosition, pShiftedPosition);
          if (i && fShowerData->fLogGauss) {
            const double planeFrontTime =
              PlaneFrontTime(showerCS, corePosition, pShiftedPosition);
            pShiftedTime =
              fShowerData->fLogGauss->GetSmearedTime(planeFrontTime, pShiftedTime);
          }
          newParticle.SetTime(pShiftedTime);

          InsertValue(fShowerData->fTimeStat, sId, pShiftedTime);

          sSim.AddParticle(newParticle);

        fShowerData->fSideParticles[sId] += n;
        nParticlesInjected += n;
        if (pStationCosine >= 0)
          fShowerData->fUpwardSideParticles[sId] += n;

        }// end loop n particles injected on side
				
      } // if side entries
			*/

    }//loop stations

  }//loop ground particles

  ostringstream info;
  info << nParticlesInjected << " particles processed.";
  INFO(info);

  if (fShowerData->fParticleIt == fShowerData->fParticlesEnd)
    OutputStats(event);

  return eSuccess;
}


void
ShowerRegenerator::OutputStats(Event& event)
{
  const StationPositionMatrix::StationInfoList& sList =
    fShowerData->fStationMatrix.GetStationList();

  int nStations = 0;

  {
    INFO("\nStation particle statistics:");

    TabularStream tab("r|.|r|r|r|.|.|.");
    tab <<              endc << "Rho"  << endc <<          endc <<           endc <<         endc
        << "wght" << endc << "wght" << endc << "wght" << endr

        << "station" << endc << "[km]" << endc << "top" << endc << "side" << endc << "up" << endc
        << "min"  << endc << "avg"  << endc << "max"  << endr << hline;

    for (StationPositionMatrix::StationInfoList::const_iterator sIt = sList.begin();
         sIt != sList.end(); ++sIt) {

      const unsigned sId = sIt->GetDetStation().GetId();

      if (fShowerData->fTopParticles[sId] ||
          fShowerData->fSideParticles[sId] ||
          fShowerData->fUpwardSideParticles[sId]) {

        ShowerData::WeightStatMap::const_iterator wIt = fShowerData->fWeightStat.find(sId);
        ShowerData::TimeStatMap::const_iterator mmIt = fShowerData->fTimeStat.find(sId);
        tab << sId << ' ' << endc
            << setprecision(4) << sIt->GetR()/km << endc
            << ' ' << fShowerData->fTopParticles[sId] << ' ' << endc
            << ' ' << fShowerData->fSideParticles[sId] << ' ' << endc
            << ' ' << fShowerData->fUpwardSideParticles[sId] << ' ' << endc
            << ' ' << setprecision(4) << wIt->second.GetMin() << endc
            << ' ' << setprecision(4) << wIt->second.GetAverage() << endc
            << ' ' << setprecision(4) << wIt->second.GetMax() << endr;
        ++nStations;
      }

    }

    tab << delr;

    if (nStations)
      cerr << tab;
    else
      cerr << "No stations were hit." << endl;
  }

  {
    INFO("\nStation timing statistics:\n"
         "abs = absolute time difference to core time\n"
         "rel = relative time to station plane-front arrival");
    TabularStream tab("r|.|.|.|.|.");
    tab <<              endc << "Rho"  << endc << "min abs"   << endc << "max abs"   << endc
        << "min rel"   << endc << "max rel"   << endr

        << "station" << endc << "[km]" << endc << "time [ns]" << endc << "time [ns]" << endc
        << "time [ns]" << endc << "time [ns]" << endr << hline;

    const ShowerSimData& simShower = event.GetSimShower();
    const TimeStamp& coreTime = simShower.GetTimeStamp();

    for (StationPositionMatrix::StationInfoList::const_iterator sIt = sList.begin();
         sIt != sList.end(); ++sIt) {

      const unsigned sId = sIt->GetDetStation().GetId();

      if (fShowerData->fTopParticles[sId] ||
          fShowerData->fSideParticles[sId] ||
          fShowerData->fUpwardSideParticles[sId]) {

        //const TimeInterval sTime = sIt->GetEvtStation().GetSimData().GetPlaneFrontTime() - coreTime;
				TimeInterval sTime;
				#if (OFFLINE_MINOR_VERSION>5)
					sTime = sIt->GetEvtStation().GetSimData().GetPlaneFrontTime() - coreTime;
				#else
					sTime = sIt->GetEvtStation().GetSimData().GetTime() - coreTime;
				#endif


        ShowerData::WeightStatMap::const_iterator wIt = fShowerData->fWeightStat.find(sId);
        ShowerData::TimeStatMap::const_iterator mmIt = fShowerData->fTimeStat.find(sId);
        tab << sId << ' ' << endc
            << setprecision(4) << sIt->GetR()/km << endc
            << ' ' << Round(10, mmIt->second.GetMin().GetInterval()/ns) << endc
            << ' ' << Round(10, mmIt->second.GetMax().GetInterval()/ns) << endc
            << ' ' << Round(10, (mmIt->second.GetMin() - sTime).GetInterval()/ns) << endc
            << ' ' << Round(10, (mmIt->second.GetMax() - sTime).GetInterval()/ns) << endr;

      }

    }

    tab << delr;

    if (nStations)
      cerr << tab;
  }

  if (nStations) {
    ostringstream info;
    info << fShowerData->fNHorizontalParticles << " horizontal particles rejected.";
    INFO(info);
  }
}


VModule::ResultFlag
ShowerRegenerator::Finish()
{
  if (fNOutOfRangeWeight) {
    ostringstream warn;
    warn << "During shower regeneration, a total of " << fNOutOfRangeWeight << " "
            "weighted particle(s) resulted in generation of unity-weight particles exceeding "
         << fWeightLimit << " in number.";
    WARNING(warn);
  }

	if(fSaveInjParticlesInfo){
		fOutputFile->cd();
		fShowerInfo->Write();
		fStationInfo->Write();
		fInjParticleInfo->Write();
		fOutputFile->Close();
	}
		

//   cout << "unityWeightCtr = " << unityWeightCtr << endl;
//   cout << "largeWeightCtr = " << largeWeightCtr << endl;

  return eSuccess;
}

  /// (Optional) special handling for particles with very large weights.
  /**
     EXPLAIN HOW THIS WORKS
     need for arguments..
   */
  boost::tuple<unsigned int, double>
  ShowerRegenerator::ParticleNumberAndWeight(const double& avgN, const int& sId)
  {
    unsigned int n;
    double weight;

    if (avgN > fWeightThreshold) {

      //      fShowerData->fWeightCounterMap[sId] += avgN;

      if (fShowerData->fWeightCounterMap[sId] <= fAccumulatedWeightLimit) {

//         if (sId == 5669)
//           cout << "a (avgN = " << avgN << ") ";

        n = RandPoisson::shoot(fRandomEngine, avgN);
        weight = 1;
      } else  {

//         if (sId == 5669)
//           cout << "b (avgN = " << avgN << ") ";

        n = 1;
        weight = avgN;
      }
    }
    else {
      n = RandPoisson::shoot(fRandomEngine, avgN);
      weight = 1;
    }

    fShowerData->fWeightCounterMap[sId] += n;
    return boost::make_tuple(n, weight);
  }

