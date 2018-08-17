#include <utl/Accumulator.h>
#include "StationPositionMatrix.h"

using namespace std;
using namespace utl;


namespace G4MuonCounterShowerRegenerator {

  void
  StationInfo::Dump()
    const
  {
    std::cerr
      << "ds " << fDetStation->GetId() << " es " << fEvtStation->GetId() << " "
         "phi " << fPhi1 << ' ' << fPhi2 << " "
         "r " << fR1 << ' ' << fR2 << " "
         "a " << fResamplingArea << "\n";
  }


  void
  StationPositionMatrix::Clear()
  {
    fR1 = -2;
    fRStep = -1;
    fStations.clear();
    fStationMatrix.clear();
    fEmptyStationInfoPtrList.clear();
  }


  void
  StationPositionMatrix::CreateMatrix(const bool useSpatialStationMatrix)
  {
    const double tpi = 2*utl::kPi;
    if (!useSpatialStationMatrix || fStations.empty()) {
      // ensure negative index for all possible r
      fR1 = -2;
      fRStep = -1;
      // "empty" list of stations will be returned
      fEmptyStationInfoPtrList.clear();
      for (StationInfoList::const_iterator it = fStations.begin();
           it != fStations.end(); ++it)
        fEmptyStationInfoPtrList.push_back(&(*it));
      return;
    }
    // find range
    Accumulator::Min<double> r1(fStations[0].GetR1());
    Accumulator::Max<double> r2(fStations[0].GetR2());
    Accumulator::Average deltaPhi;
    Accumulator::Average deltaR;
    for (unsigned int i = 1; i < fStations.size(); ++i) {
      r1(fStations[i].GetR1());
      r2(fStations[i].GetR2());
      deltaR(fStations[i].GetR2() - fStations[i].GetR1());
      const double dphi = fStations[i].GetPhi2() - fStations[i].GetPhi1();
      deltaPhi((dphi >= 0) ? dphi : tpi-dphi);
    }
    fNPhi = int(tpi / deltaPhi.GetAverage() * fPhiGranularity + 0.5);
    fPhiStep = tpi / fNPhi;
    const double dr = r2.GetMax() - r1.GetMin();
    fNR = int(dr / deltaR.GetAverage() * fRGranularity + 0.5);
    fR1 = r1.GetMin();
    fRStep = dr / fNR;

    Resize(fNPhi, fNR);

    // fill matrix with station info pointers
    for (StationInfoList::const_iterator it = fStations.begin();
         it != fStations.end(); ++it) {
      const int iPhi1 = GetPhiIndex(it->GetPhi1());
      const int iPhi2 = GetPhiIndex(it->GetPhi2());
      const int iR1 = GetRIndex(it->GetR1());
      int iR2 = GetRIndex(it->GetR2());
      if (iR2 == fNR)
        iR2 = fNR - 1;
      if (iPhi1 <= iPhi2)
        for (int i = iPhi1; i <= iPhi2; ++i)
          for (int j = iR1; j <= iR2; ++j)
          fStationMatrix[i][j].push_back(&(*it));
      else
        for (int j = iR1; j <= iR2; ++j) {
          for (int i = 0; i <= iPhi2; ++i)
            fStationMatrix[i][j].push_back(&(*it));
          for (int i = iPhi1; i < fNPhi; ++i)
            fStationMatrix[i][j].push_back(&(*it));
        }
    }
  }


  void
  StationPositionMatrix::DumpStats()
    const
  {
    Accumulator::MinMax<int> minmax(fStationMatrix[0][0].size());
    Accumulator::Average size;
    for (int i = 0; i < fNPhi; ++i)
      for (int j = 0; j < fNR; ++j) {
        const int s = fStationMatrix[i][j].size();
        minmax(s);
        size(s);
      }
    cerr
      << "\n"
         "nPhi = " << fNPhi << "\n"
         "R = " << fR1 << " : " << fR1 + fNR * fRStep << ", " << fNR << "\n"
         "matrix = " << minmax.GetMin() << ", " << size.GetAverage() << ", "
      << minmax.GetMax() << '\n';
    for (StationInfoList::const_iterator it = fStations.begin();
         it != fStations.end(); ++it)
      it->Dump();
  }


  void
  StationPositionMatrix::Resize(const int nPhi, const int nR)
  {
    fStationMatrix.clear();
    fStationMatrix.resize(nPhi);
    for (StationMatrix::iterator it = fStationMatrix.begin();
         it != fStationMatrix.end(); ++it)
      it->resize(nR);
  }

}
