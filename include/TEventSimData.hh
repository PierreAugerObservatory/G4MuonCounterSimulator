/**
* @file TEventSimData.hh
* @class TEventSimData
* @brief Define the sim data structure for each simulated event
*
* A ROOT dictionary for this class is generated by the Makefile
* @author Simone Riggi
* @date 22/12/2010
*/

#ifndef TEventSimData_h
#define TEventSimData_h 1

#include "TStationSimData.hh"
#include "TScintHit.hh"
#include "TPMTHit.hh"

#include <Rtypes.h>
#include <TObject.h>
#include <TVector3.h>

#include <iostream>
#include <vector>


class TEventSimData : public TObject {

	public:
   
		/** 
		\brief Class constructor: initialize structures.
 		*/
  	TEventSimData();
		/** 
		\brief Class destructor: free allocated memory.
 		*/
    //virtual ~TEventSimData();
		~TEventSimData();

		/**
		* \brief Add a station sim data to the list 
		*/
		void AddSimStationData(TStationSimData data) {fStationSimDataCollection.push_back(data);}
		/**
		* \brief Add a station sim data to the list 
		*/
		bool AddSimStation(TStationSimData station);
		/**
		* \brief Get sim data station of given id
		*/
		TStationSimData* GetSimStation(int id);

		/**
		* \brief Collection of ROOT station sim data
		*/
		std::vector<TStationSimData> fStationSimDataCollection;
  	
  
  	ClassDef(TEventSimData,1)
};

//## PMT HIT
#ifdef __MAKECINT__
#pragma link C++ class TPMTHit+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TPMTHit*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TPMTHit>+;
#endif

//## SCINT HIT
#ifdef __MAKECINT__
#pragma link C++ class TScintHit+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TScintHit*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TScintHit>+;
#endif

//## TRACK POINT
#ifdef __MAKECINT__
#pragma link C++ class TrackPoint+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TrackPoint*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TrackPoint>+;
#endif

//## TParticleSimData
#ifdef __MAKECINT__
#pragma link C++ class TParticleSimData+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TParticleSimData*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TParticleSimData>+;
#endif

//## TStationSimData
#ifdef __MAKECINT__
#pragma link C++ class TStationSimData+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TStationSimData*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TStationSimData>+;
#endif

//## TEventSimData
#ifdef __MAKECINT__
#pragma link C++ class TEventSimData+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TEventSimData*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TEventSimData>+;
#endif


#endif
