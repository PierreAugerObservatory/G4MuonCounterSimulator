/**
* @file Utilities.cc
* @class Utilities
* @brief Utility functions for 
*
* Useful functions to be used in the simulation, or to handle output files
* @author Dr. Simone Riggi
* @date 23/08/2010
*/


#ifndef Utilities_hh
#define Utilities_hh 1

#include "AnalysisConsts.hh"

#include <vector>
#include <algorithm>
#include <map>
#include <string>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TGraph.h>
#include <TVector3.h>

using namespace std;

class Utilities{

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Utilities();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~Utilities();

		
	public:

		/**
		* \brief Merge simulation ROOT file given in a list
		*/
    void MergeFiles(const char * fileListName);
		
		/**
		* \brief Convert photon energy (in eV) into wavelength (nm)
		*/
    static double Energy2Wavelength(double energy)
		{
			energy*= 1./Joule2eV;//convert energy in Joule
  		double lambda= 1.E+9* PlanckConstant*SpeedOfLight/energy;//in nm
  		return lambda;		
		};

		/**
		* \brief Convert photon wavelength (nm) into energy (eV)
		*/
    static double Wavelength2Energy(double lambda){	
			lambda*= 1.E-09;//convert lambda in m
  		double E= Joule2eV* PlanckConstant*SpeedOfLight/lambda;//in eV 
 		 	return E;
		};

		/**
		* \brief Order vectors and get ordering index
		*/
		template<class T> struct index_cmp{

  		index_cmp(const T arr) : arr(arr) {}
  		bool operator()(const size_t a, const size_t b) const
 			{
    		return arr[a] < arr[b];
  		}
  		const T arr;
		};

		template< class T >
			static void reorder(std::vector<T> & unordered,std::vector<size_t> const & index_map,std::vector<T> & ordered){
  			// copy for the reorder according to index_map, because unsorted may also be
  			// sorted
  			std::vector<T> copy = unordered;
  			ordered.resize(index_map.size());
  			for(int i = 0; i<index_map.size();i++)
					ordered[i] = copy[index_map[i]];
			}

		template <class T>
			static void sort(std::vector<T> & unsorted,std::vector<T> & sorted,std::vector<size_t> & index_map){
  			// Original unsorted index map
  			index_map.resize(unsorted.size());
 				for(size_t i=0;i<unsorted.size();i++)
					index_map[i] = i;
  
  			// Sort the index map, using unsorted for comparison
  			std::sort(index_map.begin(),index_map.end(),index_cmp<std::vector<T>& >(unsorted));
  			sorted.resize(unsorted.size());
  			reorder(unsorted,index_map,sorted);
			}

		static std::string ExecSystemCommand(const char* cmd);

	private:
	
		
	

};

#endif /*Utilities_hh*/
