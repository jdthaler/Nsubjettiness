//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

//NEW FILE CREATED BY TJW 12/25
//Update to move AxesFinder class and definitions into separate .cc/.hh files

#ifndef __FASTJET_CONTRIB_AXESFINDER_HH__
#define __FASTJET_CONTRIB_AXESFINDER_HH__

#include "WinnerTakeAllRecombiner.hh" //new file added by TJW 12/28
#include "MeasureFunction.hh" //new file added by TJW 12/25

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Parameter classes
//
///////

//KMeansParameters class moved from MeasureFunction.hh -- TJW 12/31
//------------------------------------------------------------------------
/// \class KmeansParameters
// Parameters that change minimization procedure. They are used in the functions that find axes by Kmeans minimization.
// They are set automatically when you choose NsubAxesMode, but can be adjusted manually as well -- comment added by TJ
class KmeansParameters {
private:
   int _n_iterations;  // Number of iterations to run  (0 for no minimization, 1 for one-pass, >>1 for global minimum)
   double _precision;  // Desired precision in axes alignment
   int _halt;          // maximum number of steps per iteration
   double _noise_range;// noise range for random initialization
   
public:
   KmeansParameters() : _n_iterations(0), _precision(0.0), _halt(0), _noise_range(0.0) {}
   KmeansParameters(const int my_n_iterations, double my_precision, int my_halt, double my_noise_range) :
   _n_iterations(my_n_iterations),  _precision(my_precision), _halt(my_halt), _noise_range(my_noise_range) {}
   int n_iterations() const { return _n_iterations;}
   double precision() const {return _precision;}
   int halt() const {return _halt;}
   double noise_range() const {return _noise_range;}
};

///////
//
// Axes Finder Options
//
///////

//AxesFinder base class definition moved from Njettiness.hh -- TJW 12/25

//------------------------------------------------------------------------
/// \class AxesFinder
// This is the base class for all axes finders. These axes are used along with the MeasureFunctions to calculate 
// tau_N. There are different implementations of axes finding that are defined in derived classes below. -- comment added by TJ 
class AxesFinder {

   protected:
      AxesFinder() {}
      
   public:
      virtual ~AxesFinder(){}

      virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets, const std::vector<fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& currentAxes) = 0;
      
};

//AxesFinder from exlusive jets and derived class definitions moved from Njettiness.hh -- TJW 12/28
//------------------------------------------------------------------------
/// \class AxesFinderFromExclusiveJetDefinition
// This class finds axes by clustering the particles and then finding the exclusive jets. This can be implemented
// with different jet algorithms. -- comment added by TJW
class AxesFinderFromExclusiveJetDefinition : public AxesFinder {

   private:
      fastjet::JetDefinition _def;
   
   public:
      AxesFinderFromExclusiveJetDefinition(fastjet::JetDefinition def) : _def(def) {}
      
      virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& currentAxes) {
         fastjet::ClusterSequence jet_clust_seq(inputs, _def);
         return jet_clust_seq.exclusive_jets(n_jets);
      }
};

//------------------------------------------------------------------------
/// \class AxesFinderFromWTA_KT
// This class finds axes by finding the exlusive jets after clustering according to a kT algorithm and a 
// winner take all recombination scheme. -- comment added by TJW
class AxesFinderFromWTA_KT : public AxesFinderFromExclusiveJetDefinition { 
   private: 
      const WinnerTakeAllRecombiner *recomb;
   public:
      AxesFinderFromWTA_KT() : AxesFinderFromExclusiveJetDefinition(
         fastjet::JetDefinition(fastjet::kt_algorithm, 
         M_PI/2.0, //TODO: Want to use maximum jet radius constant here
         recomb = new WinnerTakeAllRecombiner(), 
         fastjet::Best)) {}
      ~AxesFinderFromWTA_KT() {delete recomb;}
   };
   
//------------------------------------------------------------------------
/// \class AxesFinderFromWTA_CA
// This class finds axes by finding the exlusive jets after clustering according to a CA algorithm and a 
// winner take all recombination scheme. -- comment added by TJW
class AxesFinderFromWTA_CA : public AxesFinderFromExclusiveJetDefinition {
   private: 
      const WinnerTakeAllRecombiner *recomb;
   public:
      AxesFinderFromWTA_CA() : AxesFinderFromExclusiveJetDefinition(
         fastjet::JetDefinition(fastjet::cambridge_algorithm, 
         M_PI/2.0, //TODO: Want to use maximum jet radius constant here
         recomb = new WinnerTakeAllRecombiner(), 
         fastjet::Best)) {}
      ~AxesFinderFromWTA_CA() {delete recomb;}
};

//------------------------------------------------------------------------
/// \class AxesFinderFromKT
// This class finds axes by finding the exlusive jets after clustering according to a kT algorithm and a 
// E_scheme recombination. -- comment added by TJW   
class AxesFinderFromKT : public AxesFinderFromExclusiveJetDefinition {
   public:
      AxesFinderFromKT() : AxesFinderFromExclusiveJetDefinition(
         fastjet::JetDefinition(fastjet::kt_algorithm,
         M_PI/2.0,
         fastjet::E_scheme,
         fastjet::Best)) {}
};

//------------------------------------------------------------------------
/// \class AxesFinderFromKT
// This class finds axes by finding the exlusive jets after clustering according to a CA algorithm and a 
// E_scheme recombination. -- comment added by TJW      
class AxesFinderFromCA : public AxesFinderFromExclusiveJetDefinition { 
   public:
      AxesFinderFromCA() : AxesFinderFromExclusiveJetDefinition(
         fastjet::JetDefinition(fastjet::cambridge_algorithm,
         M_PI/2.0,
         fastjet::E_scheme,
         fastjet::Best)) {}
};


//AxesFinder from hardest jets and derived class definitions moved from Njettiness.hh -- TJW 12/28
//------------------------------------------------------------------------
/// \class AxesFinderFromHardestJetDefinition
// This class finds axes by clustering the particles and then finding the n hardest inclusive jets. 
// This can be implemented with different jet algorithms. -- comment added by TJW
class AxesFinderFromHardestJetDefinition : public AxesFinder {

   private:
      fastjet::JetDefinition _def;
   
   public:
      AxesFinderFromHardestJetDefinition(fastjet::JetDefinition def) : _def(def) {}
      
      virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& currentAxes) {
         fastjet::ClusterSequence jet_clust_seq(inputs, _def);
         std::vector<fastjet::PseudoJet> myJets = sorted_by_pt(jet_clust_seq.inclusive_jets());
         myJets.resize(n_jets);  // only keep n hardest
         return myJets;
      }      
};

//------------------------------------------------------------------------
/// \class AxesFinderFromAntiKT
// This class finds axes by finding the n hardest jets after clustering the particles according 
// to an anti kT algorithm and E_scheme. -- comment added by TJW
class AxesFinderFromAntiKT : public AxesFinderFromHardestJetDefinition {
   public:
      AxesFinderFromAntiKT(double R0) : AxesFinderFromHardestJetDefinition(fastjet::JetDefinition(fastjet::antikt_algorithm,R0,fastjet::E_scheme,fastjet::Best)) {}
};


//Manual Axes Finder class moved from Njettiness.hh -- TJW 12/28
//------------------------------------------------------------------------
/// \class AxesFinderFromUserInput
// This class allows the user to manually define the axes. -- comment added by TJW
class AxesFinderFromUserInput : public AxesFinder {

   public:
      AxesFinderFromUserInput() {}
      ~AxesFinderFromUserInput() {}
      
      virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& currentAxes) {
         assert(currentAxes.size() == (unsigned int) n_jets);
         return currentAxes;
      }
};

//This is a helper class for the Minimum Axes Finders. It is defined later.
class LightLikeAxis;                                          

/// Minimum Axes moved from Njettiness.hh -- TJW 12/28
//------------------------------------------------------------------------
/// \class AxesFinderFromKmeansMinimization
// This class finds finds axes by using Kmeans clustering to minimizaiton N-jettiness. Given a first set of 
// starting axes, it updates n times to get as close to the global minimum as possible. -- comment added by TJW
class AxesFinderFromKmeansMinimization : public AxesFinder {

   private:
      AxesFinder* _startingFinder;
      KmeansParameters _paraKmeans;
      
      MeasureFunction* _function;

      // use separate beta, R0, Rcutoff since NsubParameters is no longer used -- TJW 1/9
      // R0 removed since it is unnecessary -- TJW 1/9
      double _beta;
      double _Rcutoff;
      
   public:

      //updated constructor to use three separate parameters instead of NsubParameters in definition of DefaultNormalizedMeasure -- TJW 1/9
      AxesFinderFromKmeansMinimization(AxesFinder* startingFinder, KmeansParameters paraKmeans, double beta, double Rcutoff)
         : _startingFinder(startingFinder), _paraKmeans(paraKmeans), _beta(beta), _Rcutoff(Rcutoff) {
         //_function changed to unnormalized because minimization is independent of R0 -- TJW 1/11
         _function = new DefaultUnnormalizedMeasure(beta, Rcutoff); 
      }
      
      ~AxesFinderFromKmeansMinimization() {
         delete _startingFinder;  //TODO: Convert to smart pointers to avoid this.
         delete _function;
      }

      virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& currentAxes);

      //updated function arguments to use three separate parameters instead of NsubParameters-- TJW 1/9
      template <int N> std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes, 
                                  const std::vector <fastjet::PseudoJet> & inputJets,
                                  double beta, double Rcutoff, double precision);

      //updated function arguments to use three separate parameters instead of NsubParameters-- TJW 1/9                                  
      std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes, 
                                      const std::vector <fastjet::PseudoJet> & inputJets, double beta, double Rcutoff, double precision);
      
};

// class added by TJW 1/10
//------------------------------------------------------------------------
/// \class AxesFinderFromOnePassMinimization
// This class defines an AxesFinder that uses Kmeans minimization, but only on a single pass. It inherits from the class above and defines a set of
// KMeansParameters such that only one pass is made. -- comment added by TJW
class AxesFinderFromOnePassMinimization : public AxesFinderFromKmeansMinimization {

   // noise parameter changed to NAN -- TJW 1/13
   public:
      AxesFinderFromOnePassMinimization(AxesFinder* startingFinder, double beta, double Rcutoff) 
      : AxesFinderFromKmeansMinimization(startingFinder, KmeansParameters(1,0.0001,1000,NAN), beta, Rcutoff) {}

};

//------------------------------------------------------------------------
/// \class AxesFinderFromGeometricMinimization
// This class finds finds axes by minimizing the distance between axes and particles. Given a first set of starting axes, 
// it updates n times to get as close to the global minimum as possible. -- comment added by TJW
class AxesFinderFromGeometricMinimization : public AxesFinder {

   private:
      AxesFinder* _startingFinder;
      MeasureFunction* _function;
      double _Rcutoff;
      double _nAttempts;
      double _accuracy;

   
   public:
      AxesFinderFromGeometricMinimization(AxesFinder* startingFinder, double Rcutoff) : _startingFinder(startingFinder), _Rcutoff(Rcutoff) {
         _nAttempts = 100;
         _accuracy = 0.000000001;
         _function = new GeometricMeasure(_Rcutoff);
      }

      ~AxesFinderFromGeometricMinimization() {
         delete _startingFinder;  //TODO: Convert to smart pointers to avoid this.
         delete _function;
      }

      //definition of getAxes moved to AxesFinder.cc, since it is a large function -- TJW 12/28      
      virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& currentAxes);
};

//------------------------------------------------------------------------
/// \class LightLikeAxis
// This is a helper class for the minimum Axes Finders classes above. It creates a convenient way of defining axes
// in order to better facilitate calculations. -- comment added by TJW
class LightLikeAxis {
private:
   double _rap, _phi, _weight, _mom;
   
   double DistanceSq(double rap2, double phi2) const {
      double rap1 = _rap;
      double phi1 = _phi;
      
      double distRap = rap1-rap2;
      double distPhi = std::fabs(phi1-phi2);
      if (distPhi > M_PI) {distPhi = 2.0*M_PI - distPhi;}
      return sq(distRap) + sq(distPhi);
   }
   
   double Distance(double rap2, double phi2) const {
      return std::sqrt(DistanceSq(rap2,phi2));
   }
   
   
public:
   LightLikeAxis() : _rap(0.0), _phi(0.0), _weight(0.0), _mom(0.0) {}
   LightLikeAxis(double my_rap, double my_phi, double my_weight, double my_mom) :
   _rap(my_rap), _phi(my_phi), _weight(my_weight), _mom(my_mom) {}
   double rap() const {return _rap;}
   double phi() const {return _phi;}
   double weight() const {return _weight;}
   double mom() const {return _mom;}
   void set_rap(double my_set_rap) {_rap = my_set_rap;}
   void set_phi(double my_set_phi) {_phi = my_set_phi;}
   void set_weight(double my_set_weight) {_weight = my_set_weight;}
   void set_mom(double my_set_mom) {_mom = my_set_mom;}
   void reset(double my_rap, double my_phi, double my_weight, double my_mom) {_rap=my_rap; _phi=my_phi; _weight=my_weight; _mom=my_mom;}

   // floating function moved to LightLikeAxis class -- TJW 1/14
   fastjet::PseudoJet ConvertToPseudoJet(LightLikeAxis axis);
   
   double DistanceSq(const fastjet::PseudoJet& input) const {
      return DistanceSq(input.rap(),input.phi());
   }

   double Distance(const fastjet::PseudoJet& input) const {
      return std::sqrt(DistanceSq(input));
   }

   double DistanceSq(const LightLikeAxis& input) const {
      return DistanceSq(input.rap(),input.phi());
   }

   double Distance(const LightLikeAxis& input) const {
      return std::sqrt(DistanceSq(input));
   }

};

} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_AXESFINDER_HH__
