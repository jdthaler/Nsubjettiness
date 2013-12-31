//NEW FILE CREATED BY TJW 12/25
//Update to move AxesFinder class and definitions into separate .cc/.hh files

//legal info below copied directly from Njettiness.hh

//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-13
//  Jesse Thaler, Ken Van Tilburg, and Christopher K. Vermilion
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
      NsubParameters _paraNsub;
      
      MeasureFunction* _functor;
      
   public:
      AxesFinderFromKmeansMinimization(AxesFinder* startingFinder, KmeansParameters paraKmeans, NsubParameters paraNsub)
         : _startingFinder(startingFinder), _paraKmeans(paraKmeans), _paraNsub(paraNsub) {
         _functor = new DefaultMeasure(paraNsub);
      }
      
      ~AxesFinderFromKmeansMinimization() {
         delete _startingFinder;  //TODO: Convert to smart pointers to avoid this.
         delete _functor;
      }

      virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& currentAxes);

      template <int N> std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes, 
                                  const std::vector <fastjet::PseudoJet> & inputJets,
                                  NsubParameters paraNsub, double precision);
                                  
      std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes, 
                                      const std::vector <fastjet::PseudoJet> & inputJets, NsubParameters paraNsub, double precision);
      
};

//------------------------------------------------------------------------
/// \class AxesFinderFromGeometricMinimization
// This class finds finds axes by minimizing the distance between axes and particles. Given a first set of starting axes, 
// it updates n times to get as close to the global minimum as possible. -- comment added by TJW
class AxesFinderFromGeometricMinimization : public AxesFinder {

   private:
      AxesFinder* _startingFinder;
      MeasureFunction* _functor;
      double _Rcutoff;
      double _nAttempts;
      double _accuracy;

   
   public:
      AxesFinderFromGeometricMinimization(AxesFinder* startingFinder, double Rcutoff) : _startingFinder(startingFinder), _Rcutoff(Rcutoff) {
         _nAttempts = 100;
         _accuracy = 0.000000001;
         _functor = new GeometricMeasure(_Rcutoff);
      }

      ~AxesFinderFromGeometricMinimization() {
         delete _startingFinder;  //TODO: Convert to smart pointers to avoid this.
         delete _functor;
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