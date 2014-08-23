//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id$
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

#include "WinnerTakeAllRecombiner.hh"
#include "MeasureFunction.hh"

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
// Starting Axes Finders
//
///////
  
//------------------------------------------------------------------------
/// \class StartingAxesFinder
// This is the base class for all starting axes finders.  Starting axes finders do not need
// seed axes and typical are based on sequential jet algorithms
class StartingAxesFinder {
   
public:
   
   // This function should be overloaded, and updates the seedAxes to return new axes
   virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets,
                                                   const std::vector<fastjet::PseudoJet>& inputs) const = 0;
   
   //virtual destructor
   virtual ~StartingAxesFinder(){}
   
};
   
//------------------------------------------------------------------------
/// \class AxesFinderFromExclusiveJetDefinition
// This class finds axes by clustering the particles and then finding the exclusive jets. This can be implemented
// with different jet algorithms.
class AxesFinderFromExclusiveJetDefinition : public StartingAxesFinder {
   
public:
   AxesFinderFromExclusiveJetDefinition(fastjet::JetDefinition def)
   : _def(def) {}
   
   virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets,
                                                   const std::vector <fastjet::PseudoJet> & inputs) const {
      fastjet::ClusterSequence jet_clust_seq(inputs, _def);
      return jet_clust_seq.exclusive_jets(n_jets);
   }
   
private:
   fastjet::JetDefinition _def;

};

//------------------------------------------------------------------------
/// \class AxesFinderFromWTA_KT
// This class finds axes by finding the exlusive jets after clustering according to a kT algorithm and a 
// winner take all recombination scheme.
class AxesFinderFromWTA_KT : public AxesFinderFromExclusiveJetDefinition { 

public:
   AxesFinderFromWTA_KT()
   : AxesFinderFromExclusiveJetDefinition(
      fastjet::JetDefinition(fastjet::kt_algorithm,
      fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
      &_recomb,
      fastjet::Best)) {}
   
private:
   const WinnerTakeAllRecombiner _recomb;

};
   
//------------------------------------------------------------------------
/// \class AxesFinderFromWTA_CA
// This class finds axes by finding the exlusive jets after clustering according to a CA algorithm and a 
// winner take all recombination scheme.
class AxesFinderFromWTA_CA : public AxesFinderFromExclusiveJetDefinition {
public:
   AxesFinderFromWTA_CA()
   : AxesFinderFromExclusiveJetDefinition(
      fastjet::JetDefinition(fastjet::cambridge_algorithm, 
      fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
      &_recomb,
      fastjet::Best)) {}
   
private:
   const WinnerTakeAllRecombiner _recomb;
};


//------------------------------------------------------------------------
/// \class AxesFinderFromKT
// This class finds axes by finding the exlusive jets after clustering according to a kT algorithm and a 
// E_scheme recombination.
class AxesFinderFromKT : public AxesFinderFromExclusiveJetDefinition {
public:
   AxesFinderFromKT()
   : AxesFinderFromExclusiveJetDefinition(
      fastjet::JetDefinition(fastjet::kt_algorithm,
      fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
      fastjet::E_scheme,
      fastjet::Best)) {}
};

//------------------------------------------------------------------------
/// \class AxesFinderFromCA
// This class finds axes by finding the exlusive jets after clustering according to a CA algorithm and a 
// E_scheme recombination.
class AxesFinderFromCA : public AxesFinderFromExclusiveJetDefinition { 
public:
   AxesFinderFromCA()
   : AxesFinderFromExclusiveJetDefinition(
      fastjet::JetDefinition(fastjet::cambridge_algorithm,
                             fastjet::JetDefinition::max_allowable_R,  //maximum jet radius constant
      fastjet::E_scheme,
      fastjet::Best)) {}
};


//------------------------------------------------------------------------
/// \class AxesFinderFromHardestJetDefinition
// This class finds axes by clustering the particles and then finding the n hardest inclusive jets. 
// This can be implemented with different jet algorithms.
class AxesFinderFromHardestJetDefinition : public StartingAxesFinder {
public:
   AxesFinderFromHardestJetDefinition(fastjet::JetDefinition def)
   : _def(def) {}
   
   virtual std::vector<fastjet::PseudoJet> getAxes(int n_jets,
                                                   const std::vector <fastjet::PseudoJet> & inputs) const {
      fastjet::ClusterSequence jet_clust_seq(inputs, _def);
      std::vector<fastjet::PseudoJet> myJets = sorted_by_pt(jet_clust_seq.inclusive_jets());
      myJets.resize(n_jets);  // only keep n hardest
      return myJets;
   }
   
private:
   fastjet::JetDefinition _def;
};

//------------------------------------------------------------------------
/// \class AxesFinderFromAntiKT
// This class finds axes by finding the n hardest jets after clustering the particles according 
// to an anti kT algorithm and E_scheme.
class AxesFinderFromAntiKT : public AxesFinderFromHardestJetDefinition {
public:
   AxesFinderFromAntiKT(double R0)
   : AxesFinderFromHardestJetDefinition(
      fastjet::JetDefinition(fastjet::antikt_algorithm,
                             R0,fastjet::E_scheme,fastjet::Best)) {}
};

   

///////
//
// Refining Axes Finders
//
///////

//------------------------------------------------------------------------
/// \class RefiningAxesFinder
// This is the base class for all refining axes finders. These require seed axes, which are
// then refined using an algorithm that is specific to a particular measure.
class RefiningAxesFinder {
   
public:
   
   // This function should be overloaded, and updates the seedAxes to return new axes (should be deterministic)
   virtual std::vector<fastjet::PseudoJet> getOnePassAxes(int n_jets,
                                                          const std::vector<fastjet::PseudoJet>& inputs,
                                                          const std::vector<fastjet::PseudoJet>& seedAxes) const = 0;
   
   // Uses the value of Npass to decide whether to do one pass or multi-pass minimization
   std::vector<fastjet::PseudoJet> getAxes(int n_jets,
                                           const std::vector<fastjet::PseudoJet>& inputs,
                                           const std::vector<fastjet::PseudoJet>& seedAxes) const;
   
   void setMultiPassMode(int n_pass) {
      if (n_pass <= 1) Error("MultiPassMode can only be set with n_pass > 1.");
      _noise_range = 1.0;  //TODO:  Allow this to be changed by the user
      _Npass = n_pass;
   }
   
   void setOnePassMode() { _Npass = 1; }
   
   //virtual destructor
   virtual ~RefiningAxesFinder(){}
   
protected:
   
   // Constructor
   RefiningAxesFinder() : _Npass(1) // default to one pass in case setOnePassMode() isn't called
   {}
   
   // information for multi pass mode
   int _Npass;
   double _noise_range; // noise range for random initialization
   
   // This has to be set in each derived class so multi pass minimization knows when better axes are found
   SharedPtr<MeasureFunction> _associatedMeasure;
   
   // Does multi-pass minimization by jiggling the axes.
   std::vector<fastjet::PseudoJet> getMultiPassAxes(int n_jets,
                                                    const std::vector<fastjet::PseudoJet>& inputs,
                                                    const std::vector<fastjet::PseudoJet>& seedAxes) const;
   
   
   PseudoJet jiggle(const PseudoJet& axis) const;
};


//This is a helper class for the Minimum Axes Finders. It is defined later.
class LightLikeAxis;                                          


//------------------------------------------------------------------------
/// \class AxesFinderFromConicalMinimization
// This class defines an AxesFinder that uses Kmeans minimization, but only on a single pass.
class AxesFinderFromConicalMinimization : public RefiningAxesFinder {

public:

   // From a startingFinder, try to minimize the unnormalized_measure
   AxesFinderFromConicalMinimization(double beta, double Rcutoff)
      : _precision(0.0001), //hard coded for now
        _halt(1000), //hard coded for now
        _beta(beta),
        _Rcutoff(Rcutoff)
   {
      _associatedMeasure.reset(new ConicalUnnormalizedMeasureFunction(beta,Rcutoff));
   }

   virtual std::vector<fastjet::PseudoJet> getOnePassAxes(int n_jets,
                                                          const std::vector <fastjet::PseudoJet> & inputJets,
                                                          const std::vector<fastjet::PseudoJet>& currentAxes) const;
   
private:
   double _precision;  // Desired precision in axes alignment
   int _halt;  // maximum number of steps per iteration
   
   double _beta;
   double _Rcutoff;
   
   // convenient shorthand for squaring
   static inline double sq(double x) {return x*x;}
   
   template <int N> std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes,
                                                              const std::vector <fastjet::PseudoJet> & inputJets) const;
   
   std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes,
                                         const std::vector <fastjet::PseudoJet> & inputJets) const;

};

//------------------------------------------------------------------------
/// \class AxesFinderFromGeometricMinimization
// This class finds axes by minimizing the Lorentz dot product distance between axes and particles. Given a first set of starting axes,
// it essentially does stable cone finxing.
class AxesFinderFromGeometricMinimization : public RefiningAxesFinder {

public:
   AxesFinderFromGeometricMinimization(double beta, double Rcutoff)
   :  _nAttempts(100),
      _accuracy(0.000000001)
   {
      if (beta != 2.0) {
         throw Error("Geometric minimization is currently only defined for beta = 2.0.");
      }
      _associatedMeasure.reset(new GeometricMeasureFunction(beta,Rcutoff));
   }

   virtual std::vector<fastjet::PseudoJet> getOnePassAxes(int n_jets, const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& currentAxes) const;

private:
   double _nAttempts;
   double _accuracy;

};

//------------------------------------------------------------------------
/// \class LightLikeAxis
// This is a helper class for the minimum Axes Finders classes above. It creates a convenient way of defining axes
// in order to better facilitate calculations.
class LightLikeAxis {

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

   // return PseudoJet with information
   fastjet::PseudoJet ConvertToPseudoJet();
   
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

private:
   double _rap, _phi, _weight, _mom;
   
   double DistanceSq(double rap2, double phi2) const {
      double rap1 = _rap;
      double phi1 = _phi;
      
      double distRap = rap1-rap2;
      double distPhi = std::fabs(phi1-phi2);
      if (distPhi > M_PI) {distPhi = 2.0*M_PI - distPhi;}
      return distRap*distRap + distPhi*distPhi;
   }
   
   double Distance(double rap2, double phi2) const {
      return std::sqrt(DistanceSq(rap2,phi2));
   }
   
   
};

} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_AXESFINDER_HH__
