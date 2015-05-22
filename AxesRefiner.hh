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


#ifndef __FASTJET_CONTRIB_AXESREFINER_HH__
#define __FASTJET_CONTRIB_AXESREFINER_HH__

#include "WinnerTakeAllRecombiner.hh"
#include "MeasureDefinition.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

// BIG TODO:  Arrange things so this class is not needed.  (Put NPass in AxesDefinition, associate refining with measure)

///////
//
// AxesRefiner
//
///////

//------------------------------------------------------------------------
/// \class AxesRefiner
// This is the base class for all refining axes finders. These require seed axes, which are
// then refined using an algorithm that is specific to a particular measure.  Like
// AxesDefinition, you have to specify the number of desired axes (even though it is redundantly
// encoded in seedAxes.size()), since otherwise the user could reverse inputs and seeds.
class AxesRefiner {
   
public:
   
   // This function should be overloaded, and updates the seedAxes to return new axes (should be deterministic)
   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets,
                                                             const std::vector<fastjet::PseudoJet>& inputs,
                                                             const std::vector<fastjet::PseudoJet>& seedAxes) const = 0;

   // Uses the value of Npass to decide whether to do one pass or multi-pass minimization
   std::vector<fastjet::PseudoJet> get_axes(int n_jets,
                                            const std::vector<fastjet::PseudoJet>& inputs,
                                            const std::vector<fastjet::PseudoJet>& seedAxes) const;
   
   //virtual destructor
   virtual ~AxesRefiner(){}
   
protected:
   
   // This has to be set in each derived class so multi-pass minimization knows when better axes are found
   SharedPtr<MeasureDefinition> _associatedMeasure;
   
   // Constructor
   AxesRefiner(int nPass) : _Npass(nPass), _noise_range(1.0) //TODO:  Allow noise range to be changed by the user
   {
      if (nPass < 0) throw Error("AxesRefiner needs nPass >= 0.");
   }
   
   // information for multi pass mode
   int _Npass;
   double _noise_range; // noise range for random initialization
   

   
   // Does multi-pass minimization by jiggling the axes.
   std::vector<fastjet::PseudoJet> get_multi_pass_axes(int n_jets,
                                                       const std::vector<fastjet::PseudoJet>& inputs,
                                                       const std::vector<fastjet::PseudoJet>& seedAxes) const;
   
   PseudoJet jiggle(const PseudoJet& axis) const;

};

//This is a helper class for the Minimum Axes Finders. It is defined later.
class LightLikeAxis;                                          


//------------------------------------------------------------------------
/// \class AxesFinderFromConicalMinimization
// This class defines an AxesFinder that uses Kmeans minimization, but only on a single pass.
class ConicalAxesRefiner : public AxesRefiner {

public:

   // From a startingFinder, try to minimize the unnormalized_measure
   ConicalAxesRefiner(double beta, double Rcutoff, int nPass)
      : AxesRefiner(nPass),
        _precision(0.0001), //hard coded for now
        _halt(1000), //hard coded for now
        _beta(beta),
        _Rcutoff(Rcutoff)
   {
      _associatedMeasure.reset(new UnnormalizedCutoffMeasure(beta,Rcutoff));
   }

   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets,
                                                             const std::vector <fastjet::PseudoJet> & inputJets,
                                                             const std::vector<fastjet::PseudoJet>& currentAxes) const;
   
private:
   double _precision;  // Desired precision in axes alignment
   int _halt;  // maximum number of steps per iteration
   
   double _beta;
   double _Rcutoff;
   
   // convenient shorthand for squaring
   static inline double sq(double x) {return x*x;}
   
   // TODO:  consider getting rid of this for readability.
   template <int N> std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes,
                                                              const std::vector <fastjet::PseudoJet> & inputJets) const;
   
   std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes,
                                         const std::vector <fastjet::PseudoJet> & inputJets) const;

};

//------------------------------------------------------------------------
/// \class AxesFinderFromGeometricMinimization
// This class finds axes by minimizing the Lorentz dot product distance between axes and particles. Given a first set of starting axes,
// it essentially does stable cone finxing.
class GeometricAxesRefiner : public AxesRefiner {

public:
   GeometricAxesRefiner(double beta, double Rcutoff, int nPass)
   :  AxesRefiner(nPass),
      _nAttempts(100),
      _accuracy(0.000000001)
   {
      if (beta != 2.0) {
         throw Error("Geometric minimization is currently only defined for beta = 2.0.");
      }
      _associatedMeasure.reset(new OriginalGeometricMeasure(Rcutoff));
   }

   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets, const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& currentAxes) const;

private:
   double _nAttempts;
   double _accuracy;

};

// added by TJW

//------------------------------------------------------------------------
/// \class GeneralAxesFinder
// This class finds axes by minimizing the Lorentz dot product distance between axes and particles. Given a first set of starting axes,
// it essentially does stable cone finxing.
class GeneralAxesRefiner : public AxesRefiner {

public:
  // GeneralAxesRefiner(double beta, double Rcutoff, int nPass) 
  GeneralAxesRefiner(MeasureDefinition* associatedMeasure, int nPass) 
  : AxesRefiner(nPass),
    _nAttempts(100),
    _accuracy(0.00001) 
  {
    _associatedMeasure.reset(associatedMeasure);
  }

  virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets, const std::vector<fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& currentAxes) const;

private:
    double _nAttempts;
    double _accuracy;

   // create light-like axis -- added to easily make light-like axes, may not be necessary
   fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
      double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
      return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
   }
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
