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

#ifndef __FASTJET_CONTRIB_MEASUREDEFINITION_HH__
#define __FASTJET_CONTRIB_MEASUREDEFINITION_HH__

#include "fastjet/PseudoJet.hh"
#include <cmath>
#include <vector>
#include <list>
#include <limits>

#include "TauComponents.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{



// The following Measures are available (and their relevant arguments):
class DefaultMeasure;               // Default measure from which next classes derive (should not be called directly)
class NormalizedMeasure;            // (beta,R0)
class UnnormalizedMeasure;          // (beta)
class NormalizedCutoffMeasure;      // (beta,R0,Rcutoff)
class UnnormalizedCutoffMeasure;    // (beta,Rcutoff)

// Formerly GeometricMeasure, now no longer recommended, kept only for cross-check purposes
//class DeprecatedGeometricMeasure;         // (beta)
//class DeprecatedGeometricCutoffMeasure;   // (beta,Rcutoff)

// New measures as of v2.2
class ConicalMeasure;               // (beta,Rcutoff)
class OriginalGeometricMeasure;     // (Rcutoff)
class ModifiedGeometricMeasure;     // (Rcutoff)
class ConicalGeometricMeasure;      // (beta, gamma, Rcutoff)
class XConeMeasure;                 // (beta, Rcutoff)

   
///////
//
// MeasureDefinition
//
///////

//This is a helper class for the Minimum Axes Finders. It is defined later.
class LightLikeAxis;                                          
   
//------------------------------------------------------------------------
/// \class MeasureDefinition
// This is the base class for measure definitions.  Derived classes will calculate
// the tau_N of a jet given a specific measure and a set of axes.  The measure is
// determined by various jet and beam distances (and possible normalization factors).
class MeasureDefinition {
   
public:
   
   // Description of measure and parameters
   virtual std::string description() const = 0;
   
   // In derived classes, this should return a copy of the corresponding derived class
   virtual MeasureDefinition* create() const = 0;

   //The following five functions define the measure by which tau_N is calculated,
   //and are overloaded by the various measures below
   
   // Distanes to axes.  These are called many times, so need to be as fast as possible
   // Unless overloaded, they just call jet_numerator and beam_numerator
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return jet_numerator(particle,axis);
   }
   virtual double beam_distance_squared(const fastjet::PseudoJet& particle) const {
      return beam_numerator(particle);
   }
   
   // The actual measures used in N-(sub)jettiness
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const = 0;
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const = 0;
   
   // a possible normalization factor
   virtual double denominator(const fastjet::PseudoJet& particle) const = 0;
   
public:
   
   // Returns the tau value for a give set of particles and axes
   double result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const {
      return component_result(particles,axes).tau();
   }
   
   // Short-hand for the above function.
   inline double operator() (const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const {
      return result(particles,axes);
   }
   
   // Return all of the TauComponents for specific input particles and axes
   TauComponents component_result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const;
   
   // Create the partitioning according the jet/beam distances and stores them a TauPartition
   TauPartition get_partition(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const;

   // calculates the tau result using an existing partition
   TauComponents component_result_from_partition(const TauPartition& partition, const std::vector<fastjet::PseudoJet>& axes) const;

   
   // Create associated axes refiner (i.e. one-pass minimization routine), if available, otherwise use default
   
   // added this to retrieve one-pass axes based on measure-specific minimization scheme
   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets,
                                                             const std::vector<fastjet::PseudoJet>& inputs,
                                                             const std::vector<fastjet::PseudoJet>& seedAxes,
                                                             int nAttempts = 1000,         // cap number of iterations
                                                             double accuracy = 0.0001      // cap distance of closest approach
                                                             ) const;

   //virtual destructor
   virtual ~MeasureDefinition(){}
   
protected:
   
   //flag set by derived classes to choose whether or not to use beam/denominator
   TauMode _tau_mode;
   bool _useAxisScaling;

   // This constructor allows _tau_mode to be set by derived classes
   MeasureDefinition() : _tau_mode(UNDEFINED_SHAPE) {}

   void setTauMode(TauMode tau_mode) {
      _tau_mode = tau_mode;
   }

   void setAxisScaling(bool useAxisScaling) {
      _useAxisScaling = useAxisScaling;
   }

   // methods added to generalize energy and angle squared for different measure types
   bool has_denominator() const { return (_tau_mode == NORMALIZED_JET_SHAPE || _tau_mode == NORMALIZED_EVENT_SHAPE);}
   bool has_beam() const {return (_tau_mode == UNNORMALIZED_EVENT_SHAPE || _tau_mode == NORMALIZED_EVENT_SHAPE);}

   // create light-like axis (used in default one pass minimization
   fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
      double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
      return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
   }
   
   // shorthand for squaring
   static inline double sq(double x) {return x*x;}

};
   

///////
//
// Default Measures
//
///////
   
   
// For switching between pp and ee measure types in Default Measure
enum DefaultMeasureType {
   pt_R,       ///  use transverse momenta and boost-invariant angles,
   E_theta,    ///  use energies and angles,
   lorentz_dot, ///  use dot product inspired measure
   perp_lorentz_dot /// use conical geometric inspired measures
};

//------------------------------------------------------------------------
/// \class DefaultMeasure
// This class is the default measure as defined in the original N-subjettiness papers.
// Based on the conical measure, but with a normalization factor
// This measure is defined as the pT of the particle multiplied by deltaR
// to the power of beta. This class includes the normalization factor determined by R0
class DefaultMeasure : public MeasureDefinition {
   
public:

   virtual std::string description() const;
   
   virtual DefaultMeasure* create() const {return new DefaultMeasure(*this);}

   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return angleSquared(particle, axis);
   }
   
   virtual double beam_distance_squared(const fastjet::PseudoJet& /*particle*/) const {
      return _RcutoffSq;
   }
   
   // general definitions of energy and angle
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const{
      double jet_dist = angleSquared(particle, axis);
      if (jet_dist > 0.0) {
         return energy(particle) * std::pow(jet_dist,_beta/2.0);
      } else {
         return 0.0;
      }
   }
   
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      return energy(particle) * std::pow(_Rcutoff,_beta);
   }
   
   virtual double denominator(const fastjet::PseudoJet& particle) const {
      return energy(particle) * std::pow(_R0,_beta);
   }

   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets,
                                                             const std::vector<fastjet::PseudoJet>& inputs,
                                                             const std::vector<fastjet::PseudoJet>& seedAxes,
                                                             int nAttempts,   // cap number of iterations
                                                             double accuracy  // cap distance of closest approach
                                                             ) const;
   
protected:
   double _beta;
   double _R0;
   double _Rcutoff;
   double _RcutoffSq;   
   DefaultMeasureType _measure_type;
   
   
   // Constructor is protected so that no one tries to call this directly.
   DefaultMeasure(double beta, double R0, double Rcutoff, DefaultMeasureType measure_type = pt_R)
   : MeasureDefinition(), _beta(beta), _R0(R0), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)), _measure_type(measure_type)
   {
      setAxisScaling(true); // for those measure types that will require the default minimization
   }
   
   
   // added set measure method in case it becomes useful later
   void setDefaultMeasureType(DefaultMeasureType measure_type) {
      _measure_type = measure_type;
   }

   double energy(const PseudoJet& jet) const;
   double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;

   // added so that description will include the measure type
   std::string measure_type_name() const {
      if (_measure_type == pt_R) return "pt_R";
      else if (_measure_type == E_theta) return "E_theta";
      else if (_measure_type == lorentz_dot) return "lorentz_dot";
      else if (_measure_type == perp_lorentz_dot) return "perp_lorentz_dot";
      else return "Measure Type Undefined";
   }      

   template <int N> std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes,
                                                              const std::vector <fastjet::PseudoJet> & inputJets,
                                                              double accuracy) const;
   
   std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes,
                                         const std::vector <fastjet::PseudoJet> & inputJets,
                                         double accuracy) const;
};
   

//------------------------------------------------------------------------
/// \class NormalizedCutoffMeasure
// This measure is just a wrapper for DefaultMeasure
class NormalizedCutoffMeasure : public DefaultMeasure {

public:

   NormalizedCutoffMeasure(double beta, double R0, double Rcutoff, DefaultMeasureType measure_type = pt_R) 
   : DefaultMeasure(beta, R0, Rcutoff, measure_type) {
      setTauMode(NORMALIZED_JET_SHAPE);
   }

   virtual std::string description() const;

   virtual NormalizedCutoffMeasure* create() const {return new NormalizedCutoffMeasure(*this);}

};

//------------------------------------------------------------------------
/// \class NormalizedMeasure
// This measure is the same as NormalizedCutoffMeasure, with Rcutoff taken to infinity.
class NormalizedMeasure : public NormalizedCutoffMeasure {

public:

   NormalizedMeasure(double beta, double R0, DefaultMeasureType measure_type = pt_R)
   : NormalizedCutoffMeasure(beta, R0, std::numeric_limits<double>::max(), measure_type) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(NORMALIZED_JET_SHAPE);
   }
   
   virtual std::string description() const;

   virtual NormalizedMeasure* create() const {return new NormalizedMeasure(*this);}

};
   

//------------------------------------------------------------------------
/// \class UnnormalizedCutoffMeasure
// This class is the unnormalized conical measure. The only difference from NormalizedCutoffMeasure
// is that the denominator is defined to be 1.0 by setting _has_denominator to false.
// class UnnormalizedCutoffMeasure : public NormalizedCutoffMeasure {
class UnnormalizedCutoffMeasure : public DefaultMeasure {
   
public:
   // Since all methods are identical, UnnormalizedMeasure inherits directly
   // from NormalizedMeasure. R0 is a dummy value since the value of R0 is unecessary for this class,
   // and the "false" flag sets _has_denominator in MeasureDefinition to false so no denominator is used.
   UnnormalizedCutoffMeasure(double beta, double Rcutoff, DefaultMeasureType measure_type = pt_R)
   : DefaultMeasure(beta, std::numeric_limits<double>::quiet_NaN(), Rcutoff, measure_type) {
      setTauMode(UNNORMALIZED_EVENT_SHAPE);
   }

   virtual std::string description() const;
   
   virtual UnnormalizedCutoffMeasure* create() const {return new UnnormalizedCutoffMeasure(*this);}

};

   
//------------------------------------------------------------------------
/// \class UnnormalizedMeasure
// This measure is the same as UnnormalizedCutoffMeasure, with Rcutoff taken to infinity.
class UnnormalizedMeasure : public UnnormalizedCutoffMeasure {
   
public:
   // Since all methods are identical, UnnormalizedMeasure inherits directly
   // from NormalizedMeasure. R0 is a dummy value since the value of R0 is unecessary for this class,
   // and the "false" flag sets _has_denominator in MeasureDefinition to false so no denominator is used.
   UnnormalizedMeasure(double beta, DefaultMeasureType measure_type = pt_R)
   : UnnormalizedCutoffMeasure(beta, std::numeric_limits<double>::max(), measure_type) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(UNNORMALIZED_JET_SHAPE);
   }


   virtual std::string description() const;

   virtual UnnormalizedMeasure* create() const {return new UnnormalizedMeasure(*this);}
   
};

////------------------------------------------------------------------------
///// \class DeprecatedGeometricCutoffMeasure
//// This class is the old, incorrectly coded geometric measure.
//// It is kept in case anyone wants to check old code, but should not be used for production purposes.
//class DeprecatedGeometricCutoffMeasure : public MeasureDefinition {
//
//public:
//
//   // Please, please don't use this.
//   DeprecatedGeometricCutoffMeasure(double jet_beta, double Rcutoff)
//   :   MeasureDefinition(),
//      _jet_beta(jet_beta),
//      _beam_beta(1.0), // This is hard coded, since alternative beta_beam values were never checked.
//      _Rcutoff(Rcutoff),
//      _RcutoffSq(sq(Rcutoff)) {
//         setTauMode(UNNORMALIZED_EVENT_SHAPE);
//         setAxisScaling(false);
//         if (jet_beta != 2.0) {
//         throw Error("Geometric minimization is currently only defined for beta = 2.0.");
//      }      
//   }
//
//   virtual std::string description() const;
//   
//   virtual DeprecatedGeometricCutoffMeasure* create() const {return new DeprecatedGeometricCutoffMeasure(*this);}
//   
//   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
//      fastjet::PseudoJet lightAxis = lightFrom(axis);
//      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
//      return pseudoRsquared;
//   }
//
//   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
//      return _RcutoffSq;
//   }
//
//   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
//      fastjet::PseudoJet lightAxis = lightFrom(axis);
//      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(lightAxis.pt(),_beam_beta - 1.0);
//      return particle.pt() * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
//   }
//
//   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
//      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(particle.pt()/particle.e(),_beam_beta - 1.0);
//      return particle.pt() * weight * std::pow(_Rcutoff,_jet_beta);
//   }
//
//   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
//      return std::numeric_limits<double>::quiet_NaN();
//   }
//
//   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets,
//                                                             const std::vector<fastjet::PseudoJet>& inputs,
//                                                             const std::vector<fastjet::PseudoJet>& seedAxes,
//                                                             int nAttempts,    // cap number of iterations
//                                                             double accuracy   // cap distance of closest approach
//                                                            ) const;
//   
//protected:
//   double _jet_beta;
//   double _beam_beta;
//   double _Rcutoff;
//   double _RcutoffSq;
//
//};
//
//// ------------------------------------------------------------------------
//// / \class DeprecatedGeometricMeasure
//// Same as DeprecatedGeometricMeasureCutoffMeasure, but with Rcutoff taken to infinity.
//// NOTE:  This class should not be used for production purposes.
//class DeprecatedGeometricMeasure : public DeprecatedGeometricCutoffMeasure {
//   
//public:
//   DeprecatedGeometricMeasure(double beta)
//   : DeprecatedGeometricCutoffMeasure(beta,std::numeric_limits<double>::max()) {
//      _RcutoffSq = std::numeric_limits<double>::max();
//      setTauMode(UNNORMALIZED_JET_SHAPE);
//   }
//
//   virtual std::string description() const;
//
//   virtual DeprecatedGeometricMeasure* create() const {return new DeprecatedGeometricMeasure(*this);}
//};


//------------------------------------------------------------------------
/// \class ConicalMeasure
// Same as UnnormalizedCutoffMeasure, but using the new default one-pass minimization algorithm.
// Axes are also made to be light-like to ensure sensible behavior
// Intended to be used as an event shape.
class ConicalMeasure : public MeasureDefinition {
   
public:
   ConicalMeasure(double beta, double Rcutoff)
   :   MeasureDefinition(), _beta(beta), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
      setTauMode(UNNORMALIZED_EVENT_SHAPE);
      setAxisScaling(true);
   }
   
   virtual std::string description() const;
   
   virtual ConicalMeasure* create() const {return new ConicalMeasure(*this);}
   
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      PseudoJet lightAxis = lightFrom(axis);
      return particle.squared_distance(lightAxis);
   }
   
   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      return _RcutoffSq;
   }
   
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      PseudoJet lightAxis = lightFrom(axis);
      double jet_dist = particle.squared_distance(lightAxis)/_RcutoffSq;
      double jet_perp = particle.perp();
      
      if (_beta == 2.0) {
         return jet_perp * jet_dist;
      } else {
         return jet_perp * pow(jet_dist,_beta/2.0);
      }
   }
   
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      return particle.perp();
   }
   
   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
protected:
   double _beta;
   double _Rcutoff;
   double _RcutoffSq;
};
   


//------------------------------------------------------------------------
/// \class OriginalGeometricMeasure
// This class is the original (and hopefully now correctly coded) geometric measure.
// This measure is defined by the Lorentz dot product between
// the particle and the axis.  This class does not include normalization of tau_N.
// New in Nsubjettiness v2.2
// NOTE: This is defined differently from the DeprecatedGeometricMeasure above.
class OriginalGeometricMeasure : public MeasureDefinition {

public:
   OriginalGeometricMeasure(double Rcutoff)
   :   MeasureDefinition(),
     _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
         setAxisScaling(false);
      }

   virtual std::string description() const;
   
   virtual OriginalGeometricMeasure* create() const {return new OriginalGeometricMeasure(*this);}
   
   // use default jet_distance_squared and beam_distance_squared

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return dot_product(lightFrom(axis), particle)/_RcutoffSq;
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet beam_a(1,0,0,1);
      fastjet::PseudoJet beam_b(-1,0,0,1);
      double min_perp = std::min(dot_product(beam_a, particle),dot_product(beam_b, particle));
      return min_perp;
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
protected:
   double _Rcutoff;
   double _RcutoffSq;

};


//------------------------------------------------------------------------
/// \class ModifiedGeometricMeasure
// This class is the Modified geometric measure.  This jet measure is defined by the Lorentz dot product between
// the particle and the axis, as in the Original Geometric Measure. The beam measure is defined differently from 
// the above OriginalGeometric to allow for more conical jets. New in Nsubjettiness v2.2
class ModifiedGeometricMeasure : public MeasureDefinition {

public:
   ModifiedGeometricMeasure(double Rcutoff)
   :   MeasureDefinition(),
     _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
         setAxisScaling(false);
      }

   virtual std::string description() const;
   
   virtual ModifiedGeometricMeasure* create() const {return new ModifiedGeometricMeasure(*this);}
   
   // use default jet_distance_squared and beam_distance_squared
   
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return dot_product(lightFrom(axis), particle)/_RcutoffSq;
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet lightParticle = lightFrom(particle);
      return 0.5*particle.mperp()*lightParticle.pt();
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
protected:
   double _Rcutoff;
   double _RcutoffSq;

};

// ------------------------------------------------------------------------
// / \class ConicalGeometricMeasure
// This class is the Conical Geometric measure.  This measure is defined by the Lorentz dot product between
// the particle and the axis normalized by the axis and particle pT, as well as a factor of cosh(y) to vary
// the rapidity depepdence of the beam. New in Nsubjettiness v2.2, and the basis for the XCone jet algorithm
class ConicalGeometricMeasure : public MeasureDefinition {

public:
   ConicalGeometricMeasure(double jet_beta, double beam_gamma, double Rcutoff)
   :   MeasureDefinition(),
      _jet_beta(jet_beta), _beam_gamma(beam_gamma), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)){
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
         setAxisScaling(true);
      }

   virtual std::string description() const;
   
   virtual ConicalGeometricMeasure* create() const {return new ConicalGeometricMeasure(*this);}
   
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      return pseudoRsquared;
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      return _RcutoffSq;
   }

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      double jet_dist = jet_distance_squared(particle,axis)/_RcutoffSq;
      if (jet_dist > 0.0) {
         fastjet::PseudoJet lightParticle = lightFrom(particle);
         double weight = (_beam_gamma == 1.0) ? 1.0 : std::pow(0.5*lightParticle.pt(),_beam_gamma - 1.0);
         return particle.pt() * weight * std::pow(jet_dist,_jet_beta/2.0);
      } else {
         return 0.0;
      }
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet lightParticle = lightFrom(particle);
      double weight = (_beam_gamma == 1.0) ? 1.0 : std::pow(0.5*lightParticle.pt(),_beam_gamma - 1.0);
      return particle.pt() * weight;
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
      
protected:
   double _jet_beta;
   double _beam_gamma;
   double _Rcutoff;
   double _RcutoffSq;
   
};

// ------------------------------------------------------------------------
// / \class XConeCutoffMeasure
// This class is the XCone Measure.  This is the default measure for use with the
// XCone algorithm. It is identical to the conical geometric measure but with gamma = 1.0.
class XConeMeasure : public ConicalGeometricMeasure {

public:
   XConeMeasure(double jet_beta, double R)
   :   ConicalGeometricMeasure(jet_beta,
                               1.0, // beam_gamma, hard coded
                               R    // Rcutoff scale
                               ) { }

   virtual std::string description() const;

   virtual XConeMeasure* create() const {return new XConeMeasure(*this);}

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

#endif  // __FASTJET_CONTRIB_MEASUREDEFINITION_HH__
