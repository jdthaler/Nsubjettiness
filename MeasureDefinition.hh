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

// for switching between pp and ee measure types -- TJW
enum MeasureType {
   pt_R,       ///  use transverse momenta and boost-invariant angles,
   E_theta,    ///  use energies and angles,
   lorentz_dot ///  use dot product inspired measure -- do we need this? -- TJW
};
 

// The following Measures are available (and their relevant arguments):
class DefaultMeasure;            // Default Conical measure from which all conical classes derive
class NormalizedMeasure;         // (beta,R0)
class UnnormalizedMeasure;       // (beta)
class DeprecatedGeometricMeasure;          // (beta) -- gone TJW
class NormalizedCutoffMeasure;   // (beta,R0,Rcutoff)
class UnnormalizedCutoffMeasure; // (beta,Rcutoff)
class DeprecatedGeometricCutoffMeasure;    // (beta,Rcutoff) -- gone TJW

class OriginalGeometricMeasure; // TESTING -- TJW
class ModifiedGeometricMeasure; // TESTING -- TJW
class ConicalGeometricCutoffMeasure; // TESTING -- TJW
class ConicalGeometricMeasure; // TESTING -- TJW
class XConeCutoffMeasure; // TESTING -- TJW
class XConeMeasure; // TESTING -- TJW
///////
//
// MeasureDefinition
//
///////

class AxesRefiner;  // Needed to avoid clashes with .hh files

   
//------------------------------------------------------------------------
/// \class MeasureDefinition
// This is the base class for measure definitions.  Derived classes will calculates
// the tau_N of a jet given a specific measure and a set of axes.  The measure is
// determined by various jet and beam distances (and possible normalization factors).
class MeasureDefinition {
   
public:
   
   // Description of measure and parameters
   virtual std::string description() const = 0;
   
   // In derived classes, this should return a copy of the corresponding derived class
   virtual MeasureDefinition* create() const = 0;

   //The follwoing five functions define the measure by which tau_N is calculated,
   //and are overloaded by the various measures below
   
   // Distanes to axes.  These are called many times, so need to be as fast as possible
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const = 0;
   virtual double beam_distance_squared(const fastjet::PseudoJet& particle) const = 0;

   // axes_numerator added for use in minimization in AxesRefiner -- TJW
   virtual double axes_numerator(const fastjet::PseudoJet& /*particle*/, const fastjet::PseudoJet& /*axis*/) const {
      return 1.0;
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

   // Create associated axes refiner (i.e. minimization routine), if available.
   // virtual AxesRefiner* createAxesRefiner(int /*nPass*/) const = 0;
   virtual AxesRefiner* createAxesRefiner(int /*nPass*/) const;
   AxesRefiner* createAxesRefiner() const { return createAxesRefiner(1);} // no argument means one pass
   
   // shorthand for squaring
   // TODO:  This should be moved to special function file?
   static inline double sq(double x) {return x*x;}

   //virtual destructor
   virtual ~MeasureDefinition(){}
   
protected:
   
   //flag set by derived classes to choose whether or not to use beam/denominator
   TauMode _tau_mode;

   // removed old constructor so that users cannot change the value of TauMode -- TJW
   // added constructor with MeasureType argument so user will be able to set measure type when they need to -- TJW

   // This constructor allows _has_denominator to be set by derived classes
   // MeasureDefinition(TauMode tau_mode) : _tau_mode(tau_mode) {}
   MeasureDefinition() : _tau_mode(UNDEFINED_SHAPE) {}

   void setTauMode(TauMode tau_mode) {
      _tau_mode = tau_mode;
   }

   // new methods added to generalize energy and angle squared for different measure types -- TJW

   bool has_denominator() const { return (_tau_mode == NORMALIZED_JET_SHAPE || _tau_mode == NORMALIZED_EVENT_SHAPE);}
   bool has_beam() const {return (_tau_mode == UNNORMALIZED_EVENT_SHAPE || _tau_mode == NORMALIZED_EVENT_SHAPE);}

};


// updated all classes to remove TauMode argument in constructor, added default MeasureType argument to constructor -- TJW 
// flipped inheritance of NormalizedCutoffMeasure and NormalizedMeasure to avoid issue of squaring a maximum value -- TJW
// updated instances of perp() and jet_distance_squared() to account for more general metrics -- TJW

// Added by TJW
//------------------------------------------------------------------------
/// \class DefaultMeasure
// This class is the default measure, based on the conical measure.
// This measure is defined as the pT of the particle multiplied by deltaR
// to the power of beta. This class includes the normalization factor determined by R0
class DefaultMeasure : public MeasureDefinition {
   
public:
   
   DefaultMeasure(double beta, double R0, double Rcutoff, MeasureType measure_type = pt_R)
   : MeasureDefinition(), _beta(beta), _R0(R0), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)), _measure_type(measure_type) {}

   virtual std::string description() const;
   
   virtual DefaultMeasure* create() const {return new DefaultMeasure(*this);}

   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return particle.squared_distance(axis);
   }
   
   virtual double beam_distance_squared(const fastjet::PseudoJet& /*particle*/) const {
      // return sq(_Rcutoff);
      return _RcutoffSq;
   }
   
   // updated for new general definitions of energy and angle -- TJW
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const{
      // return particle.perp() * std::pow(jet_distance_squared(particle,axis),_beta/2.0);
      return energy(particle) * std::pow(angleSquared(particle, axis), _beta/2.0);
   }
   
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      // return particle.perp() * std::pow(_Rcutoff,_beta);
      return energy(particle) * std::pow(_Rcutoff,_beta);
   }
   
   virtual double denominator(const fastjet::PseudoJet& particle) const {
      // return particle.perp() * std::pow(_R0,_beta);
      return energy(particle) * std::pow(_R0,_beta);
   }

   // The minimization routine is ConicalAxesRefiner
   AxesRefiner* createAxesRefiner(int nPass) const;
   
protected:
   double _beta;
   double _R0;
   double _Rcutoff;
   double _RcutoffSq;   
   MeasureType _measure_type;

   // added set measure method in case it becomes useful later -- TJW
   void setMeasureType(MeasureType measure_type) {
      _measure_type = measure_type;
   }

   double energy(const PseudoJet& jet) const;
   double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;

   // added so that description will include the measure type -- TJW
   std::string measure_type_name() const {
      if (_measure_type == pt_R) return "pt_R";
      else if (_measure_type == E_theta) return "E_theta";
      else if (_measure_type == lorentz_dot) return "lorentz_dot";
      else return "Measure Type Undefined";
   }      

   
};
   
class NormalizedCutoffMeasure : public DefaultMeasure {

public:

   NormalizedCutoffMeasure(double beta, double R0, double Rcutoff, MeasureType measure_type = pt_R) 
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

   NormalizedMeasure(double beta, double R0, MeasureType measure_type = pt_R)
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
   UnnormalizedCutoffMeasure(double beta, double Rcutoff, MeasureType measure_type = pt_R)
   : DefaultMeasure(beta, std::numeric_limits<double>::quiet_NaN(), Rcutoff, measure_type) {
      setTauMode(UNNORMALIZED_EVENT_SHAPE);
   }

   virtual std::string description() const;
   
   virtual UnnormalizedCutoffMeasure* create() const {return new UnnormalizedCutoffMeasure(*this);}

};

   
//------------------------------------------------------------------------
/// \class UnnormalizedMeasure
// This measure is the same as UnnormalizedCutoffMeasure, with Rcutoff taken
// to infinity.
class UnnormalizedMeasure : public UnnormalizedCutoffMeasure {
   
public:
   // Since all methods are identical, UnnormalizedMeasure inherits directly
   // from NormalizedMeasure. R0 is a dummy value since the value of R0 is unecessary for this class,
   // and the "false" flag sets _has_denominator in MeasureDefinition to false so no denominator is used.
   UnnormalizedMeasure(double beta, MeasureType measure_type = pt_R)
   : UnnormalizedCutoffMeasure(beta, std::numeric_limits<double>::max(), measure_type) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(UNNORMALIZED_JET_SHAPE);
   }


   virtual std::string description() const;

   virtual UnnormalizedMeasure* create() const {return new UnnormalizedMeasure(*this);}
   
};

//------------------------------------------------------------------------
/// \class GeometricCutoffMeasure
// This class is the geometric measure.  This measure is defined by the Lorentz dot product between
// the particle and the axis.  This class does not include normalization of tau_N.
// NOTE:  This class should not be used for production purposes.
class DeprecatedGeometricCutoffMeasure : public MeasureDefinition {

public:
   // Right now, we are hard coded for beam_beta = 1.0, but that will need to change
   DeprecatedGeometricCutoffMeasure(double jet_beta, double Rcutoff)
   :   MeasureDefinition(),
      _jet_beta(jet_beta), _beam_beta(1.0), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
      }

   virtual std::string description() const;
   
   virtual DeprecatedGeometricCutoffMeasure* create() const {return new DeprecatedGeometricCutoffMeasure(*this);}
   
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      return pseudoRsquared;
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      // return sq(_Rcutoff);
      return _RcutoffSq;
   }

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(lightAxis.pt(),_beam_beta - 1.0);
      return particle.pt() * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(particle.pt()/particle.e(),_beam_beta - 1.0);
      return particle.pt() * weight * std::pow(_Rcutoff,_jet_beta);
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
   // The minimization routine is GeometricAxesRefiner
   AxesRefiner* createAxesRefiner(int nPass) const;
   
protected:
   double _jet_beta;
   double _beam_beta;
   double _Rcutoff;
   double _RcutoffSq;
   
   // create light-like axis
   fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
      double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
      return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
   }

};

// ------------------------------------------------------------------------
// / \class GeometricMeasure
// Same as GeometricCutoffMeasure, but with Rcutoff taken to infinity.
// NOTE:  This class should not be used for production purposes.
class DeprecatedGeometricMeasure : public DeprecatedGeometricCutoffMeasure {
   
public:
   DeprecatedGeometricMeasure(double beta)
   : DeprecatedGeometricCutoffMeasure(beta,std::numeric_limits<double>::max()) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(UNNORMALIZED_JET_SHAPE);
   }

   virtual std::string description() const;

   virtual DeprecatedGeometricMeasure* create() const {return new DeprecatedGeometricMeasure(*this);}
};

// Classes below added by TJW

//------------------------------------------------------------------------
/// \class GeometricCutoffMeasure
// This class is the geometric measure.  This measure is defined by the Lorentz dot product between
// the particle and the axis.  This class does not include normalization of tau_N.
// New in Nsubjettiness version (latest)
// NOTE: This is defined differently from the DeprecatedGeometricMeasure above. This is the correct 
// Geometric measure and should be used for all production purposes.
// NOTE:  This class is in flux and should not be used for production purposes.
class OriginalGeometricMeasure : public MeasureDefinition {

public:
   OriginalGeometricMeasure(double Rcutoff)
   :   MeasureDefinition(),
     _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
      }

   virtual std::string description() const;
   
   virtual OriginalGeometricMeasure* create() const {return new OriginalGeometricMeasure(*this);}
   
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      return pseudoRsquared;
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      // return sq(_Rcutoff);
      return _RcutoffSq;
   }

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return dot_product(lightFrom(axis), particle);
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet beam_a(1,0,0,1);
      fastjet::PseudoJet beam_b(-1,0,0,1);
      double min_perp = std::min(dot_product(beam_a, particle),dot_product(beam_b, particle));
      return _RcutoffSq*min_perp;
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
      
protected:
   double _Rcutoff;
   double _RcutoffSq;
   
   // create light-like axis
   fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
      double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
      return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
   }

};


//------------------------------------------------------------------------
/// \class ModifiedGeometricMeasure
// This class is the Modified geometric measure.  This jet measure is defined by the Lorentz dot product between
// the particle and the axis, as in the Original Geometric Measure. The beam measure is defined differently from 
// the above OriginalGeometric to allow for more conical jets. New in Nsubjettiness version (latest)
// NOTE:  This class is in flux and should not be used for production purposes.
class ModifiedGeometricMeasure : public MeasureDefinition {

public:
   ModifiedGeometricMeasure(double Rcutoff)
   :   MeasureDefinition(),
     _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
      }

   virtual std::string description() const;
   
   virtual ModifiedGeometricMeasure* create() const {return new ModifiedGeometricMeasure(*this);}
   
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      return pseudoRsquared;
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      // return sq(_Rcutoff);
      return _RcutoffSq;
   }

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return dot_product(lightFrom(axis), particle);
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet lightParticle = lightFrom(particle);
      return 0.5*particle.mperp()*_RcutoffSq*lightParticle.pt();
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
      
protected:
   double _Rcutoff;
   double _RcutoffSq;
   
   // create light-like axis
   fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
      double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
      return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
   }

};

// ------------------------------------------------------------------------
// / \class ConicalGeometricCutoffMeasure
// This class is the Conical geometric measure.  This measure is defined by the Lorentz dot product between
// the particle and the axis normalized by the axis and particle pT, as well as a factor of cosh(y) to vary
// the rapidity depepdence of the beam. New in Nsubjettiness version (latest)
// NOTE:  This class is in flux and should not be used for production purposes.
class ConicalGeometricMeasure : public MeasureDefinition {

public:
   ConicalGeometricMeasure(double jet_beta, double jet_gamma, double Rcutoff)
   :   MeasureDefinition(),
      _jet_beta(jet_beta), _jet_gamma(jet_gamma), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
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
      fastjet::PseudoJet lightParticle = lightFrom(particle);
      double weight = (_jet_gamma == 1.0) ? 1.0 : std::pow(0.5*lightParticle.pt(),_jet_gamma - 1.0);
      return particle.pt() * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet lightParticle = lightFrom(particle);
      double weight = (_jet_gamma == 1.0) ? 1.0 : std::pow(0.5*lightParticle.pt(),_jet_gamma - 1.0);
      return particle.pt() * weight * std::pow(_Rcutoff,_jet_beta);
   }

   // defined according to G(n,p) in the paper -- TJW
   virtual double axes_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      double pseudoMomentum = dot_product(lightFrom(axis),particle) + 0.0000001*particle.E(); 
      return (double)jet_numerator(particle, axis)/pseudoMomentum;
   }


   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
      
protected:
   double _jet_beta;
   double _jet_gamma;
   double _Rcutoff;
   double _RcutoffSq;
   
   // create light-like axis
   fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
      double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
      return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
   }

};

// ------------------------------------------------------------------------
// / \class XConeCutoffMeasure
// This class is the XCone Measure.  This is the default measure for use with the
// XCone algorithm. It is identical to the conical geomtric measure but with gamma = 1.0.
// NOTE:  This class is in flux and should not be used for production purposes.
class XConeMeasure : public ConicalGeometricMeasure {

public:
   XConeMeasure(double jet_beta, double Rcutoff)
   :   ConicalGeometricMeasure(jet_beta, 1.0, Rcutoff) { }

   virtual std::string description() const;

   virtual XConeMeasure* create() const {return new XConeMeasure(*this);}

};

} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_MEASUREDEFINITION_HH__
