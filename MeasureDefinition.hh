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
class NormalizedMeasure;         // (beta,R0)
class UnnormalizedMeasure;       // (beta)
class GeometricMeasure;          // (beta)
class NormalizedCutoffMeasure;   // (beta,R0,Rcutoff)
class UnnormalizedCutoffMeasure; // (beta,Rcutoff)
class GeometricCutoffMeasure;    // (beta,Rcutoff)

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
   virtual double axes_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const = 0;
   
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
   virtual AxesRefiner* createAxesRefiner(int /*nPass*/) const {return NULL;}
   AxesRefiner* createAxesRefiner() const { return createAxesRefiner(1);} // no argument means one pass
   
   // shorthand for squaring
   // TODO:  This should be moved to special function file?
   static inline double sq(double x) {return x*x;}

   //virtual destructor
   virtual ~MeasureDefinition(){}
   
protected:
   
   //flag set by derived classes to choose whether or not to use beam/denominator
   TauMode _tau_mode;
   MeasureType _measure_type;

   // removed old constructor so that users cannot change the value of TauMode -- TJW
   // added constructor with MeasureType argument so user will be able to set measure type when they need to -- TJW

   // This constructor allows _has_denominator to be set by derived classes
   // MeasureDefinition(TauMode tau_mode) : _tau_mode(tau_mode) {}
   MeasureDefinition() : _tau_mode(UNDEFINED_SHAPE), _measure_type(pt_R) {}
   MeasureDefinition(MeasureType measure_type) : _tau_mode(UNDEFINED_SHAPE), _measure_type(measure_type) {}

   void setTauMode(TauMode tau_mode) {
      _tau_mode = tau_mode;
   }

   // added set measure method in case it becomes useful later -- TJW
   void setMeasureType(MeasureType measure_type) {
      _measure_type = measure_type;
   }

   // new methods added to generalize energy and angle squared for different measure types -- TJW
   double energy(const PseudoJet& jet) const;
   double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;

   bool has_denominator() const { return (_tau_mode == NORMALIZED_JET_SHAPE || _tau_mode == NORMALIZED_EVENT_SHAPE);}
   bool has_beam() const {return (_tau_mode == UNNORMALIZED_EVENT_SHAPE || _tau_mode == NORMALIZED_EVENT_SHAPE);}

   // added so that description will include the measure type -- TJW
   std::string measure_type_name() const {
      if (_measure_type == pt_R) return "pt_R";
      else if (_measure_type == E_theta) return "E_theta";
      else if (_measure_type == lorentz_dot) return "lorentz_dot";
      else return "Measure Type Undefined";
   }      

};


// updated all classes to remove TauMode argument in constructor, added default MeasureType argument to constructor -- TJW 
// flipped inheritance of NormalizedCutoffMeasure and NormalizedMeasure to avoid issue of squaring a maximum value -- TJW
// updated instances of perp() and jet_distance_squared() to account for more general metrics -- TJW

//------------------------------------------------------------------------
/// \class NormalizedCutoffMeasure
// This class is the default measure, based on the conical measure.
// This measure is defined as the pT of the particle multiplied by deltaR
// to the power of beta. This class includes the normalization factor determined by R0
class NormalizedCutoffMeasure : public MeasureDefinition {
   
public:
   
   NormalizedCutoffMeasure(double beta, double R0, double Rcutoff/*, TauMode tau_mode = NORMALIZED_EVENT_SHAPE*/, MeasureType measure_type = pt_R)
   : MeasureDefinition(measure_type), _beta(beta), _R0(R0), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
      setTauMode(NORMALIZED_EVENT_SHAPE);
   }
   
   virtual std::string description() const;
   
   virtual NormalizedCutoffMeasure* create() const {return new NormalizedCutoffMeasure(*this);}

   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return particle.squared_distance(axis);
   }
   
   virtual double beam_distance_squared(const fastjet::PseudoJet& /*particle*/) const {
      // return sq(_Rcutoff);
      return _RcutoffSq;
   }
   
   // temporarily trivially defined -- TJW
   virtual double axes_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return 1.0;
   };

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
   virtual AxesRefiner* createAxesRefiner(int nPass) const;
   
protected:
   double _beta;
   double _R0;
   double _Rcutoff;
   double _RcutoffSq;   
   
};
   
//------------------------------------------------------------------------
/// \class NormalizedMeasure
// This measure is the same as NormalizedCutoffMeasure, with Rcutoff taken to infinity.
class NormalizedMeasure : public NormalizedCutoffMeasure {

public:

   NormalizedMeasure(double beta, double R0/*, TauMode tau_mode = NORMALIZED_JET_SHAPE*/, MeasureType measure_type = pt_R)
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
class UnnormalizedCutoffMeasure : public NormalizedCutoffMeasure {
   
public:
   // Since all methods are identical, UnnormalizedMeasure inherits directly
   // from NormalizedMeasure. R0 is a dummy value since the value of R0 is unecessary for this class,
   // and the "false" flag sets _has_denominator in MeasureDefinition to false so no denominator is used.
   UnnormalizedCutoffMeasure(double beta, double Rcutoff/*,TauMode tau_mode = UNNORMALIZED_EVENT_SHAPE*/, MeasureType measure_type = pt_R)
   : NormalizedCutoffMeasure(beta, std::numeric_limits<double>::quiet_NaN(), Rcutoff, measure_type) {
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
   UnnormalizedMeasure(double beta/*,TauMode tau_mode = UNNORMALIZED_JET_SHAPE*/, MeasureType measure_type = pt_R)
   : UnnormalizedCutoffMeasure(beta, std::numeric_limits<double>::max(), measure_type) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(UNNORMALIZED_JET_SHAPE);
   }


   virtual std::string description() const;

   virtual UnnormalizedMeasure* create() const {return new UnnormalizedMeasure(*this);}
   
};
      
// haven't put in new definitions of e+e- variables for Geometric Measure because the naive redefinition didn't work -- TJW

//------------------------------------------------------------------------
/// \class GeometricCutoffMeasure
// This class is the geometric measure.  This measure is defined by the Lorentz dot product between
// the particle and the axis.  This class does not include normalization of tau_N.
// NOTE:  This class is in flux and should not be used for production purposes.
class GeometricCutoffMeasure : public MeasureDefinition {

public:
   // Right now, we are hard coded for beam_beta = 1.0, but that will need to change
   GeometricCutoffMeasure(double jet_beta, double Rcutoff/*, TauMode tau_mode = UNNORMALIZED_EVENT_SHAPE*/, MeasureType measure_type = pt_R)
   :   MeasureDefinition(measure_type),
      _jet_beta(jet_beta), _beam_beta(1.0), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
      }

   virtual std::string description() const;
   
   virtual GeometricCutoffMeasure* create() const {return new GeometricCutoffMeasure(*this);}
   
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      // double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(energy(lightAxis)*energy(particle));
      return pseudoRsquared;
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      // return sq(_Rcutoff);
      return _RcutoffSq;
   }

   // temporarily trivially defined -- TJW
   virtual double axes_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const  {
      return 1.0;
   }

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      // double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(lightAxis.pt(),_beam_beta - 1.0);
      // return particle.pt() * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(energy(lightAxis),_beam_beta - 1.0);
      return energy(particle) * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      // double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(particle.pt()/particle.e(),_beam_beta - 1.0);
      // return particle.pt() * weight * std::pow(_Rcutoff,_jet_beta);
      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(energy(particle)/energy(particle),_beam_beta - 1.0);
      return energy(particle) * weight * std::pow(_Rcutoff,_jet_beta);
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
   // The minimization routine is GeometricAxesRefiner
   virtual AxesRefiner* createAxesRefiner(int nPass) const;
   
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

   
//------------------------------------------------------------------------
/// \class GeometricMeasure
// Same as GeometricCutoffMeasure, but with Rcutoff taken to infinity.
class GeometricMeasure : public GeometricCutoffMeasure {
   
public:
   GeometricMeasure(double beta/*, TauMode tau_mode = UNNORMALIZED_JET_SHAPE*/, MeasureType measure_type = pt_R)
   : GeometricCutoffMeasure(beta,std::numeric_limits<double>::max(),measure_type) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(UNNORMALIZED_JET_SHAPE);
   }

   virtual std::string description() const;

   virtual GeometricMeasure* create() const {return new GeometricMeasure(*this);}
};

// Classes below added by TJW

// ------------------------------------------------------------------------
// / \class ConicalGeometricCutoffMeasure
// This class is the Conical geometric measure.  This measure is defined by the Lorentz dot product between
// the particle and the axis normalized by the axis and particle pT, as well as a factor of cosh(y) to vary
// the rapidity depepdence of the beam.
// NOTE:  This class is in flux and should not be used for production purposes.
class ConicalGeometricCutoffMeasure : public MeasureDefinition {

public:
   ConicalGeometricCutoffMeasure(double jet_beta, double jet_gamma, double Rcutoff/*, TauMode tau_mode = UNNORMALIZED_EVENT_SHAPE*/, MeasureType measure_type = pt_R)
   :   MeasureDefinition(measure_type),
      _jet_beta(jet_beta), _jet_gamma(jet_gamma), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
         setTauMode(UNNORMALIZED_EVENT_SHAPE);
      }

   virtual std::string description() const;
   
   virtual ConicalGeometricCutoffMeasure* create() const {return new ConicalGeometricCutoffMeasure(*this);}
   
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      // double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(energy(lightAxis)*energy(particle));
      return pseudoRsquared;
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      // return sq(_Rcutoff);
      return _RcutoffSq;
   }

   // defined according to G(n,p) in the paper -- TJW
   virtual double axes_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double jet_distance_param = (_jet_beta == 2.0) ? 1.0 : std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0 - 1.0);
      double weight = (_jet_gamma == 1.0) ? 1.0 : std::pow(lightAxis.pt(),_jet_gamma - 1.0);
      return (1/lightAxis.pt()) * weight * jet_distance_param;
   }

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double weight = (_jet_gamma == 1.0) ? 1.0 : std::pow(lightAxis.pt(),_jet_gamma - 1.0);
      return particle.pt() * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
      // double weight = (_jet_gamma == 1.0) ? 1.0 : std::pow(energy(lightAxis),_jet_gamma - 1.0);
      // return energy(particle) * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      double weight = (_jet_gamma == 1.0) ? 1.0 : std::pow(particle.pt()/particle.e(),_jet_gamma - 1.0);
      return particle.pt() * weight * std::pow(_Rcutoff,_jet_beta);
      // double weight = (_jet_gamma == 1.0) ? 1.0 : std::pow(energy(particle)/energy(particle),_jet_gamma - 1.0);
      // return energy(particle) * weight * std::pow(_Rcutoff,_jet_beta);
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
   // The minimization routine is GeometricAxesRefiner
   virtual AxesRefiner* createAxesRefiner(int nPass) const;
   
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

//------------------------------------------------------------------------
/// \class GeometricMeasure
// Same as GeometricCutoffMeasure, but with Rcutoff taken to infinity.
class ConicalGeometricMeasure : public ConicalGeometricCutoffMeasure {
   
public:
   ConicalGeometricMeasure(double jet_beta, double jet_gamma/*, TauMode tau_mode = UNNORMALIZED_JET_SHAPE*/, MeasureType measure_type = pt_R)
   : ConicalGeometricCutoffMeasure(jet_beta, jet_gamma, std::numeric_limits<double>::max(),measure_type) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(UNNORMALIZED_JET_SHAPE);
   }

   virtual std::string description() const;

   virtual ConicalGeometricMeasure* create() const {return new ConicalGeometricMeasure(*this);}
};

// ------------------------------------------------------------------------
// / \class XConeCutoffMeasure
// This class is the XCone Measure.  This is the default measure for use with the
// XCone algorithm. It is identical to the conical geomtric measure but with gamma = 1.0.
// NOTE:  This class is in flux and should not be used for production purposes.
class XConeCutoffMeasure : public ConicalGeometricCutoffMeasure {

public:
   XConeCutoffMeasure(double jet_beta, double Rcutoff/*, TauMode tau_mode = UNNORMALIZED_EVENT_SHAPE*/, MeasureType measure_type = pt_R)
   :   ConicalGeometricCutoffMeasure(jet_beta, 1.0, Rcutoff, measure_type) { }

   virtual XConeCutoffMeasure* create() const {return new XConeCutoffMeasure(*this);}

};

class XConeMeasure : public ConicalGeometricMeasure {

public:
   XConeMeasure(double jet_beta/*, TauMode tau_mode = UNNORMALIZED_EVENT_SHAPE*/, MeasureType measure_type = pt_R)
   :   ConicalGeometricMeasure(jet_beta, 1.0, measure_type) { }

   virtual XConeMeasure* create() const {return new XConeMeasure(*this);}
};


   
} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_MEASUREDEFINITION_HH__
