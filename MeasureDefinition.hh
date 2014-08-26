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
class NormalizedMeasure;         // (beta,R0)
class UnnormalizedMeasure;       // (beta)
class GeometricMeasure;          // (beta)
class NormalizedCutoffMeasure;   // (beta,R0,Rcutoff)
class UnnormalizedCutoffMeasure; // (beta,Rcutoff)
class GeometricCutoffMeasure;    // (beta,Rcutoff)
   
   
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
   //bool set by derived classes to choose whether or not to use the denominator
   bool _has_denominator;
   bool _has_beam;
   
   // This constructor allows _has_denominator to be set by derived classes
   MeasureDefinition(bool has_denominator = true, bool has_beam = true) : _has_denominator(has_denominator), _has_beam(has_beam) {
   }
   
   
};

//------------------------------------------------------------------------
/// \class NormalizedCutoffMeasure
// This class is the default measure, based on the conical measure.
// This measure is defined as the pT of the particle multiplied by deltaR
// to the power of beta. This class includes the normalization factor determined by R0
class NormalizedCutoffMeasure : public MeasureDefinition {
   
public:
   
   NormalizedCutoffMeasure(double beta, double R0, double Rcutoff, bool normalized = true)
   : MeasureDefinition(normalized), _beta(beta), _R0(R0), _Rcutoff(Rcutoff) {}
   
   virtual std::string description() const;
   
   virtual NormalizedCutoffMeasure* create() const {return new NormalizedCutoffMeasure(*this);}

   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return particle.squared_distance(axis);
   }
   
   virtual double beam_distance_squared(const fastjet::PseudoJet& /*particle*/) const {
      return sq(_Rcutoff);
   }
   
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const{
      return particle.perp() * std::pow(jet_distance_squared(particle,axis),_beta/2.0);
   }
   
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      return particle.perp() * std::pow(_Rcutoff,_beta);
   }
   
   virtual double denominator(const fastjet::PseudoJet& particle) const {
      return particle.perp() * std::pow(_R0,_beta);
   }

   // The minimization routine is ConicalAxesRefiner
   virtual AxesRefiner* createAxesRefiner(int nPass) const;
   
protected:
   double _beta;
   double _R0;
   double _Rcutoff;
   
   
};
   

//------------------------------------------------------------------------
/// \class NormalizedMeasure
// This measure is the same as NormalizedCutoffMeasure, with Rcutoff taken
// to infinity.
class NormalizedMeasure : public NormalizedCutoffMeasure {

public:
   NormalizedMeasure(double beta, double R0)
   : NormalizedCutoffMeasure(beta,R0,std::numeric_limits<double>::max()) {}
   
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
   UnnormalizedCutoffMeasure(double beta, double Rcutoff)
   : NormalizedCutoffMeasure(beta, std::numeric_limits<double>::quiet_NaN(), Rcutoff, false) {}

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
   UnnormalizedMeasure(double beta)
   : UnnormalizedCutoffMeasure(beta, std::numeric_limits<double>::max()) {}

   virtual std::string description() const;

   virtual UnnormalizedMeasure* create() const {return new UnnormalizedMeasure(*this);}
   
};
   


//------------------------------------------------------------------------
/// \class GeometricCutoffMeasure
// This class is the geometic measure.  This measure is defined by the Lorentz dot product between
// the particle and the axis.  This class does not include normalization of tau_N.
// NOTE:  This class is in flux and should not be used for production purposes.
class GeometricCutoffMeasure : public MeasureDefinition {

public:
   // Right now, we are hard coded for beam_beta = 1.0, but that will need to change
   GeometricCutoffMeasure(double jet_beta, double Rcutoff)
   :   MeasureDefinition(false), // doesn't have denominator
      _jet_beta(jet_beta), _beam_beta(1.0), _Rcutoff(Rcutoff) {}

   virtual std::string description() const;
   
   virtual GeometricCutoffMeasure* create() const {return new GeometricCutoffMeasure(*this);}
   
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      return pseudoRsquared;
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      return sq(_Rcutoff);
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
   virtual AxesRefiner* createAxesRefiner(int nPass) const;
   
protected:
   double _jet_beta;
   double _beam_beta;
   double _Rcutoff;
   
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
   GeometricMeasure(double beta)
   : GeometricCutoffMeasure(beta,std::numeric_limits<double>::max()) {}

   virtual std::string description() const;

   virtual GeometricMeasure* create() const {return new GeometricMeasure(*this);}
};
   
   
} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_MEASUREDEFINITION_HH__
