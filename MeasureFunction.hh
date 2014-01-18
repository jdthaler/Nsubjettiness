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

#ifndef __FASTJET_CONTRIB_MEASUREFUNCTION_HH__
#define __FASTJET_CONTRIB_MEASUREFUNCTION_HH__

#include "fastjet/PseudoJet.hh"
#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

inline double sq(double x) {return x*x;}

   ///////
//
// Measure Function
//
///////

/// \class TauComponents
// This class creates a wrapper for the various tau/subtau values calculated in Njettiness. This class allows Njettiness access to these variables
// without ever having to do the calculation itself. It takes in subtau numerators and tau denominator from MeasureFunction
// and outputs tau numerator, and normalized tau and subtau.
// TODO:  Consider merging with NjettinessExtras.  Add axes information.
class TauComponents {
   private:

      // these values are input in the constructor
      std::vector<double> _subtaus_numerator;
      double _tau_denominator;
      bool _has_denominator; //added so that TauComponents knows if denominator is used or not

      // these values are derived from above values
      std::vector<double> _subtaus_normalized;
      double _tau_numerator;
      double _tau_normalized;


   public: 
      // empty constructor necessary to initialize tau_components in Njettiness
      // later set correctly in Njettiness::getTau function
      TauComponents() {
         _subtaus_numerator.resize(1, 0.0);
         _tau_denominator = 0;
         _tau_numerator = 0;
         _subtaus_normalized.resize(1, 0.0);
         _tau_normalized = 0;
         _has_denominator = false;
      }

      // This constructor takes input vector and double and calculates all necessary tau components
      TauComponents(std::vector<double> subtaus_numerator, double tau_denominator, bool has_denominator)
         : _subtaus_numerator(subtaus_numerator), _tau_denominator(tau_denominator), _has_denominator(has_denominator) {
         if (!_has_denominator) assert(tau_denominator == 1); //make sure that tau_denominator is 1 if _has_denominator is flagged
         _tau_numerator = 0.0;
         _tau_normalized = 0.0;
         _subtaus_normalized.resize(_subtaus_numerator.size(),0.0);
         for (unsigned j = 0; j < _subtaus_numerator.size(); j++) {
            _subtaus_normalized[j] = _subtaus_numerator[j]/_tau_denominator;
            _tau_numerator += _subtaus_numerator[j];
            _tau_normalized += _subtaus_normalized[j];
         }
      }


      // return values
      std::vector<double> subtaus_numerator() const { return _subtaus_numerator; }
      double tau_denominator() const { return _tau_denominator; }
      bool has_denominator() const { return _has_denominator; }
   
      double tau_numerator() const { return _tau_numerator; }
      std::vector<double> subtaus_normalized() const { return _subtaus_normalized; }
      double tau_normalized() const { return _tau_normalized; }

      //separate function for ease of user
      double tau() const { return _tau_normalized; }
   
};

//------------------------------------------------------------------------
/// \class MeasureFunction
// This class calculates the tau_N of a jet given a specific measure. 
// It serves as a base class for calculating tau_N according to different measures that are implemented 
// in further derived classes, but does not define a particular measure itself.
class MeasureFunction {

   protected:
      //bool set by derived classes to choose whether or not to use the denominator
      bool _has_denominator; 

      // This constructor allows _has_denominator to be set by derived classes
      MeasureFunction(bool has_denominator = true) : _has_denominator(has_denominator) {} 

      // helper function to calculate unnormalized subTau values
      std::vector<double> subtaus_numerator(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes);

      // helper function to calculate normalization factor for tau and subTau
      double tau_denominator(const std::vector <fastjet::PseudoJet>& particles);

   public:
      virtual ~MeasureFunction(){}

      //These functions define the measure by which tau_N is calculated,
      //and are overloaded by the various measures below
      virtual bool do_cluster(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      virtual double distance(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      
      virtual double numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      virtual double denominator(const fastjet::PseudoJet& particle) = 0;

      // Return all of the necessary TauComponents for specific input particles and axes
      TauComponents result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes);

};


/// \class DefaultNormalizedMeasure
// This class is the default measure, inheriting from the class above. This class will calculate tau_N 
// of a jet according to this measure. This measure is defined as the pT of the particle multiplied by deltaR 
// to the power of beta. This class includes the normalization factor determined by R0
class DefaultNormalizedMeasure : public MeasureFunction {

   private:
      double _beta;
      double _R0;
      double _Rcutoff;

   public:

      DefaultNormalizedMeasure(double beta, double R0, double Rcutoff, bool normalized = true)
      : MeasureFunction(normalized), _beta(beta), _R0(R0), _Rcutoff(Rcutoff) {}

      virtual bool do_cluster(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         return (distance(particle,axis) <= _Rcutoff);
      }

      virtual double distance(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         return std::sqrt(particle.squared_distance(axis));
      }

      virtual double numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         double deltaR = std::sqrt(particle.squared_distance(axis));
         if (deltaR > _Rcutoff) deltaR = _Rcutoff;
         return particle.perp() * std::pow(deltaR,_beta);
      }

      virtual double denominator(const fastjet::PseudoJet& particle) {
         return particle.perp() * std::pow(_R0,_beta);
      }

};

//------------------------------------------------------------------------
/// \class DefaultUnnormalizedMeasure
// This class is the unnormalized default measure, inheriting from the class above. The only difference from above
// is that the denominator is defined to be 1.0 by setting _has_denominator to false.
class DefaultUnnormalizedMeasure : public DefaultNormalizedMeasure {

   public:
      // Since all methods are identical, UnnormalizedMeasure inherits directly from NormalizedMeasure. R0 is defaulted to NAN since the value of R0 is unecessary for this class.
      // the "false" flag sets _has_denominator in MeasureFunction to false so no denominator is used.
      DefaultUnnormalizedMeasure(double beta, double Rcutoff) : DefaultNormalizedMeasure(beta, NAN, Rcutoff, false) {}

      
};

//------------------------------------------------------------------------
/// \class GeometricMeasure
// This class is the geometic measure, inheriting from the class above. This class will calculate tau_N 
// of a jet according to this measure. This measure is defined by the Lorentz dot product between
// the particle and the axis. This class includes normalization of tau_N.
class GeometricMeasure : public MeasureFunction {

   private:
      double _Rcutoff;

      // create light-like axis
      fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
         double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
         return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
      }
      
      // get Q value
      double qValueOf(const fastjet::PseudoJet& input) const {
         return lightFrom(input).pt();
      }

   public:
      GeometricMeasure(double Rcutoff) : _Rcutoff(Rcutoff) {}
            
      virtual bool do_cluster(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         return (distance(particle,axis) <= _Rcutoff);
      }
      
      virtual double distance(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         double axisValue = 2.0*dot_product(lightFrom(axis),particle)/qValueOf(axis);
         double beamValue = particle.pt();
         double pseudoR = std::sqrt(axisValue/beamValue);
         return pseudoR;
      }

      virtual double numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         double pseudoR = distance(particle,axis);

         if (pseudoR > _Rcutoff) {
            pseudoR = _Rcutoff;
         }
         return particle.pt()*sq(pseudoR);
      }

      virtual double denominator(const fastjet::PseudoJet& particle) {
         return 1.0;
      }
};


} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_MEASUREFUNCTION_HH__
