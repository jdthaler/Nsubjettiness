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
//Update to move MeasureFunction class and definitions into separate .cc/.hh files

#ifndef __FASTJET_CONTRIB_MEASUREFUNCTION_HH__
#define __FASTJET_CONTRIB_MEASUREFUNCTION_HH__

#include "fastjet/PseudoJet.hh"
#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

inline double sq(double x) {return x*x;}

//classes moved from Njettiness in order to avoid cross-references in AxesFinder and MeasureFunctor -- TJW 12/28

// NsubParameters class removed since it is no longer necessary -- TJW 1/9
// NsubGeometricParameters class removed since it is no longer necessary -- TJW 1/10

//KmeansParameters class moved to AxesFinder.hh -- TJW 12/31

///////
//
// Measure Function
//
///////

//name of MeasureFunctor changed to MeasureFunction -- TJW 12/25
//base MeasureFunction class definition moved from Nsubjettiness.hh -- TJW 12/25

// class added by TJW 1/14
//------------------------------------------------------------------------
/// \class TauComponents
// This class creates a wrapper for the various tau/subtau values calculated in Njettiness. This class allows Njettiness access to these variables
// without ever having to do the calculation itself. It takes in subtau numerators and tau denominator from MeasureFunction
// and outputs tau numerator, and normalized tau and subtau.
class TauComponents {
   private:

      std::vector<double> _subtaus_numerator;
      double _tau_numerator;
      double _tau_denominator;

      std::vector<double> _subtaus_normalized;
      double _tau_normalized;

   public: 
      // empty constructor necessary to initialize tau_components in Njettiness; set correctly in Njettiness::getTau function -- TJW 1/15
      TauComponents() {
         _subtaus_numerator.resize(1, 0.0);
         _tau_denominator = 0;
         _tau_numerator = 0;
         _subtaus_normalized.resize(1, 0.0);
         _tau_normalized = 0;
      }

      TauComponents(std::vector<double> subtaus_numerator, double tau_denominator) : _subtaus_numerator(subtaus_numerator), _tau_denominator(tau_denominator) {
         _tau_numerator = 0.0;
         _tau_normalized = 0.0;
         _subtaus_normalized.resize(_subtaus_numerator.size(),0.0);
         for (unsigned j = 0; j < _subtaus_numerator.size(); j++) {
            _subtaus_normalized[j] = _subtaus_numerator[j]/_tau_denominator;
            _tau_numerator += _subtaus_numerator[j];
            _tau_normalized += _subtaus_normalized[j];
         }
      }

      std::vector<double> subtaus_numerator() { return _subtaus_numerator; }
      double tau_numerator() { return _tau_numerator; }
      double tau_denominator() { return _tau_denominator; }

      std::vector<double> subtaus_normalized() { return _subtaus_normalized; }
      double tau_normalized() { return _tau_normalized; }

};

//------------------------------------------------------------------------
/// \class MeasureFunction
// This class calculates the tau_N of a jet given a specific measure. 
// It serves as a base class for calculating tau_N according to different measures that are implemented 
// in further derived classes, but does not define a particular measure itself. -- comment added by TJ
class MeasureFunction {

   protected:
      //bool set by derived classes to choose whether or not to use the denominator -- TJW 1/7
      bool _has_denominator; 

      //MeasureFunction() {} //removed since it is overloaded with constructor below -- TJW 1/7

      //new constructor to allow _has_denominator to be set by derived classes -- TJW 1/7
      MeasureFunction(bool has_denominator = true) : _has_denominator(has_denominator) {} 

      //function to calculate unnormalized subTau values
      std::vector<double> subtaus_numerator(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes);

      //function to calculate normalization factor for tau and subTau
      double tau_denominator(const std::vector <fastjet::PseudoJet>& particles);

   public:
      virtual ~MeasureFunction(){}

      //These functions define the measure by which tau_N is calculated. -- comment added by TJW
      virtual bool do_cluster(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      virtual double distance(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      
      virtual double numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      virtual double denominator(const fastjet::PseudoJet& particle) = 0;

      // Functions below call the virtual functions to implement the desired measures
      
      // The calculation of these values has been moved over to TauComponents; result() will return all necessary TauComponents -- TJW 1/15
      //function to calculate unnormalized taus
      // double tau_numerator(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes);

      //return normalized  subTaus
      // std::vector<double> subtaus_normalized(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes); 
      
      //returns normalized tau
      // double tau_normalized(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes);
      
      // double tau(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) {
      //    return tau_normalized(particles,axes);
      // }

      // new function to return all of the necessary TauComponents for specific input particles and axes -- TJW 1/15
      TauComponents result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) {
         std::vector<double> _subtaus_numerator = subtaus_numerator(particles, axes);
         double _tau_denominator = tau_denominator(particles);
         TauComponents _tau_components(_subtaus_numerator, _tau_denominator);
         return _tau_components;
      }

};


// moved from Njettiness.hh -- TJW 12/28
// name changed to DefaultNormalizedMeasure -- TJW 1/7
// updated to use three separate parameters in constructor instead of NsubParameters, changed all constructor definitions in other files -- TJW 1/9
//------------------------------------------------------------------------
/// \class DefaultNormalizedMeasure
// This class is the default measure, inheriting from the class above. This class will calculate tau_N 
// of a jet according to this measure. This measure is defined as the pT of the particle multiplied by deltaR 
// to the power of beta. This class includes normalization -- comment added by TJW
class DefaultNormalizedMeasure : public MeasureFunction {

   private:
      double _beta;
      double _R0;
      double _Rcutoff;

   public:

      //new constructor to use three separate parameters intead of NsubParameters -- TJW 1/9
      DefaultNormalizedMeasure(double beta, double R0, double Rcutoff, bool normalized = true) : MeasureFunction(normalized), _beta(beta), _R0(R0), _Rcutoff(Rcutoff) {}
      //MeasureFunction(true) tells MeasureFunction to use the denominator, since this value is normalized. -- TJW

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

// Class added by TJW 1/7
//------------------------------------------------------------------------
/// \class DefaultUnnormalizedMeasure
// This class is the unnormalized default measure, inheriting from the class above. The only difference from above is that the denominator is defined to be
// 1.0 by setting _has_denominator to false. -- comment added by TJW
class DefaultUnnormalizedMeasure : public DefaultNormalizedMeasure {

   public:
      //removed NsubParameters from DefaultNormalizedMeasure constructor -- TJW 1/9
      DefaultUnnormalizedMeasure(double beta, double Rcutoff) : DefaultNormalizedMeasure(beta, NAN, Rcutoff, false) {}
      // Since all methods are identical, UnnormalizedMeasure inherits directly from NormalizedMeasure. R0 is defaulted to NAN since the value of R0 is unecessary for this class. 
      // the "false" flag sets _has_denominator in MeasureFunction to false so no denominator is used. -- TJW
      
};

// moved from Njettiness.hh -- TJW 12/28
//------------------------------------------------------------------------
/// \class GeometricMeasure
// This class is the geometic measure, inheriting from the class above. This class will calculate tau_N 
// of a jet according to this measure. This measure is defined purely geometrically as the distance between
// the particle and the axis. This class includes normalization of tau_N. -- comment added by TJW
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
