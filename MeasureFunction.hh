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

///////
//
// Parameter classes
//
///////
//classes moved from Njettiness in order to avoid cross-references in AxesFinder and MeasureFunctor -- TJW 12/28

//------------------------------------------------------------------------
/// \class NsubParameters
// These parameters define Njettiness. These characteristic parameters are used to define how tau_N is calculated
// using MeasureFunction -- comment added by TJ 
class NsubParameters {
private:
   double _beta;  // angular weighting exponent
   double _R0;    // characteristic jet radius (for normalization)
   double _Rcutoff;  // Cutoff scale for cone jet finding (default is large number such that no boundaries are used)
   
public:
   NsubParameters(const double mybeta, const double myR0, const double myRcutoff=10000.0) :
   _beta(mybeta), _R0(myR0), _Rcutoff(myRcutoff) {}
   double beta() const {return _beta;}
   double R0() const {return _R0;}
   double Rcutoff() const {return _Rcutoff;}
};
//------------------------------------------------------------------------
/// \class NsubGeometricParameters
// Parameters that define GeometricMeasure. This parameter is used in the classes that find axes through
// geometric minimization. -- comment added by TJ
class NsubGeometricParameters {
private:
   double _Rcutoff;  // Cutoff scale for cone jet finding (default is large number such that no boundaries are used)
   
public:
   NsubGeometricParameters(const double myRcutoff=10000.0) : _Rcutoff(myRcutoff) {}
   double Rcutoff() const {return _Rcutoff;}
};

//KmeansParameters class moved to AxesFinder.hh -- TJW 12/31

///////
//
// Measure Function
//
///////

//name of MeasureFunctor changed to MeasureFunction -- TJW 12/25
//base MeasureFunction class definition moved from Nsubjettiness.hh -- TJW 12/25

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

   public:
      virtual ~MeasureFunction(){}

      //These functions define the measure by which tau_N is calculated. -- comment added by TJW
      virtual bool do_cluster(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      virtual double distance(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      
      virtual double numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) = 0;
      virtual double denominator(const fastjet::PseudoJet& particle) = 0;

      // Functions below call the virtual functions to implement the desired measures

      //function to calculate unnormalized subTau values
      std::vector<double> subtaus_numerator(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes);

      //function to calculate normalization factor for tau and subTau
      double tau_denominator(const std::vector <fastjet::PseudoJet>& particles);
      
      // These functions are built out of the above functions. 
      
      //function to calculate unnormalized taus
      double tau_numerator(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes);

      //return normalized  subTaus
      std::vector<double> subtaus_normalized(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes); 
      
      //returns normalized tau
      double tau_normalized(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes);
      
      double tau(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) {
         return tau_normalized(particles,axes);
      }

};

// moved from Njettiness.hh -- TJW 12/28
// name changed to DefaultNormalizedMeasure -- TJW 1/7
//------------------------------------------------------------------------
/// \class DefaultNormalizedMeasure
// This class is the default measure, inheriting from the class above. This class will calculate tau_N 
// of a jet according to this measure. This measure is defined as the pT of the particle multiplied by deltaR 
// to the power of beta. This class includes normalization -- comment added by TJW
class DefaultNormalizedMeasure : public MeasureFunction {

   private:
      NsubParameters _paraNsub;

   public:
      DefaultNormalizedMeasure(NsubParameters paraNsub, bool normalized = true) : MeasureFunction(normalized), _paraNsub(paraNsub) {}
      //MeasureFunction(true) tells MeasureFunction to use the denominator, since this value is normalized. -- TJW

      virtual bool do_cluster(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         return (distance(particle,axis) <= _paraNsub.Rcutoff());
      }

      virtual double distance(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         return std::sqrt(particle.squared_distance(axis));
      }

      virtual double numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) {
         double deltaR = std::sqrt(particle.squared_distance(axis));
         if (deltaR > _paraNsub.Rcutoff()) deltaR = _paraNsub.Rcutoff();
         return particle.perp() * std::pow(deltaR,_paraNsub.beta());
      }

      virtual double denominator(const fastjet::PseudoJet& particle) {
         return particle.perp() * std::pow(_paraNsub.R0(),_paraNsub.beta());
      }

};

// Class added by TJW 1/7
//------------------------------------------------------------------------
/// \class DefaultUnnormalizedMeasure
// This class is the unnormalized default measure, inheriting from the class above. The only difference from above is that the denominator is defined to be
// 1.0 by setting _has_denominator to false. -- comment added by TJW
class DefaultUnnormalizedMeasure : public DefaultNormalizedMeasure {

   public:
      DefaultUnnormalizedMeasure(double beta, double Rcutoff) : DefaultNormalizedMeasure(NsubParameters(beta, NAN, Rcutoff), false) {}
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
