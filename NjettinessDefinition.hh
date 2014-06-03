//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: Njettiness.hh 653 2014-06-02 08:49:04Z jthaler $
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

#ifndef __FASTJET_CONTRIB_NJETTINESS_DEFINITION_HH__
#define __FASTJET_CONTRIB_NJETTINESS_DEFINITION_HH__


#include "MeasureFunction.hh"
#include "AxesFinder.hh"

#include "fastjet/PseudoJet.hh"
#include <fastjet/LimitedWarning.hh>

#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
   
// Eventually this file will contain an AxesDefinition as well as an NjettinessDefinition that combines an AxesDefinition with a MeasureDefinition.  Also, eventually these will report their names
   
///////
//
// MeasureDefinition
//
///////

//The MeasureDefinition is a way to set the MeasureMode and
//parameters used by Njettiness, with errors in the number of parameters given at
//compile time (instead of at run time).  The MeasureDefintion knows which core objects
//to call to make the measurement, as well as properties of the MeasureFunction
class MeasureDefinition {
   
public:

   virtual ~MeasureDefinition() {};
   
   // In derived classes, this should return a copy of the corresponding
   // derived class
   virtual MeasureDefinition* copy() const = 0;
   
   //Return the MeasureFunction corresponding to this definition
   virtual MeasureFunction* associatedMeasureFunction() const = 0;
   
   //Return the AxesFinder that should be used for one-pass minimization
   virtual AxesFinder* associatedOnePassAxesFinder(AxesFinder* startingFinder) const = 0;

   //Specify whether multi-pass minimization makes sense, and if so, return the AxesFinder that should be used for multi-pass minimization
   virtual bool canDoMultiPassMinimization() const = 0;
      virtual AxesFinder* associatedMultiPassAxesFinder(AxesFinder* startingFinder) const = 0;
   
};

// The normalized measure, with two parameters: beta and R0
class NormalizedMeasure : public MeasureDefinition {
private:
   double _beta;
   double _R0;
   
public:
   NormalizedMeasure(double beta, double R0)
   : _beta(beta), _R0(R0) {}
   
   virtual MeasureDefinition* copy() const {return new NormalizedMeasure(*this);}
   
   virtual MeasureFunction* associatedMeasureFunction() const { return new DefaultNormalizedMeasureFunction(_beta,_R0,std::numeric_limits<double>::max()); }
   virtual AxesFinder* associatedOnePassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromOnePassMinimization(startingFinder, _beta, std::numeric_limits<double>::max());}
   virtual AxesFinder* associatedMultiPassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromKmeansMinimization(startingFinder, _beta, std::numeric_limits<double>::max(),100);}

   virtual bool canDoMultiPassMinimization() const { return true; }
};
   
// The unnormalized measure, with just one parameter: beta
class UnnormalizedMeasure : public MeasureDefinition {
private:
   double _beta;
   
public:
   UnnormalizedMeasure(double beta)
   : _beta(beta) {}

   virtual MeasureDefinition* copy() const {return new UnnormalizedMeasure(*this);}
   
   virtual MeasureFunction* associatedMeasureFunction() const { return new DefaultUnnormalizedMeasureFunction(_beta,std::numeric_limits<double>::max());}
   virtual AxesFinder* associatedOnePassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromOnePassMinimization(startingFinder, _beta, std::numeric_limits<double>::max());}
   virtual AxesFinder* associatedMultiPassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromKmeansMinimization(startingFinder, _beta, std::numeric_limits<double>::max(),100);}

   virtual bool canDoMultiPassMinimization() const { return true; }

};
   
   
// The geometric  measure, with 1 parameter: beta
// This measure is still evolving and shouldn't be used for any critial applications yet
class GeometricMeasure : public MeasureDefinition {
private:
   double _beta;
   
public:
   GeometricMeasure(double beta)
   : _beta(beta) {}
   
   virtual MeasureDefinition* copy() const {return new GeometricMeasure(*this);}
   
   virtual MeasureFunction* associatedMeasureFunction() const { return new GeometricMeasureFunction(_beta,std::numeric_limits<double>::max());}
   virtual AxesFinder* associatedOnePassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromGeometricMinimization(startingFinder, _beta, std::numeric_limits<double>::max());}
   virtual AxesFinder* associatedMultiPassAxesFinder(AxesFinder* /*startingFinder*/) const {
      throw Error("GeometricMeasure does not support multi-pass minimization.");
      return NULL;
   }
   
   virtual bool canDoMultiPassMinimization() const { return false; }
};


// The normalized cutoff measure, with 3 parameters: beta, R0, Rcutoff
class NormalizedCutoffMeasure : public MeasureDefinition {
private:
   double _beta;
   double _R0;
   double _Rcutoff;
   
public:
   NormalizedCutoffMeasure(double beta, double R0, double Rcutoff)
   : _beta(beta), _R0(R0), _Rcutoff(Rcutoff) {}

   virtual MeasureDefinition* copy() const {return new NormalizedCutoffMeasure(*this);}
   
   virtual MeasureFunction* associatedMeasureFunction() const { return new DefaultNormalizedMeasureFunction(_beta,_R0,_Rcutoff); }
   virtual AxesFinder* associatedOnePassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromOnePassMinimization(startingFinder, _beta, _Rcutoff);}
   virtual AxesFinder* associatedMultiPassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromKmeansMinimization(startingFinder, _beta, std::numeric_limits<double>::max(),100);}
  
   virtual bool canDoMultiPassMinimization() const { return true; }

};

// The unnormalized cutoff measure, with 2 parameters: beta, Rcutoff
class UnnormalizedCutoffMeasure : public MeasureDefinition {
private:
   double _beta;
   double _Rcutoff;
   
public:
   UnnormalizedCutoffMeasure(double beta, double Rcutoff)
   : _beta(beta), _Rcutoff(Rcutoff) {}
   
   virtual MeasureDefinition* copy() const {return new UnnormalizedCutoffMeasure(*this);}

   virtual MeasureFunction* associatedMeasureFunction() const { return new DefaultUnnormalizedMeasureFunction(_beta,_Rcutoff);}
   virtual AxesFinder* associatedOnePassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromOnePassMinimization(startingFinder, _beta, _Rcutoff);}
   virtual AxesFinder* associatedMultiPassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromKmeansMinimization(startingFinder, _beta, std::numeric_limits<double>::max(),100);}

   virtual bool canDoMultiPassMinimization() const { return true; }

};

// The Geometric  measure, with 2 parameters: beta, Rcutoff
// This measure is still evolving and shouldn't be used for any critial applications yet
class GeometricCutoffMeasure : public MeasureDefinition {
private:
   double _beta;
   double _Rcutoff;
   
public:
   GeometricCutoffMeasure(double beta, double Rcutoff)
   : _beta(beta), _Rcutoff(Rcutoff) {}

   virtual MeasureDefinition* copy() const {return new GeometricCutoffMeasure(*this);}
   
   virtual MeasureFunction* associatedMeasureFunction() const { return new GeometricMeasureFunction(_beta,_Rcutoff);}
   virtual AxesFinder* associatedOnePassAxesFinder(AxesFinder* startingFinder) const { return new AxesFinderFromGeometricMinimization(startingFinder, _beta, _Rcutoff);}
   virtual AxesFinder* associatedMultiPassAxesFinder(AxesFinder* /*startingFinder*/) const {
      throw Error("GeometricCutoffMeasure does not support multi-pass minimization.");
      return NULL;
   }
   
   virtual bool canDoMultiPassMinimization() const { return false; }

};
      
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

