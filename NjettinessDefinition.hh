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

#ifndef __FASTJET_CONTRIB_NJETTINESS_DEFINITION_HH__
#define __FASTJET_CONTRIB_NJETTINESS_DEFINITION_HH__


#include "MeasureFunction.hh"
#include "AxesFinder.hh"

#include "fastjet/PseudoJet.hh"
#include <fastjet/LimitedWarning.hh>

#include <iomanip>
#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
   
// Eventually this file might contains an NjettinessDefinition that combines an
// AxesDefinition with a MeasureDefinition. It's not clear how useful that would be, though
  
// The following Measures are available (and the relevant arguments):
class NormalizedMeasure;         // (beta,R0)
class UnnormalizedMeasure;       // (beta)
class GeometricMeasure;          // (beta)
class NormalizedCutoffMeasure;   // (beta,R0,Rcutoff)
class UnnormalizedCutoffMeasure; // (beta,Rcutoff)
class GeometricCutoffMeasure;    // (beta,Rcutoff)

// The following Axes are available (and the relevant arguments, if needed)
class KT_Axes;
class CA_Axes;
class AntiKT_Axes;   // (R0)
class WTA_KT_Axes;
class WTA_CA_Axes;
class Manual_Axes;
class OnePass_KT_Axes;
class OnePass_CA_Axes;
class OnePass_AntiKT_Axes;   // (R0)
class OnePass_WTA_KT_Axes;
class OnePass_WTA_CA_Axes;
class OnePass_Manual_Axes;
class MultiPass_Axes;

// Below are just technical implementations of the variable axes and measures.
  
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

   // Description of measure and parameters
   virtual std::string description() const = 0;
   
   // In derived classes, this should return a copy of the corresponding
   // derived class
   virtual MeasureDefinition* create() const = 0;
   
   //Return the MeasureFunction corresponding to this definition
   SharedPtr<MeasureFunction> measureFunction() const { return _measureFunction;}
   
   //Return the AxesFinder that should be used for minimization
   SharedPtr<RefiningAxesFinder> refiningAxesFinder() const {return _refiningAxesFinder;}
      
   virtual ~MeasureDefinition(){};
   
protected:
   
   SharedPtr<MeasureFunction> _measureFunction;
   SharedPtr<RefiningAxesFinder> _refiningAxesFinder;
   
};

// The normalized measure, with two parameters: beta and R0
class NormalizedMeasure : public MeasureDefinition {
   
public:
   NormalizedMeasure(double beta, double R0)
   : _beta(beta), _R0(R0) {
      _measureFunction.reset(new ConicalNormalizedMeasureFunction(_beta,_R0,std::numeric_limits<double>::max()));
      _refiningAxesFinder.reset(new AxesFinderFromConicalMinimization(_beta, std::numeric_limits<double>::max()));
   }
   
   virtual std::string description() const;

   virtual NormalizedMeasure* create() const {return new NormalizedMeasure(*this);}

private:
   double _beta;
   double _R0;
};
   
// The unnormalized measure, with just one parameter: beta
class UnnormalizedMeasure : public MeasureDefinition {

public:
   UnnormalizedMeasure(double beta)
   : _beta(beta) {
      _measureFunction.reset(new ConicalUnnormalizedMeasureFunction(_beta,std::numeric_limits<double>::max()));
      _refiningAxesFinder.reset(new AxesFinderFromConicalMinimization(_beta, std::numeric_limits<double>::max()));
   }

   virtual UnnormalizedMeasure* create() const {return new UnnormalizedMeasure(*this);}
   
   virtual std::string description() const;

private:
   double _beta;
};
   
   
// The geometric  measure, with 1 parameter: beta
// This measure is still evolving and shouldn't be used for any critial applications yet
class GeometricMeasure : public MeasureDefinition {
   
public:
   GeometricMeasure(double beta)
   : _beta(beta) {
      _measureFunction.reset(new GeometricMeasureFunction(_beta,std::numeric_limits<double>::max()));
      _refiningAxesFinder.reset(new AxesFinderFromGeometricMinimization(_beta, std::numeric_limits<double>::max()));
   }
   
   virtual GeometricMeasure* create() const {return new GeometricMeasure(*this);}
   
   virtual std::string description() const;
   
private:
   double _beta;

};


// The normalized cutoff measure, with 3 parameters: beta, R0, Rcutoff
class NormalizedCutoffMeasure : public MeasureDefinition {

public:
   NormalizedCutoffMeasure(double beta, double R0, double Rcutoff)
   : _beta(beta), _R0(R0), _Rcutoff(Rcutoff) {
      _measureFunction.reset(new ConicalNormalizedMeasureFunction(_beta,_R0,_Rcutoff));
      _refiningAxesFinder.reset(new AxesFinderFromConicalMinimization(_beta, _Rcutoff));
   }

   virtual std::string description() const;
   
   virtual NormalizedCutoffMeasure* create() const {return new NormalizedCutoffMeasure(*this);}
   
private:
   double _beta;
   double _R0;
   double _Rcutoff;
   
};

// The unnormalized cutoff measure, with 2 parameters: beta, Rcutoff
class UnnormalizedCutoffMeasure : public MeasureDefinition {
   
public:
   UnnormalizedCutoffMeasure(double beta, double Rcutoff)
   : _beta(beta), _Rcutoff(Rcutoff) {
      _measureFunction.reset(new ConicalUnnormalizedMeasureFunction(_beta,_Rcutoff));
      _refiningAxesFinder.reset(new AxesFinderFromConicalMinimization(_beta, _Rcutoff));
   }

   virtual std::string description() const;
   
   virtual UnnormalizedCutoffMeasure* create() const {return new UnnormalizedCutoffMeasure(*this);}

private:
   double _beta;
   double _Rcutoff;
   
};

// The Geometric  measure, with 2 parameters: beta, Rcutoff
// This measure is still evolving and shouldn't be used for any critial applications yet
class GeometricCutoffMeasure : public MeasureDefinition {
   
public:
   GeometricCutoffMeasure(double beta, double Rcutoff)
   : _beta(beta), _Rcutoff(Rcutoff) {
      _measureFunction.reset(new GeometricMeasureFunction(_beta,_Rcutoff));
      _refiningAxesFinder.reset(new AxesFinderFromGeometricMinimization(_beta, _Rcutoff));
   }

   virtual GeometricCutoffMeasure* create() const {return new GeometricCutoffMeasure(*this);}

   virtual std::string description() const;

private:
   double _beta;
   double _Rcutoff;
   
};
   
   
///////
//
// AxesDefinition
//
///////
  
  
//Analogous to MeasureDefinition, AxesDefinition defines which AxesFinder to use
//At the moment, most AxesDefinition do not have an arugment (except the Anti-KT ones)
class AxesDefinition {
   
public:

   enum AxesRefiningMode { // define the cases of zero pass and one pass for convenience
      NO_REFINING,
      ONE_PASS,
      MULTI_PASS,
   };
   
   AxesDefinition(AxesRefiningMode mode, int nPass = 0) : _refiningMode(mode), _Npass(nPass) {}
   
   // description of axes (and any parameters)
   virtual std::string short_description() const = 0;
   virtual std::string description() const = 0;
   
   // This has to be defined in all derived classes
   virtual AxesDefinition* create() const = 0;

   
   // These describe how the axes finder works
   bool givesRandomizedResults() const {
      return (_refiningMode == MULTI_PASS);
   }
   
   bool needsManualAxes() const {
      return (_startingAxesFinder() == NULL); // if there is no starting axes finder
   }
   
   // return starting Axes finder already defined
   SharedPtr<StartingAxesFinder> startingAxesFinder() const { return _startingAxesFinder; }
   
   AxesRefiningMode refiningMode() const { return _refiningMode; }
   
   int nPass() const { return _Npass; }

   virtual ~AxesDefinition() {};
   
protected:
   
   SharedPtr<StartingAxesFinder> _startingAxesFinder;
   AxesRefiningMode _refiningMode;
   int _Npass;
   
};
  
  
// kt axes
class KT_Axes : public AxesDefinition {
public:
   KT_Axes() : AxesDefinition(NO_REFINING) {
      _startingAxesFinder.reset(new AxesFinderFromKT());
   }

   virtual std::string short_description() const {
      return "KT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "KT Axes";
      return stream.str();
   };
   
   virtual KT_Axes* create() const {return new KT_Axes(*this);}

};

// ca axes
class CA_Axes : public AxesDefinition {
public:
   CA_Axes() : AxesDefinition(NO_REFINING) {
      _startingAxesFinder.reset(new AxesFinderFromCA());
   }

   virtual std::string short_description() const {
      return "CA";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "CA Axes";
      return stream.str();
   };
   
   virtual CA_Axes* create() const {return new CA_Axes(*this);}
   
};

// anti-kt axes, one parameter R0 is subjet radius
class AntiKT_Axes : public AxesDefinition {

public:
   AntiKT_Axes(double R0 = 0.2) : AxesDefinition(NO_REFINING), _R0(R0) {
      _startingAxesFinder.reset(new AxesFinderFromAntiKT(_R0));
   }

   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "AKT" << _R0;
      return stream.str();
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Anti-KT Axes (R0 = " << _R0 << ")";
      return stream.str();
   };
   
   virtual AntiKT_Axes* create() const {return new AntiKT_Axes(*this);}
   
private:
   double _R0;

};

// winner-take-all recombination with kt axes
class WTA_KT_Axes : public AxesDefinition {
public:
   WTA_KT_Axes() : AxesDefinition(NO_REFINING) {
      _startingAxesFinder.reset(new AxesFinderFromWTA_KT());
   }

   virtual std::string short_description() const {
      return "WTA KT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Winner-Take-All KT Axes";
      return stream.str();
   };
   
   virtual WTA_KT_Axes* create() const {return new WTA_KT_Axes(*this);}

   
};

// winner-take-all recombination with CA axes
class WTA_CA_Axes : public AxesDefinition {
public:
   WTA_CA_Axes() : AxesDefinition(NO_REFINING) {
      _startingAxesFinder.reset(new AxesFinderFromWTA_CA());
   }

   virtual std::string short_description() const {
      return "WTA CA";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Winner-Take-All CA Axes";
      return stream.str();
   };
   
   virtual WTA_CA_Axes* create() const {return new WTA_CA_Axes(*this);}


};
   
// Onepass minimization from kt axes
class OnePass_KT_Axes : public AxesDefinition {
public:
   OnePass_KT_Axes() : AxesDefinition(ONE_PASS) {
      _startingAxesFinder.reset(new AxesFinderFromKT());
   }
   
   virtual std::string short_description() const {
      return "OnePass KT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from KT Axes";
      return stream.str();
   };
   
   virtual OnePass_KT_Axes* create() const {return new OnePass_KT_Axes(*this);}
   

};

// Onepass minimization from CA axes
class OnePass_CA_Axes : public AxesDefinition {
public:
   OnePass_CA_Axes() : AxesDefinition(ONE_PASS) {
      _startingAxesFinder.reset(new AxesFinderFromCA());
   }

   virtual std::string short_description() const {
      return "OnePass CA";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from CA Axes";
      return stream.str();
   };
   
   virtual OnePass_CA_Axes* create() const {return new OnePass_CA_Axes(*this);}


};

// Onepass minimization from AntiKT axes, one parameter R0
class OnePass_AntiKT_Axes : public AxesDefinition {

public:
   OnePass_AntiKT_Axes(double R0) : AxesDefinition(ONE_PASS), _R0(R0) {
      _startingAxesFinder.reset(new AxesFinderFromAntiKT(_R0));
   }
   
   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "OnePassAKT" << _R0;
      return stream.str();
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Anti-KT Axes (R0 = " << _R0 << ")";
      return stream.str();
   };
   
   virtual OnePass_AntiKT_Axes* create() const {return new OnePass_AntiKT_Axes(*this);}



private:
   double _R0;
};

// Onepass minimization from winner-take-all kt axes
class OnePass_WTA_KT_Axes : public AxesDefinition {
public:
   OnePass_WTA_KT_Axes() : AxesDefinition(ONE_PASS) {
      _startingAxesFinder.reset(new AxesFinderFromWTA_KT());
   }
   
   virtual std::string short_description() const {
      return "OnePass WTA KT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Winner-Take-All KT Axes";
      return stream.str();
   };
   
   virtual OnePass_WTA_KT_Axes* create() const {return new OnePass_WTA_KT_Axes(*this);}
   

};


// Onepass minimization from winner-take-all CA axes
class OnePass_WTA_CA_Axes : public AxesDefinition {
   
public:
   OnePass_WTA_CA_Axes() : AxesDefinition(ONE_PASS) {
      _startingAxesFinder.reset(new AxesFinderFromWTA_CA());
   }

   virtual std::string short_description() const {
      return "OnePass WTA CA";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Winner-Take-All CA Axes";
      return stream.str();
   };
   
   virtual OnePass_WTA_CA_Axes* create() const {return new OnePass_WTA_CA_Axes(*this);}
   


   
};
   
// set axes manually
class Manual_Axes : public AxesDefinition {
public:
   Manual_Axes() : AxesDefinition(NO_REFINING) {
      _startingAxesFinder.reset(); // NULL
   }
   
   virtual std::string short_description() const {
      return "Manual";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Manual Axes";
      return stream.str();
   };
   
   virtual Manual_Axes* create() const {return new Manual_Axes(*this);}


};

// one pass minimization from manual starting point
class OnePass_Manual_Axes : public AxesDefinition {
public:
   OnePass_Manual_Axes() : AxesDefinition(ONE_PASS) {
      _startingAxesFinder.reset(); // NULL
   }

   virtual std::string short_description() const {
      return "OnePass Manual";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Manual Axes";
      return stream.str();
   };
   
   virtual OnePass_Manual_Axes* create() const {return new OnePass_Manual_Axes(*this);}

   
};
   
// multi-pass minimization from kT starting point
class MultiPass_Axes : public AxesDefinition {

public:
   MultiPass_Axes(unsigned int Npass) : AxesDefinition(MULTI_PASS,Npass) {
      _startingAxesFinder.reset(new AxesFinderFromKT());
   }

   virtual std::string short_description() const {
      return "MultiPass";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Multi-Pass Axes (Npass = " << _Npass << ")";
      return stream.str();
   };
   
   virtual MultiPass_Axes* create() const {return new MultiPass_Axes(*this);}

   virtual bool givesRandomizedResults() const {return true;}
   
};
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

