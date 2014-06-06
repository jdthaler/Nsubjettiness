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

#include <iomanip>
#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
   
// Eventually this file might contains an NjettinessDefinition that combines an
// AxesDefinition with a MeasureDefinition. It's not clear how useful that would be, though
   
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
   
   // Description of measure and parameters
   virtual std::string description() const = 0;
   
   // In derived classes, this should return a copy of the corresponding
   // derived class
   virtual MeasureDefinition* copy() const = 0;
   
   //Return the MeasureFunction corresponding to this definition
   virtual SharedPtr<MeasureFunction> createMeasureFunction() const = 0;
   
   //Return the AxesFinder that should be used for one-pass minimization
   virtual SharedPtr<AxesFinder> createOnePassAxesFinder() const = 0;

   //Specify whether multi-pass minimization makes sense, and if so, return the AxesFinder that should be used for multi-pass minimization
   virtual bool supportsMultiPassMinimization() const = 0;
   virtual SharedPtr<AxesFinder> createMultiPassAxesFinder( unsigned int) const = 0;
   
};

// The normalized measure, with two parameters: beta and R0
class NormalizedMeasure : public MeasureDefinition {
   
public:
   NormalizedMeasure(double beta, double R0)
   : _beta(beta), _R0(R0) {}
   
   virtual std::string description() const;

   virtual SharedPtr<MeasureFunction> createMeasureFunction() const {
      return SharedPtr<MeasureFunction>(new DefaultNormalizedMeasureFunction(_beta,_R0,std::numeric_limits<double>::max()));
   }

   virtual NormalizedMeasure* copy() const {return new NormalizedMeasure(*this);}
   
   virtual SharedPtr<AxesFinder> createOnePassAxesFinder() const {
      return SharedPtr<AxesFinder>(new AxesFinderFromOnePassMinimization(_beta, std::numeric_limits<double>::max()));
   }
   virtual SharedPtr<AxesFinder> createMultiPassAxesFinder( unsigned int Npass) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromKmeansMinimization(_beta, std::numeric_limits<double>::max(),Npass));
   }

   virtual bool supportsMultiPassMinimization() const { return true; }
   
private:
   double _beta;
   double _R0;
};
   
// The unnormalized measure, with just one parameter: beta
class UnnormalizedMeasure : public MeasureDefinition {

public:
   UnnormalizedMeasure(double beta)
   : _beta(beta) {}

   virtual UnnormalizedMeasure* copy() const {return new UnnormalizedMeasure(*this);}
   
   virtual std::string description() const;
   
   virtual SharedPtr<MeasureFunction> createMeasureFunction() const {
      return SharedPtr<MeasureFunction>(new DefaultUnnormalizedMeasureFunction(_beta,std::numeric_limits<double>::max()));
   }

   
   virtual SharedPtr<AxesFinder> createOnePassAxesFinder() const {
      return SharedPtr<AxesFinder>(new AxesFinderFromOnePassMinimization(_beta, std::numeric_limits<double>::max()));
   }
   virtual SharedPtr<AxesFinder> createMultiPassAxesFinder( unsigned int Npass) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromKmeansMinimization(_beta, std::numeric_limits<double>::max(),Npass));
   }

   virtual bool supportsMultiPassMinimization() const { return true; }

private:
   double _beta;
};
   
   
// The geometric  measure, with 1 parameter: beta
// This measure is still evolving and shouldn't be used for any critial applications yet
class GeometricMeasure : public MeasureDefinition {
   
public:
   GeometricMeasure(double beta)
   : _beta(beta) {}
   
   virtual GeometricMeasure* copy() const {return new GeometricMeasure(*this);}
   
   virtual std::string description() const;
   
   virtual SharedPtr<MeasureFunction> createMeasureFunction() const {
      return SharedPtr<MeasureFunction>(new GeometricMeasureFunction(_beta,std::numeric_limits<double>::max()));
   }
   
   virtual SharedPtr<AxesFinder> createOnePassAxesFinder() const {
      return SharedPtr<AxesFinder>(new AxesFinderFromGeometricMinimization(_beta, std::numeric_limits<double>::max()));
   }
   virtual SharedPtr<AxesFinder> createMultiPassAxesFinder(unsigned int) const {
      throw Error("GeometricMeasure does not support multi-pass minimization.");
      return SharedPtr<AxesFinder>();
   }
   
   virtual bool supportsMultiPassMinimization() const { return false; }

private:
   double _beta;

};


// The normalized cutoff measure, with 3 parameters: beta, R0, Rcutoff
class NormalizedCutoffMeasure : public MeasureDefinition {

public:
   NormalizedCutoffMeasure(double beta, double R0, double Rcutoff)
   : _beta(beta), _R0(R0), _Rcutoff(Rcutoff) {}

   virtual std::string description() const;
   
   virtual NormalizedCutoffMeasure* copy() const {return new NormalizedCutoffMeasure(*this);}
   
   virtual SharedPtr<MeasureFunction> createMeasureFunction() const {
      return SharedPtr<MeasureFunction>(new DefaultNormalizedMeasureFunction(_beta,_R0,_Rcutoff));
   }
   
   virtual SharedPtr<AxesFinder> createOnePassAxesFinder() const {
      return SharedPtr<AxesFinder>(new AxesFinderFromOnePassMinimization(_beta, _Rcutoff));
   }
   virtual SharedPtr<AxesFinder> createMultiPassAxesFinder( unsigned int Npass) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromKmeansMinimization(_beta, std::numeric_limits<double>::max(),Npass));
   }
  
   virtual bool supportsMultiPassMinimization() const { return true; }

private:
   double _beta;
   double _R0;
   double _Rcutoff;
   
};

// The unnormalized cutoff measure, with 2 parameters: beta, Rcutoff
class UnnormalizedCutoffMeasure : public MeasureDefinition {
   
public:
   UnnormalizedCutoffMeasure(double beta, double Rcutoff)
   : _beta(beta), _Rcutoff(Rcutoff) {}

   virtual std::string description() const;
   
   virtual UnnormalizedCutoffMeasure* copy() const {return new UnnormalizedCutoffMeasure(*this);}

   virtual SharedPtr<MeasureFunction> createMeasureFunction() const {
      return SharedPtr<MeasureFunction>(new DefaultUnnormalizedMeasureFunction(_beta,_Rcutoff));
   }
   
   virtual SharedPtr<AxesFinder> createOnePassAxesFinder() const {
      return SharedPtr<AxesFinder>(new AxesFinderFromOnePassMinimization(_beta, _Rcutoff));
   }
   virtual SharedPtr<AxesFinder> createMultiPassAxesFinder( unsigned int Npass) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromKmeansMinimization(_beta, std::numeric_limits<double>::max(),Npass));
   }

   virtual bool supportsMultiPassMinimization() const { return true; }

private:
   double _beta;
   double _Rcutoff;
   
};

// The Geometric  measure, with 2 parameters: beta, Rcutoff
// This measure is still evolving and shouldn't be used for any critial applications yet
class GeometricCutoffMeasure : public MeasureDefinition {
   
public:
   GeometricCutoffMeasure(double beta, double Rcutoff)
   : _beta(beta), _Rcutoff(Rcutoff) {}

   virtual GeometricCutoffMeasure* copy() const {return new GeometricCutoffMeasure(*this);}

   virtual std::string description() const;
   
   virtual SharedPtr<MeasureFunction> createMeasureFunction() const {
      return SharedPtr<MeasureFunction>(new GeometricMeasureFunction(_beta,_Rcutoff));
   }
   
   virtual SharedPtr<AxesFinder> createOnePassAxesFinder() const {
      return SharedPtr<AxesFinder>(new AxesFinderFromGeometricMinimization(_beta, _Rcutoff));
   }
   virtual SharedPtr<AxesFinder> createMultiPassAxesFinder(unsigned int) const {
      throw Error("GeometricCutoffMeasure does not support multi-pass minimization.");
      return SharedPtr<AxesFinder>();
   }
   
   virtual bool supportsMultiPassMinimization() const { return false; }

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
   
   virtual ~AxesDefinition() {};
   
   // description of axes (and any parameters)
   virtual std::string short_description() const = 0;
   virtual std::string description() const = 0;
   
   // In derived classes, this should return a copy of the corresponding
   // derived class
   virtual AxesDefinition* copy() const = 0;
   
   // These describe how the axes finder works
   virtual bool givesRandomizedResults() const = 0;
   virtual bool supportsManualAxes() const = 0;
   
   //for best testing
   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition &) const = 0;
   virtual SharedPtr<AxesFinder> createFinishingAxesFinder(const MeasureDefinition &) const {
      return SharedPtr<AxesFinder>();  //By default, nothing.
   };

};

// kt axes
class KT_Axes : public AxesDefinition {
public:
   KT_Axes() {}

   virtual std::string short_description() const {
      return "KT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "KT Axes";
      return stream.str();
   };
   
   virtual KT_Axes* copy() const {return new KT_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition &) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromKT());
   }
   
};

// ca axes
class CA_Axes : public AxesDefinition {
public:
   CA_Axes() {}

   virtual std::string short_description() const {
      return "CA";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "CA Axes";
      return stream.str();
   };
   
   virtual CA_Axes* copy() const {return new CA_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition &) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromCA());
   }
   
};

// anti-kt axes, one parameter R0 is subjet radius
class AntiKT_Axes : public AxesDefinition {

public:
   AntiKT_Axes(double R0 = 0.2): _R0(R0) {}

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
   
   virtual AntiKT_Axes* copy() const {return new AntiKT_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition &) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromAntiKT(_R0));
   }
   
   
private:
   double _R0;

};

// winner-take-all recombination with kt axes
class WTA_KT_Axes : public AxesDefinition {
public:
   WTA_KT_Axes() {}

   virtual std::string short_description() const {
      return "WTA KT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Winner-Take-All KT Axes";
      return stream.str();
   };
   
   virtual WTA_KT_Axes* copy() const {return new WTA_KT_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition &) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromWTA_KT());
   }

   
};

// winner-take-all recombination with CA axes
class WTA_CA_Axes : public AxesDefinition {
public:
   WTA_CA_Axes() {}

   virtual std::string short_description() const {
      return "WTA CA";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Winner-Take-All CA Axes";
      return stream.str();
   };
   
   virtual WTA_CA_Axes* copy() const {return new WTA_CA_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition &) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromWTA_CA());
   }
   
};
   
// Onepass minimization from kt axes
class OnePass_KT_Axes : public AxesDefinition {
public:
   OnePass_KT_Axes() {}
   
   virtual std::string short_description() const {
      return "OnePass KT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from KT Axes";
      return stream.str();
   };
   
   virtual OnePass_KT_Axes* copy() const {return new OnePass_KT_Axes(*this);}
   
   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition & ) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromKT());
   }
   virtual SharedPtr<AxesFinder> createFinishingAxesFinder(const MeasureDefinition & measure_def) const {
      return measure_def.createOnePassAxesFinder();
   }

};

// Onepass minimization from CA axes
class OnePass_CA_Axes : public AxesDefinition {
public:
   OnePass_CA_Axes() {}

   virtual std::string short_description() const {
      return "OnePass CA";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from CA Axes";
      return stream.str();
   };
   
   virtual OnePass_CA_Axes* copy() const {return new OnePass_CA_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition & ) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromCA());
   }
   virtual SharedPtr<AxesFinder> createFinishingAxesFinder(const MeasureDefinition & measure_def) const {
      return measure_def.createOnePassAxesFinder();
   }
};

// Onepass minimization from AntiKT axes, one parameter R0
class OnePass_AntiKT_Axes : public AxesDefinition {

public:
   OnePass_AntiKT_Axes(double R0): _R0(R0) {}
   
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
   
   virtual OnePass_AntiKT_Axes* copy() const {return new OnePass_AntiKT_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition & ) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromAntiKT(_R0));
   }
   virtual SharedPtr<AxesFinder> createFinishingAxesFinder(const MeasureDefinition & measure_def) const {
      return measure_def.createOnePassAxesFinder();
   }

private:
   double _R0;
};

// Onepass minimization from winner-take-all kt axes
class OnePass_WTA_KT_Axes : public AxesDefinition {
public:
   OnePass_WTA_KT_Axes() {}
   
   virtual std::string short_description() const {
      return "OnePass WTA KT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Winner-Take-All KT Axes";
      return stream.str();
   };
   
   virtual OnePass_WTA_KT_Axes* copy() const {return new OnePass_WTA_KT_Axes(*this);}
   
   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition & ) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromWTA_KT());
   }
   virtual SharedPtr<AxesFinder> createFinishingAxesFinder(const MeasureDefinition & measure_def) const {
      return measure_def.createOnePassAxesFinder();
   }
};


// Onepass minimization from winner-take-all CA axes
class OnePass_WTA_CA_Axes : public AxesDefinition {
public:
   OnePass_WTA_CA_Axes() {}

   virtual std::string short_description() const {
      return "OnePass WTA CA";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Winner-Take-All CA Axes";
      return stream.str();
   };
   
   virtual OnePass_WTA_CA_Axes* copy() const {return new OnePass_WTA_CA_Axes(*this);}
   

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return false;}

   
   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition & ) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromWTA_CA());
   }
   virtual SharedPtr<AxesFinder> createFinishingAxesFinder(const MeasureDefinition & measure_def) const {
      return measure_def.createOnePassAxesFinder();
   }

   
};
   
// set axes manually
class Manual_Axes : public AxesDefinition {
public:
   Manual_Axes() {}
   
   virtual std::string short_description() const {
      return "Manual";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Manual Axes";
      return stream.str();
   };
   
   virtual Manual_Axes* copy() const {return new Manual_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return true;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition &) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromUserInput());
   }
};

// one pass minimization from manual starting point
class OnePass_Manual_Axes : public AxesDefinition {
public:
   OnePass_Manual_Axes() {}

   virtual std::string short_description() const {
      return "OnePass Manual";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Manual Axes";
      return stream.str();
   };
   
   virtual OnePass_Manual_Axes* copy() const {return new OnePass_Manual_Axes(*this);}

   virtual bool givesRandomizedResults() const {return false;}
   virtual bool supportsManualAxes() const {return true;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition & ) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromUserInput());
   }
   virtual SharedPtr<AxesFinder> createFinishingAxesFinder(const MeasureDefinition & measure_def) const {
      return measure_def.createOnePassAxesFinder();
   }

   
};
   
// multi-pass minimization from kT starting point
class MultiPass_Axes : public AxesDefinition {

public:
   MultiPass_Axes(unsigned int Npass) : _Npass(Npass) {}

   virtual std::string short_description() const {
      return "MultiPass";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Multi-Pass Axes (Npass = " << _Npass << ")";
      return stream.str();
   };
   
   virtual MultiPass_Axes* copy() const {return new MultiPass_Axes(*this);}

   virtual bool givesRandomizedResults() const {return true;}
   virtual bool supportsManualAxes() const {return false;}

   virtual SharedPtr<AxesFinder> createStartingAxesFinder(const MeasureDefinition & ) const {
      return SharedPtr<AxesFinder>(new AxesFinderFromKT());
   }
   virtual SharedPtr<AxesFinder> createFnishingAxesFinder(const MeasureDefinition & measure_def) const {
      return measure_def.createMultiPassAxesFinder(_Npass);
   }
   
private:
   unsigned int _Npass;
   
};
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

