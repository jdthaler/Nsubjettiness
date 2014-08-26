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

#ifndef __FASTJET_CONTRIB_AXES_DEFINITION_HH__
#define __FASTJET_CONTRIB_AXES_DEFINITION_HH__


#include "MeasureDefinition.hh"
#include "AxesRefiner.hh"

#include "fastjet/PseudoJet.hh"
#include <fastjet/LimitedWarning.hh>

#include <iomanip>
#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
   
// The following AxesDefinitions are currently available (and the relevant arguments, if needed)
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
class MultiPass_Axes;  // (NPass)
  
   
///////
//
// AxesDefinition
//
///////

//------------------------------------------------------------------------
/// \class AxesDefinition
// This is the base class for axes definitions.  Note that there is a difference
// between an AxesDefinition and an AxesRefiner.  An AxesDefinition finds axes
// without requiring seeds, whereas an AxesRefiner requires seeds.
// The AxesDefinitions are typically based on sequential jet algorithms.
// At the moment, most AxesDefinition do not have an arugment (except the Anti-KT ones)
class AxesDefinition {
   
public:
   
   // This function should be overloaded in all derived classes, and defines how to find the starting axes
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector<fastjet::PseudoJet>& inputs) const = 0;
   
   // description of AxesDefinitions (and any parameters)
   virtual std::string short_description() const = 0;
   virtual std::string description() const = 0;
   
   // This has to be defined in all derived classes
   // Allows these to be copied around.
   virtual AxesDefinition* create() const = 0;
   
public:
   
   // Finding axes, including any required minimziation.  Note that in order to do >0 passes,
   // we need information from the measure about how to do the appropriate minimization.
   std::vector<fastjet::PseudoJet> get_axes(int n_jets,
                                            const std::vector<fastjet::PseudoJet>& inputs,
                                            const MeasureDefinition * measure = NULL) const {
      if (_needsManualAxes) throw Error("AxesDefinition:  You can't get_axes in Manual Mode");
      
      if (_Npass == 0) {
         return get_starting_axes(n_jets,inputs);
      } else {
         if (measure == NULL) throw Error("AxesDefinition:  Minimization requires specifying a MeasureDefinition.");
         SharedPtr<AxesRefiner> refiner(measure->createAxesRefiner(_Npass));
         return refiner->get_axes(n_jets,inputs,get_starting_axes(n_jets,inputs));
      }
   }
   
   // short-hand for the get_axes function.  Useful when trying to write optimized code.
   inline std::vector<fastjet::PseudoJet> operator() (int n_jets,
                                               const std::vector<fastjet::PseudoJet>& inputs,
                                               const MeasureDefinition * measure = NULL) const {
      return get_axes(n_jets,inputs,measure);
   }
   
   // define the cases of zero pass and one pass for convenience
   enum AxesRefiningEnum {
      NO_REFINING = 0,
      ONE_PASS = 1,
      MULTI_PASS = 100,
   };
   
   // This is a flag that is used externally to decide how to do minimization
   int nPass() const { return _Npass; }

   bool givesRandomizedResults() const {
      return (_Npass > 1);
   }
   
   bool needsManualAxes() const {
      return _needsManualAxes; // if there is no starting axes finder
   }
   
   // destructor
   virtual ~AxesDefinition() {};
   
protected:
   
   // Constructor that requires knowing the number of passes
   AxesDefinition(int nPass) : _Npass(nPass), _needsManualAxes(false) {
      if (nPass < 0) throw Error("AxesDefinition requires a nPass >= 0.");
   }
   
   int _Npass; // number of passes (0 = no refining, 1 = one-pass, >1 multi-pass)
   bool _needsManualAxes; // special case of manual axes
};
  

//------------------------------------------------------------------------
/// \class ExclusiveJetAxes
// This class finds axes by clustering particles with an exclusive jet definition.
// This can be implemented with different jet algorithms, and can be called by the user
class ExclusiveJetAxes : public AxesDefinition {
   
public:
   ExclusiveJetAxes(fastjet::JetDefinition def, int nPass = NO_REFINING)  // default to no minimization
   : AxesDefinition(nPass), _def(def) {}
   
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector <fastjet::PseudoJet> & inputs) const {
      fastjet::ClusterSequence jet_clust_seq(inputs, _def);
      return jet_clust_seq.exclusive_jets(n_jets);
   }
   
   virtual std::string short_description() const { return "ExclAxes";}
   virtual std::string description() const { return "ExclAxes: " + _def.description();}
   
   virtual ExclusiveJetAxes* create() const {return new ExclusiveJetAxes(*this);}

private:
   fastjet::JetDefinition _def;
   
};
   
//------------------------------------------------------------------------
/// \class HardestJetAxes
// This class finds axes by running an inclusive algorithm and then finding the n hardest jets.
// This can be implemented with different jet algorithms, and can be called by the user
class HardestJetAxes : public AxesDefinition {
public:
   HardestJetAxes(fastjet::JetDefinition def, int nPass)
   : AxesDefinition(nPass), _def(def) {}
   
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector <fastjet::PseudoJet> & inputs) const {
      fastjet::ClusterSequence jet_clust_seq(inputs, _def);
      std::vector<fastjet::PseudoJet> myJets = sorted_by_pt(jet_clust_seq.inclusive_jets());
      myJets.resize(n_jets);  // only keep n hardest
      return myJets;
   }
   
   virtual std::string short_description() const { return "HardAxes";}
   virtual std::string description() const { return "HardAxes: " + _def.description();}
   
   virtual HardestJetAxes* create() const {return new HardestJetAxes(*this);}
   
private:
   fastjet::JetDefinition _def;
};
   
   
   
//------------------------------------------------------------------------
/// \class KT_Axes
// Axes from kT algorithm with E_scheme recombination.
class KT_Axes : public ExclusiveJetAxes {
public:
   KT_Axes(int nPass = NO_REFINING)
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::kt_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             fastjet::E_scheme,
                                             fastjet::Best),
                      nPass) {}

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

//------------------------------------------------------------------------
/// \class CA_Axes
// Axes from CA algorithm with E_scheme recombination.
class CA_Axes : public ExclusiveJetAxes {
public:
   CA_Axes(int nPass = NO_REFINING)
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::cambridge_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             fastjet::E_scheme,
                                             fastjet::Best),
                      nPass) {}

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

   
//------------------------------------------------------------------------
/// \class AntiKT_Axes
// Axes from anti-kT algorithm and E_scheme.
// The one parameter R0 is subjet radius
class AntiKT_Axes : public HardestJetAxes {

public:
   AntiKT_Axes(double R0, int nPass = NO_REFINING)
   : HardestJetAxes(fastjet::JetDefinition(fastjet::antikt_algorithm,
                                           R0,
                                           fastjet::E_scheme,
                                           fastjet::Best),
                    nPass), _R0(R0) {}

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
   
protected:
   double _R0;

};

//------------------------------------------------------------------------
/// \class WTA_KT_Axes
// Axes from kT algorithm and winner-take-all recombination
class WTA_KT_Axes : public ExclusiveJetAxes {
public:
   WTA_KT_Axes(int nPass = NO_REFINING)
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::kt_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             &_recomb,
                                             fastjet::Best),
                      nPass) {}

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

private:
   const WinnerTakeAllRecombiner _recomb;

   
};
   
//------------------------------------------------------------------------
/// \class WTA_CA_Axes
// Axes from CA algorithm and winner-take-all recombination
class WTA_CA_Axes : public ExclusiveJetAxes {
public:
   WTA_CA_Axes(int nPass = NO_REFINING)
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::cambridge_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             &_recomb,
                                             fastjet::Best),
                      nPass) {}

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
   
private:
   const WinnerTakeAllRecombiner _recomb;


};
   
//------------------------------------------------------------------------
/// \class OnePass_KT_Axes
// Onepass minimization from kt axes
class OnePass_KT_Axes : public KT_Axes {
public:
   OnePass_KT_Axes() : KT_Axes(ONE_PASS) {}
   
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

//------------------------------------------------------------------------
/// \class OnePass_CA_Axes
// Onepass minimization from CA axes
class OnePass_CA_Axes : public CA_Axes {
public:
   OnePass_CA_Axes() : CA_Axes(ONE_PASS) {}

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
   
//------------------------------------------------------------------------
/// \class OnePass_AntiKT_Axes
// Onepass minimization from AntiKT axes, one parameter R0
class OnePass_AntiKT_Axes : public AntiKT_Axes {

public:
   OnePass_AntiKT_Axes(double R0) : AntiKT_Axes(R0, ONE_PASS) {}
   
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

};

//------------------------------------------------------------------------
/// \class OnePass_WTA_KT_Axes
// Onepass minimization from winner-take-all kt axes
class OnePass_WTA_KT_Axes : public WTA_KT_Axes {
public:
   OnePass_WTA_KT_Axes() : WTA_KT_Axes(ONE_PASS) {}
   
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

//------------------------------------------------------------------------
/// \class OnePass_WTA_CA_Axes
// Onepass minimization from winner-take-all CA axes
class OnePass_WTA_CA_Axes : public WTA_CA_Axes {
   
public:
   OnePass_WTA_CA_Axes() : WTA_CA_Axes(ONE_PASS) {}

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
   
//------------------------------------------------------------------------
/// \class Manual_Axes
// set axes manually
class Manual_Axes : public AxesDefinition {
public:
   Manual_Axes(int nPass = NO_REFINING) : AxesDefinition(nPass) {
      _needsManualAxes = true;
   }
   
   // dummy function since this is manual mode
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int,
                                                             const std::vector<fastjet::PseudoJet>&) const;

   
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

//------------------------------------------------------------------------
/// \class OnePass_Manual_Axes
// one pass minimization from manual starting point
class OnePass_Manual_Axes : public Manual_Axes {
public:
   OnePass_Manual_Axes() : Manual_Axes(ONE_PASS) {}

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
   
//------------------------------------------------------------------------
/// \class MultiPass_Axes
// multi-pass minimization from kT starting point
class MultiPass_Axes : public KT_Axes {

public:
   MultiPass_Axes(unsigned int Npass) : KT_Axes(Npass) {}

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
   
};
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

