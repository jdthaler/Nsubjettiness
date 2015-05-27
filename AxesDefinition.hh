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
#include "WinnerTakeAllRecombiner.hh"

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
class AntiKT_Axes;         // (R0)
class WTA_KT_Axes;
class WTA_CA_Axes;
class WTA_GenKT_Axes;      // (p, R0)
class GenET_GenKT_Axes;    // (delta, p, R0)
class Manual_Axes;
   
class OnePass_KT_Axes;
class OnePass_CA_Axes;
class OnePass_AntiKT_Axes;       // (R0)
class OnePass_WTA_KT_Axes;
class OnePass_WTA_CA_Axes;
class OnePass_WTA_GenKT_Axes;    // (p, R0)
class OnePass_GenET_GenKT_Axes;  // (delta, p, R0)
class OnePass_Manual_Axes;
   
class MultiPass_Axes;            // (NPass) (currently only defined for KT_Axes)
class MultiPass_Manual_Axes;     // (NPass)

class Comb_WTA_GenKT_Axes;       // (p, R0, nExtra)
class Comb_GenET_GenKT_Axes;     // (delta, p, R0, nExtra)

///////
//
// AxesDefinition
//
///////

//------------------------------------------------------------------------
/// \class AxesDefinition
// This is the base class for axes definitions.  A generic AxesDefinition
// first finds a set of seed axes, then (if desired) uses measure information
// (from MeasureFunction) to refine those axes starting from those seed axes.
// The AxesDefinitions are typically based on sequential jet algorithms.
class AxesDefinition {
   
public:
   
   // This function should be overloaded in all derived classes, and defines how to find the seed axes
   // If desired, the measure information (which might be NULL) can be used to test multiple axes choices, but should
   // not be used for iterative refining (since that is the job of MeasureFunction).
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector<fastjet::PseudoJet>& inputs,
                                                             const MeasureDefinition * measure) const = 0;
   
   
   
   
   // description of AxesDefinitions (and any parameters)
   virtual std::string short_description() const = 0;
   virtual std::string description() const = 0;
   
   // This has to be defined in all derived classes
   // Allows these to be copied around.
   virtual AxesDefinition* create() const = 0;
   
public:
   
   // Starting from seeds, refine axes using one or more passes.
   // Note that in order to do >0 passes,
   // we need information from the measure about how to do the appropriate minimization.
   std::vector<fastjet::PseudoJet> get_refined_axes(int n_jets,
                                                    const std::vector<fastjet::PseudoJet>& inputs,
                                                    const std::vector<fastjet::PseudoJet>& seedAxes,
                                                    const MeasureDefinition * measure = NULL) const {

      assert(n_jets == (int)seedAxes.size()); //added int casting to get rid of compiler warning
      
      // if (_needsManualAxes) throw Error("AxesDefinition:  You can't get_axes in Manual Mode");
      
      if (_Npass == 0) {
         // no refining, just use seeds
         return seedAxes;
      } else if (_Npass == 1) {
         if (measure == NULL) throw Error("AxesDefinition:  One-pass minimization requires specifying a MeasureDefinition.");
         
         // do one pass minimum using measure definition
         return measure->get_one_pass_axes(n_jets, inputs, seedAxes,_nAttempts,_accuracy);
      } else {
         if (measure == NULL) throw Error("AxesDefinition:  Minimization requires specifying a MeasureDefinition.");
         return get_multi_pass_axes(n_jets, inputs, seedAxes, measure);
      }
   }
   
   // Combining get_starting_axes with get_refined_axes.
   // In Njettiness class, these two steps are done separately in order to store
   // seed axes information.
   std::vector<fastjet::PseudoJet> get_axes(int n_jets,
                                            const std::vector<fastjet::PseudoJet>& inputs,
                                            const MeasureDefinition * measure = NULL) const {
      std::vector<fastjet::PseudoJet> seedAxes = get_starting_axes(n_jets, inputs, measure);
      return get_refined_axes(n_jets,inputs,seedAxes,measure);
   }

   
   // short-hand for the get_axes function.  Useful when trying to write optimized code.
   inline std::vector<fastjet::PseudoJet> operator() (int n_jets,
                                               const std::vector<fastjet::PseudoJet>& inputs,
                                               const MeasureDefinition * measure = NULL) const {
      return get_axes(n_jets,inputs,measure);
   }
   
   // define the cases of zero pass and one pass for convenience
   enum AxesRefiningEnum {
      UNDEFINED_REFINE = -1, // added to create a default value
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
   
   // Allow user to change number of passes.  Also used internally to set nPass
   void setNPass(int nPass,
                 int nAttempts = 1000,
                 double accuracy  = 0.0001,
                 double noise_range = 1.0 // only needed for MultiPass minimization
                 )
   {
      _Npass = nPass;
      _nAttempts = nAttempts;
      _accuracy = accuracy;
      _noise_range = noise_range;
      if (nPass < 0) throw Error("AxesDefinition requires a nPass >= 0");
   }
   
   // destructor
   virtual ~AxesDefinition() {};
   
protected:
   
   // Default constructor contains no information.  Number of passes has to be set
   // manually by derived classes using setNPass function
   AxesDefinition() : _Npass(UNDEFINED_REFINE),
                     _nAttempts(0),
                     _accuracy(0.0),
                     _noise_range(0.0),
                     _needsManualAxes(false) {}

   // Does multi-pass minimization by jiggling the axes.
   std::vector<fastjet::PseudoJet> get_multi_pass_axes(int n_jets,
                                                       const std::vector<fastjet::PseudoJet>& inputs,
                                                       const std::vector<fastjet::PseudoJet>& seedAxes,
                                                       const MeasureDefinition* measure) const;
   
   PseudoJet jiggle(const PseudoJet& axis) const;
   
   int _Npass; // number of passes (0 = no refining, 1 = one-pass, >1 multi-pass)
   int _nAttempts; // number of attempts per pass
   double _accuracy;  // accuracy goal per pass
   double _noise_range; // defaults to 1.0 (noise in rapidity/phi)
   bool _needsManualAxes; // special case of manual axes
};
  
//------------------------------------------------------------------------
/// \class ExclusiveJetAxes
// This class finds axes by clustering particles with an exclusive jet definition.
// This can be implemented with different jet algorithms, and can be called by the user
class ExclusiveJetAxes : public AxesDefinition {
   
public:
   ExclusiveJetAxes(fastjet::JetDefinition def)
   : AxesDefinition(), _def(def) {
      setNPass(NO_REFINING);    // default to no minimization
   }
   
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector <fastjet::PseudoJet> & inputs,
                                                             const MeasureDefinition * ) const {
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
/// \class ExclusiveCombinatorialJetAxes
// This class finds axes by clustering particles with an exclusive jet definition.
// It takes an extra number of jets (specificed by the user), and then finds the set of N that minimizes N-jettiness
// This can be implemented with different jet algorithms, and can be called by the user
// WARNING: If one wants to be guarenteed that results improve by increasing nExtra, then one should use
// WTA Recombination schemes
class ExclusiveCombinatorialJetAxes : public AxesDefinition {
   
public:
   ExclusiveCombinatorialJetAxes(fastjet::JetDefinition def, int nExtra = 0)
   : AxesDefinition(), _def(def), _nExtra(nExtra) {
      if (nExtra < 0) throw Error("Need nExtra >= 0");
      setNPass(NO_REFINING);   // default to no minimization
   }
   
    virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets, 
                                                           const std::vector<fastjet::PseudoJet> & inputs,
                                                           const MeasureDefinition *measure) const {
      int starting_number = n_jets + _nExtra;
      fastjet::ClusterSequence jet_clust_seq(inputs, _def);
      std::vector<fastjet::PseudoJet> starting_axes = jet_clust_seq.exclusive_jets(starting_number);

      std::vector<fastjet::PseudoJet> final_axes;

      // check so that no computation time is wasted if there are no extra axes
      if (_nExtra == 0) final_axes = starting_axes;

      else {

        // define string of 1's based on number of desired jets
        std::string bitmask(n_jets, 1);
        // expand the array size to the total number of jets with extra 0's at the end, makes string easy to permute
        bitmask.resize(starting_number, 0); 

        double min_tau = std::numeric_limits<double>::max();
        do {
          std::vector<fastjet::PseudoJet> temp_axes;

          // only take an axis if it is listed as true (1) in the string
          for (int i = 0; i < (int)starting_axes.size(); ++i) {
            if (bitmask[i]) temp_axes.push_back(starting_axes[i]);
          }

          double temp_tau = measure->result(inputs, temp_axes);
          if (temp_tau < min_tau) {
            min_tau = temp_tau;
            final_axes = temp_axes;
          }

          // permutes string of 1's and 0's according to next lexicographic ordering and returns true
          // continues to loop through all possible lexicographic orderings
          // returns false and breaks the loop when there are no more possible orderings
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      return final_axes;
    }

   virtual std::string short_description() const { return "ExclCombAxes";}
   virtual std::string description() const { return "ExclCombAxes: " + _def.description();}
   
   virtual ExclusiveCombinatorialJetAxes* create() const {return new ExclusiveCombinatorialJetAxes(*this);}

private:
   fastjet::JetDefinition _def;
   int _nExtra;
};
   
//------------------------------------------------------------------------
/// \class HardestJetAxes
// This class finds axes by running an inclusive algorithm and then finding the n hardest jets.
// This can be implemented with different jet algorithms, and can be called by the user
class HardestJetAxes : public AxesDefinition {
public:
   HardestJetAxes(fastjet::JetDefinition def)
   : AxesDefinition(), _def(def) {
      setNPass(NO_REFINING);    // default to no minimization
   }
   
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector <fastjet::PseudoJet> & inputs,
                                                             const MeasureDefinition * ) const {
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
   KT_Axes()
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::kt_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             fastjet::E_scheme,
                                             fastjet::Best)
                      ) {
      setNPass(NO_REFINING);
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

//------------------------------------------------------------------------
/// \class CA_Axes
// Axes from CA algorithm with E_scheme recombination.
class CA_Axes : public ExclusiveJetAxes {
public:
   CA_Axes()
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::cambridge_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             fastjet::E_scheme,
                                             fastjet::Best)
                      ) {
      setNPass(NO_REFINING);
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

   
//------------------------------------------------------------------------
/// \class AntiKT_Axes
// Axes from anti-kT algorithm and E_scheme.
// The one parameter R0 is subjet radius
class AntiKT_Axes : public HardestJetAxes {

public:
   AntiKT_Axes(double R0)
   : HardestJetAxes(fastjet::JetDefinition(fastjet::antikt_algorithm,
                                           R0,
                                           fastjet::E_scheme,
                                           fastjet::Best)
                    ), _R0(R0) {
      setNPass(NO_REFINING);
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
   
protected:
   double _R0;

};

// Added Jet Definition wrapper class to avoid issue of genKT FastJet bug
// Now using this for all AxesDefinition with a manual recombiner to use the delete_recombiner_when_unused function
class JetDefinitionWrapper {

public: 
   
   JetDefinitionWrapper(JetAlgorithm jet_algorithm_in, double R_in, double xtra_param_in, const JetDefinition::Recombiner *recombiner) {
      jet_def = fastjet::JetDefinition(jet_algorithm_in, R_in, xtra_param_in);
      jet_def.set_recombiner(recombiner);
      jet_def.delete_recombiner_when_unused();        // added to prevent memory leaks
   }

   //additional constructor so that normal jet algorithms can also be called
   JetDefinitionWrapper(JetAlgorithm jet_algorithm_in, double R_in, const JetDefinition::Recombiner *recombiner, fastjet::Strategy strategy_in) {
      jet_def = fastjet::JetDefinition(jet_algorithm_in, R_in, recombiner, strategy_in);
      jet_def.delete_recombiner_when_unused();
   }

   JetDefinition getJetDef() {
      return jet_def;
   }

private:
   JetDefinition jet_def;
};

//------------------------------------------------------------------------
/// \class WTA_KT_Axes
// Axes from kT algorithm and winner-take-all recombination
class WTA_KT_Axes : public ExclusiveJetAxes {
public:
   WTA_KT_Axes()
   : ExclusiveJetAxes(JetDefinitionWrapper(fastjet::kt_algorithm,
                                          fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                          _recomb = new WinnerTakeAllRecombiner(), // Needs to be explicitly declared
                                          fastjet::Best).getJetDef()
                      ) {
      setNPass(NO_REFINING);
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

private:
   const WinnerTakeAllRecombiner *_recomb; 

};
   
//------------------------------------------------------------------------
/// \class WTA_CA_Axes
// Axes from CA algorithm and winner-take-all recombination
class WTA_CA_Axes : public ExclusiveJetAxes {
public:
   WTA_CA_Axes()
   : ExclusiveJetAxes(JetDefinitionWrapper(fastjet::cambridge_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             _recomb = new WinnerTakeAllRecombiner(),
                                             fastjet::Best).getJetDef()) {
    setNPass(NO_REFINING);
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
   
private:
   const WinnerTakeAllRecombiner *_recomb;

};


//------------------------------------------------------------------------
/// \class WTA_GenKT_Axes
// Axes from a general KT algorithm with a Winner Take All Recombiner 
// Requires the power of the KT algorithm to be used and the radius parameter
class WTA_GenKT_Axes : public ExclusiveJetAxes {

public:
   WTA_GenKT_Axes(double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveJetAxes(JetDefinitionWrapper(fastjet::genkt_algorithm,
                                            R0,
                                            p,
                                            _recomb = new WinnerTakeAllRecombiner()
                                            ).getJetDef()), _p(p), _R0(R0) {
        setNPass(NO_REFINING);
    }

   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "WTA, GenKT Axes";
      return stream.str();
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "General KT (p = " << _p << "), Winner-Take-All Recombiner, R0 = " << _R0;
      return stream.str();
   };
   
   virtual WTA_GenKT_Axes* create() const {return new WTA_GenKT_Axes(*this);}
   
protected:
   double _p;
   double _R0;
   const WinnerTakeAllRecombiner *_recomb; 
};
   
//------------------------------------------------------------------------
/// \class GenET_GenKT_Axes
// Class using general KT algorithm with a more general recombination scheme
// Requires power of KT algorithm, power of recombination weights, and radius parameter
class GenET_GenKT_Axes : public ExclusiveJetAxes {

public:
   GenET_GenKT_Axes(double delta, double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveJetAxes((JetDefinitionWrapper(fastjet::genkt_algorithm, R0, p, _recomb = new GeneralEtSchemeRecombiner(delta))).getJetDef() ),
    _delta(delta), _p(p), _R0(R0) {
        setNPass(NO_REFINING);
    }

   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "GenET, GenKT Axes" ;
      return stream.str();
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      // TODO: if _delta is huge, change to "WTA"
      << "General Recombiner (delta = " << _delta << "), " << "General KT (p = " << _p << ") Axes, R0 = " << _R0;
      return stream.str();
   };
   
   virtual GenET_GenKT_Axes* create() const {return new GenET_GenKT_Axes(*this);}
   
protected:
   double _delta;
   double _p;
   double _R0;
   const GeneralEtSchemeRecombiner *_recomb;
};

//------------------------------------------------------------------------
/// \class OnePass_KT_Axes
// Onepass minimization from kt axes
class OnePass_KT_Axes : public KT_Axes {
public:
   OnePass_KT_Axes() : KT_Axes() {
      setNPass(ONE_PASS);
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

//------------------------------------------------------------------------
/// \class OnePass_CA_Axes
// Onepass minimization from CA axes
class OnePass_CA_Axes : public CA_Axes {
public:
   OnePass_CA_Axes() : CA_Axes() {
      setNPass(ONE_PASS);
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
   
//------------------------------------------------------------------------
/// \class OnePass_AntiKT_Axes
// Onepass minimization from AntiKT axes, one parameter R0
class OnePass_AntiKT_Axes : public AntiKT_Axes {

public:
   OnePass_AntiKT_Axes(double R0) : AntiKT_Axes(R0) {
      setNPass(ONE_PASS);
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

};

//------------------------------------------------------------------------
/// \class OnePass_WTA_KT_Axes
// Onepass minimization from winner-take-all kt axes
class OnePass_WTA_KT_Axes : public WTA_KT_Axes {
public:
   OnePass_WTA_KT_Axes() : WTA_KT_Axes() {
      setNPass(ONE_PASS);
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

//------------------------------------------------------------------------
/// \class OnePass_WTA_CA_Axes
// Onepass minimization from winner-take-all CA axes
class OnePass_WTA_CA_Axes : public WTA_CA_Axes {
   
public:
   OnePass_WTA_CA_Axes() : WTA_CA_Axes() {
      setNPass(ONE_PASS);
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

// ------------------------------------------------------------------------
/// \class OnePass_WTA_GenKT_Axes
// Onepass minimization from winner-take-all, General KT Axes 
class OnePass_WTA_GenKT_Axes : public WTA_GenKT_Axes {
   
public:
   OnePass_WTA_GenKT_Axes(double p, double R0 = fastjet::JetDefinition::max_allowable_R) : WTA_GenKT_Axes(p, R0) {
      setNPass(ONE_PASS);
   }

   virtual std::string short_description() const {
      return "OnePass WTA GenKT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Winner-Take-All, General KT";
      return stream.str();
   };
   
   virtual OnePass_WTA_GenKT_Axes* create() const {return new OnePass_WTA_GenKT_Axes(*this);}
};

// ------------------------------------------------------------------------
/// \class OnePass_GenET_GenKT_Axes
// Onepass minimization from General Recomb, General KT axes
class OnePass_GenET_GenKT_Axes : public GenET_GenKT_Axes {
   
public:
   OnePass_GenET_GenKT_Axes(double delta, double p, double R0 = fastjet::JetDefinition::max_allowable_R) : GenET_GenKT_Axes(delta, p, R0) {
      setNPass(ONE_PASS);
   }

   virtual std::string short_description() const {
      return "OnePass GenET, GenKT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from General Recomb, General KT";
      return stream.str();
   };
   
   virtual OnePass_GenET_GenKT_Axes* create() const {return new OnePass_GenET_GenKT_Axes(*this);}
};


//------------------------------------------------------------------------
/// \class Manual_Axes
// set axes manually
class Manual_Axes : public AxesDefinition {
public:
   Manual_Axes() : AxesDefinition() {
      setNPass(NO_REFINING);
      _needsManualAxes = true;
   }
   
   // dummy function since this is manual mode
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int,
                                                             const std::vector<fastjet::PseudoJet>&,
                                                             const MeasureDefinition *) const;

   
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
   OnePass_Manual_Axes() : Manual_Axes() {
      setNPass(ONE_PASS);
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
   
//------------------------------------------------------------------------
/// \class MultiPass_Axes
// multi-pass minimization from kT starting point
class MultiPass_Axes : public KT_Axes {

public:
   MultiPass_Axes(unsigned int Npass) : KT_Axes() {
      setNPass(Npass);
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
   
};

//------------------------------------------------------------------------
/// \class MultiPass_Manual_Axes
// multi-pass minimization from kT starting point
class MultiPass_Manual_Axes : public Manual_Axes {

public:
   MultiPass_Manual_Axes(unsigned int Npass) : Manual_Axes() {
      setNPass(Npass);
   }

   virtual std::string short_description() const {
      return "MultiPass Manual";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Multi-Pass Manual Axes (Npass = " << _Npass << ")";
      return stream.str();
   };
   
   virtual MultiPass_Manual_Axes* create() const {return new MultiPass_Manual_Axes(*this);}
   
};


//------------------------------------------------------------------------
/// \class Comb_WTA_KT_Axes
// Axes from kT algorithm and winner-take-all recombination
// Requires nExtra parameter and returns set of N that minimizes N-jettiness

class Comb_WTA_GenKT_Axes : public ExclusiveCombinatorialJetAxes {
public:
   Comb_WTA_GenKT_Axes(int nExtra, double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveCombinatorialJetAxes((JetDefinitionWrapper(fastjet::genkt_algorithm, R0, p, _recomb = new WinnerTakeAllRecombiner())).getJetDef(), nExtra),
    _p(p), _R0(R0) {
        setNPass(NO_REFINING);
    }

   virtual std::string short_description() const {
      return "N Choose M WTA GenKT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "N Choose M General KT (p = " << _p << "), Winner-Take-All Recombiner, R0 = " << _R0;
      return stream.str();
   };
   
   virtual Comb_WTA_GenKT_Axes* create() const {return new Comb_WTA_GenKT_Axes(*this);}

private:
   double _p;
   double _R0;
   const WinnerTakeAllRecombiner *_recomb; 
};
   
//------------------------------------------------------------------------
/// \class Comb_GenET_KT_Axes
// Axes from kT algorithm and General Et scheme recombination
// Requires nExtra parameter and returns set of N that minimizes N-jettiness
class Comb_GenET_GenKT_Axes : public ExclusiveCombinatorialJetAxes {
public:
   Comb_GenET_GenKT_Axes(int nExtra, double delta, double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveCombinatorialJetAxes((JetDefinitionWrapper(fastjet::genkt_algorithm, R0, p, _recomb = new GeneralEtSchemeRecombiner(delta))).getJetDef(), nExtra),
    _delta(delta), _p(p), _R0(R0) {
        setNPass(NO_REFINING);
    }

   virtual std::string short_description() const {
      return "N Choose M WTA GenKT";
   };
   
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "N Choose M General KT (p = " << _p << "), General Et-Scheme Recombiner (delta = " << _delta << "), R0 = " << _R0;
      return stream.str();
   };
   
   virtual Comb_GenET_GenKT_Axes* create() const {return new Comb_GenET_GenKT_Axes(*this);}

private:
   double _delta;
   double _p;
   double _R0;
   const GeneralEtSchemeRecombiner *_recomb; 
};
   

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

