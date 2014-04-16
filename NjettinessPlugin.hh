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

#ifndef __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__
#define __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__

#include "Njettiness.hh"
#include "MeasureFunction.hh"
#include "AxesFinder.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib {

//------------------------------------------------------------------------
/// \class NjettinessExtras
// This class contains the same information as Njettiness, but redoes it in terms of the ClusterSequence::Extras class.
// This is done in order to help improve the interface for the main NjettinessPlugin class.
// TODO:  This class should probably be merged with TauComponents, since both have access
// to similar information
class NjettinessExtras : public ClusterSequence::Extras {
   private:
   
      TauComponents _tau_components;
      std::vector<fastjet::PseudoJet> _jets;
      std::vector<fastjet::PseudoJet> _axes;
      
      int labelOf(const fastjet::PseudoJet& jet) const {
         int thisJet = -1;
         for (unsigned int i = 0; i < _jets.size(); i++) {
            if (_jets[i].cluster_hist_index() == jet.cluster_hist_index()) {
               thisJet = i;
               break;
            }
         }
         return thisJet;
      }
      
   public:
      NjettinessExtras(TauComponents tau_components, std::vector<fastjet::PseudoJet> jets, std::vector<fastjet::PseudoJet> axes) : _tau_components(tau_components), _jets(jets), _axes(axes) {}
      
      double totalTau() const {return _tau_components.tau();}
      std::vector<double> subTaus() const {return _tau_components.jet_pieces();}
      std::vector<fastjet::PseudoJet> jets() const {return _jets;}
      std::vector<fastjet::PseudoJet> axes() const {return _axes;}
      
      double totalTau(const fastjet::PseudoJet& jet) const {
         return _tau_components.tau();
      }
      double subTau(const fastjet::PseudoJet& jet) const {
         if (labelOf(jet) == -1) return NAN;
         return _tau_components.jet_pieces()[labelOf(jet)];
      }
      
      double beamTau() const {
         return _tau_components.beam_piece();
      }
      
      fastjet::PseudoJet axis(const fastjet::PseudoJet& jet) const {
         return _axes[labelOf(jet)];
      }

      bool has_njettiness_extras(const fastjet::PseudoJet& jet) const {
         return (labelOf(jet) >= 0);
      }

};

inline const NjettinessExtras * njettiness_extras(const fastjet::PseudoJet& jet) {
   const ClusterSequence * myCS = jet.associated_cluster_sequence();   
   if (myCS == NULL) return NULL;
   const NjettinessExtras* extras = dynamic_cast<const NjettinessExtras*>(myCS->extras());   
   return extras;   
}

inline const NjettinessExtras * njettiness_extras(const fastjet::ClusterSequence& myCS) {
   const NjettinessExtras* extras = dynamic_cast<const NjettinessExtras*>(myCS.extras());   
   return extras;   
}

/// The Njettiness jet algorithm
/**
 * An exclusive jet finder that identifies N jets; first N axes are found, then
 * particles are assigned to the nearest (DeltaR) axis and for each axis the
 * corresponding jet is simply the four-momentum sum of these particles.
 *
 * Axes can be found in several ways, specified by the AxesMode argument.  The
 * recommended choices are
 *
 * kt_axes              : exclusive kT
 * wta_kt_axes          : exclusive kT with winner-take-all-recombination
 * onepass_kt_axes      : one-pass minimization seeded by kt (pretty good)
 * onepass_wta_kt_axes  : one-pass minimization seeded by wta_kt
 *
 * For the unnormalized_measure, N-jettiness is defined as:
 *
 * tau_N = Sum_{all particles i} p_T^i min((DR_i1)^beta, (DR_i2)^beta, ...)
 *
 *   DR_ij is the distance sqrt(Delta_phi^2 + Delta_rap^2) between particle i
 *   and jet j.
 * 
 * The normalized_meausure include an extra parameter R0, and the various cutoff
 * measures include an Rcutoff, which effectively defines an angular cutoff
 * similar in effect to a cone-jet radius.
 *
 */

class NjettinessPlugin : public JetDefinition::Plugin {
public:

   NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode,
                    double para1 = NAN,
                    double para2 = NAN,
                    double para3 = NAN,
                    double para4 = NAN);

   // Old constructor for backwards compatibility with v1.0,
   // where normalized_cutoff_measure was the only option
   NjettinessPlugin(int N,
                    Njettiness::AxesMode mode,
                    double beta,
                    double R0,
                    double Rcutoff=std::numeric_limits<double>::max());


   // The things that are required by base class.
   virtual std::string description () const;
   virtual double R() const {return -1.0;} // TODO: make this not stupid
   virtual void run_clustering(ClusterSequence&) const;

   virtual ~NjettinessPlugin() {}

private:

   int _N;
   mutable Njettiness _njettinessFinder; // TODO:  should muck with this so run_clustering can be const without this mutable

};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__