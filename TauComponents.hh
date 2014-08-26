//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: MeasureFunction.hh 742 2014-08-23 15:43:29Z jthaler $
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

#ifndef __FASTJET_CONTRIB_TAUCOMPONENTS_HH__
#define __FASTJET_CONTRIB_TAUCOMPONENTS_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <cmath>
#include <vector>
#include <list>
#include <limits>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

// Classes defined in this file.
class TauExtras;
class TauComponents;
class TauPartition;
class NjettinessExtras;

   
///////
//
// TauComponents
//
///////
   
/// \class TauComponents
// This class creates a wrapper for the various tau/subtau values calculated in Njettiness. This class allows Njettiness access to these variables
// without ever having to do the calculation itself. It takes in subtau numerators and tau denominator from MeasureFunction
// and outputs tau numerator, and normalized tau and subtau.
class TauComponents {
   
public:
   
   // empty constructor necessary to initialize tau_components in Njettiness
   // later set correctly in Njettiness::getTau function
   TauComponents() {
      _jet_pieces_numerator.resize(1, 0.0);
      _beam_piece_numerator = 0.0;
      _denominator = 0;
      _numerator = 0;
      _jet_pieces.resize(1, 0.0);
      _beam_piece = 0.0;
      _tau = 0;
      _has_denominator = false;
      _has_beam = false;
   }
   
   // This constructor takes input vector and double and calculates all necessary tau components
   TauComponents(std::vector<double> jet_pieces_numerator, double beam_piece_numerator, double denominator, bool has_denominator, bool has_beam)
   :  _jet_pieces_numerator(jet_pieces_numerator),
      _beam_piece_numerator(beam_piece_numerator),
      _denominator(denominator),
      _has_denominator(has_denominator),
      _has_beam(has_beam) {
      
      if (!_has_denominator) assert(_denominator == 1.0); //make sure no effect from _denominator if _has_denominator is false
      if (!_has_beam) assert (_beam_piece_numerator == 0.0); //make sure no effect from _beam_piece_numerator if _has_beam is false
      
      _numerator = _beam_piece_numerator;
      _jet_pieces.resize(_jet_pieces_numerator.size(),0.0);
      for (unsigned j = 0; j < _jet_pieces_numerator.size(); j++) {
         _jet_pieces[j] = _jet_pieces_numerator[j]/_denominator;
         _numerator += _jet_pieces_numerator[j];
      }
      
      _beam_piece = _beam_piece_numerator/_denominator;
      _tau = _numerator/_denominator;
   }
   
   
   // return values
   std::vector<double> jet_pieces_numerator() const { return _jet_pieces_numerator; }
   double beam_piece_numerator() const { return _beam_piece_numerator; }
   double denominator() const { return _denominator; }
   double numerator() const { return _numerator; }

   bool has_denominator() const { return _has_denominator; }
   bool has_beam() const { return _has_beam; }
   
   std::vector<double> jet_pieces() const { return _jet_pieces; }
   double beam_piece() const { return _beam_piece; }
   double tau() const { return _tau; }

private:
   
   // these values are input in the constructor
   std::vector<double> _jet_pieces_numerator;
   double _beam_piece_numerator;
   double _denominator;
   bool _has_denominator; //added so that TauComponents knows if denominator is used or not
   bool _has_beam; //added so that TauComponents knows if beam regions is used or not
   
   // these values are derived from above values
   std::vector<double> _jet_pieces;
   double _beam_piece;
   double _numerator;
   double _tau;
   
};

///////
//
// TauPartition
//
///////

// Class for storing partitioning information.
class TauPartition {

public:
   // empty constructor
   TauPartition() {}
   
   TauPartition(int n_jet) {
      _jets_list.resize(n_jet);
      _jets_partition.resize(n_jet);
   }
   
   void push_back_jet(int jet_num, const PseudoJet& part_to_add, int part_index) {
      _jets_list[jet_num].push_back(part_index);
      _jets_partition[jet_num].push_back(part_to_add);
   }
   
   void push_back_beam(const PseudoJet& part_to_add, int part_index) {
      _beam_list.push_back(part_index);
      _beam_partition.push_back(part_to_add);
   }
   
   PseudoJet jet(int jet_num) const { return join(_jets_partition.at(jet_num)); }
   PseudoJet beam() const { return join(_beam_partition);}
   
   std::vector<PseudoJet> jets() const {
      std::vector<PseudoJet> jets;
      for (int i = 0; i < _jets_partition.size(); i++) {
         jets.push_back(jet(i));
      }
      return jets;
   }
   
   const std::list<int> & jet_list(int jet_num) const { return _jets_list.at(jet_num);}
   const std::list<int> & beam_list() const { return _beam_list;}
   const std::vector<std::list<int> > & jets_list() const { return _jets_list;}
   
private:
   
   std::vector<std::list<int> > _jets_list;
   std::list<int> _beam_list;
  
   std::vector<std::vector<PseudoJet> > _jets_partition;
   std::vector<PseudoJet> _beam_partition;
   
};

///////
//
// TauExtras
//
///////

//------------------------------------------------------------------------
/// \class TauExtras
class TauExtras{

public:
   TauExtras(TauComponents & components, TauPartition & partition, std::vector<PseudoJet> & axes, std::vector<PseudoJet> & seed_axes)
   :  _components(components),
      _partition(partition),
      _axes(axes),
      _seed_axes(seed_axes) {
      _jets = partition.jets();  // this is to avoid multiple calculations, since TauPartition calculates it on the fly
   }
   
   const TauComponents & components() const {return _components;}
   const TauPartition & partition() const {return _partition;}
   
   const std::vector<PseudoJet> & jets() const {return _jets;}
   const std::vector<PseudoJet> & axes() const {return _axes;}
   const std::vector<PseudoJet> & seedAxes() const {return _seed_axes;}
   
   
private:
   
   TauComponents _components;
   TauPartition _partition;

   std::vector<PseudoJet> _jets;

   
   std::vector<PseudoJet> _axes;
   std::vector<PseudoJet> _seed_axes;

};
   
   
///////
//
// NjettinessExtras
//
///////
// This class contains the same information as Njettiness, but redoes it in terms of the ClusterSequence::Extras class.
// This is done in order to help improve the interface for the main NjettinessPlugin class.
// TODO:  This class should probably be merged with TauComponents, since both have access
// to similar information
class NjettinessExtras : public ClusterSequence::Extras {
   
public:
   NjettinessExtras(TauComponents tau_components, std::vector<fastjet::PseudoJet> jets, std::vector<fastjet::PseudoJet> axes) : _tau_components(tau_components), _jets(jets), _axes(axes) {}
   
   double totalTau() const {return _tau_components.tau();}
   std::vector<double> subTaus() const {return _tau_components.jet_pieces();}
   std::vector<fastjet::PseudoJet> jets() const {return _jets;}
   std::vector<fastjet::PseudoJet> axes() const {return _axes;}
   
   double totalTau(const fastjet::PseudoJet& /*jet*/) const {
      return _tau_components.tau();
   }
   
   double subTau(const fastjet::PseudoJet& jet) const {
      if (labelOf(jet) == -1) return std::numeric_limits<double>::quiet_NaN(); // nonsense
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


   
   
   
} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_TAUCOMPONENTS_HH__
