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
#include "fastjet/WrappedStructure.hh"


#include <cmath>
#include <vector>
#include <list>
#include <limits>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

// Classes defined in this file.
class TauComponents;
class TauPartition;
class NjettinessExtras;


// Specified whether there are beam or denominators
enum TauMode {
   DEFAULT = -1, // Added so that constructor would default to some value -- TJW 
   UNNORMALIZED_JET_SHAPE = 0,
   NORMALIZED_JET_SHAPE = 1,
   UNNORMALIZED_EVENT_SHAPE = 2,
   NORMALIZED_EVENT_SHAPE = 3,
};
   
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
   TauComponents() {}
   
   // This constructor takes input vector and double and calculates all necessary tau components
   TauComponents(TauMode tau_mode,
                 const std::vector<double> & jet_pieces_numerator,
                 double beam_piece_numerator,
                 double denominator,
                 const std::vector<PseudoJet> & jets,
                 const std::vector<PseudoJet> & axes
                 );
   
   // test for denominator/beams
   bool has_denominator() const;
   bool has_beam() const;
   
   // return values
   double tau() const { return _tau; }
   std::vector<double> jet_pieces() const { return _jet_pieces; }
   double beam_piece() const { return _beam_piece; }
   
   std::vector<double> jet_pieces_numerator() const { return _jet_pieces_numerator; }
   double beam_piece_numerator() const { return _beam_piece_numerator; }
   double numerator() const { return _numerator; }
   double denominator() const { return _denominator; }

   // information on jets and axes
   PseudoJet total_jet() const { return _total_jet;}
   std::vector<PseudoJet> jets() const { return _jets;}
   std::vector<PseudoJet> axes() const { return _axes;}

   class StructureType;
   
protected:
   
   // Defines whether there is a beam or denominator
   TauMode _tau_mode;
   
   // these values are input in the constructor
   std::vector<double> _jet_pieces_numerator;
   double _beam_piece_numerator;
   double _denominator;
   
   // these values are derived from above values
   std::vector<double> _jet_pieces;
   double _beam_piece;
   double _numerator;
   double _tau;
   
   // Information about jets and axes
   PseudoJet _total_jet;
   std::vector<PseudoJet> _jets;
   std::vector<PseudoJet> _axes;
   
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
   
   // Make partition of size to hold n_jet partitions
   TauPartition(int n_jet) {
      _jets_list.resize(n_jet);
      _jets_partition.resize(n_jet);
   }
   
   // add a particle to the jet
   void push_back_jet(int jet_num, const PseudoJet& part_to_add, int part_index) {
      _jets_list[jet_num].push_back(part_index);
      _jets_partition[jet_num].push_back(part_to_add);
   }
   
   // add a particle to the beam
   void push_back_beam(const PseudoJet& part_to_add, int part_index) {
      _beam_list.push_back(part_index);
      _beam_partition.push_back(part_to_add);
   }
   
   // return jet and beam regions
   PseudoJet jet(int jet_num) const { return join(_jets_partition.at(jet_num)); }
   PseudoJet beam() const { return join(_beam_partition);}
   
   std::vector<PseudoJet> jets() const {
      std::vector<PseudoJet> jets;
      for (unsigned int i = 0; i < _jets_partition.size(); i++) { // made unsigned int to get rid of compiler warnings -- TJW
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
// NjettinessExtras
//
///////
// This class contains the same information as TauComponents, but adds additional ways of linking up
// the jets found in the ClusterSequence::Extras class.
// This is done in order to help improve the interface for the main NjettinessPlugin class.
   
class NjettinessExtras : public ClusterSequence::Extras, public TauComponents {
   
public:
   NjettinessExtras(TauComponents tau_components,
                    std::vector<int> cluster_hist_indices)
   : TauComponents(tau_components), _cluster_hist_indices(cluster_hist_indices) {}
   
   
   
   // Ask for tau of the whole event, but by querying a jet
   double tau(const fastjet::PseudoJet& /*jet*/) const {return _tau;}

   // Ask for tau of an individual jet
   double tau_piece(const fastjet::PseudoJet& jet) const {
      if (labelOf(jet) == -1) return std::numeric_limits<double>::quiet_NaN(); // nonsense
      return _jet_pieces[labelOf(jet)];
   }

   fastjet::PseudoJet axis(const fastjet::PseudoJet& jet) const {
      return _axes[labelOf(jet)];
   }
   
   bool has_njettiness_extras(const fastjet::PseudoJet& jet) const {
      return (labelOf(jet) >= 0);
   }
   
private:
   
   // Store cluster history indices to link up with ClusterSequence
   std::vector<int> _cluster_hist_indices;
   
   // Figure out which jet things belonged to
   int labelOf(const fastjet::PseudoJet& jet) const {
      int thisJet = -1;
      for (unsigned int i = 0; i < _jets.size(); i++) {
         if (_cluster_hist_indices[i] == jet.cluster_hist_index()) {
            thisJet = i;
            break;
         }
      }
      return thisJet;
   }

public:
   
   // These are old methods for gaining this information
   // The recommended interface is given in TauComponents
   
   double totalTau() const {return _tau;}
   std::vector<double> subTaus() const {return _jet_pieces;}
   
   double totalTau(const fastjet::PseudoJet& /*jet*/) const {
      return _tau;
   }
   
   double subTau(const fastjet::PseudoJet& jet) const {
      if (labelOf(jet) == -1) return std::numeric_limits<double>::quiet_NaN(); // nonsense
      return _jet_pieces[labelOf(jet)];
   }
   
   double beamTau() const {
      return _beam_piece;
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

///////
//
// TauComponents::StructureType
//
///////
   
// Small wrapped structure to store tau information
// TODO:  Can these be auto-joined?
class TauComponents::StructureType : public WrappedStructure {

public:
   StructureType(const PseudoJet& j) :
      WrappedStructure(j.structure_shared_ptr())
   {}
   
   double tau_piece() const { return _tau_piece; }

   // alternative call, though might be confusing
   double tau() const { return _tau_piece; }

private:
   friend class TauComponents;
   double _tau_piece;
};
   
   
   
   
} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_TAUCOMPONENTS_HH__
