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


#include "AxesRefiner.hh"
#include "MeasureDefinition.hh"

#include <iomanip>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

///////
//
// Measure Function
//
///////

   
std::string NormalizedMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Normalized Measure (beta = " << _beta << ", R0 = " << _R0 << ")";
   return stream.str();
};

std::string UnnormalizedMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Unnormalized Measure (beta = " << _beta << ", in GeV)";
   return stream.str();
};

std::string GeometricMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Geometric Measure (beta = " << _jet_beta << ", in GeV)";
   return stream.str();
};

std::string NormalizedCutoffMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Normalized Cutoff Measure (beta = " << _beta << ", R0 = " << _R0 << ", Rcut = " << _Rcutoff << ")";
   return stream.str();
};

std::string UnnormalizedCutoffMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Unnormalized Cutoff Measure (beta = " << _beta << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
};

std::string GeometricCutoffMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Geometric Cutoff Measure (beta = " << _jet_beta << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
};
   
   
   
AxesRefiner* NormalizedCutoffMeasure::createAxesRefiner(int nPass) const {
   return (new ConicalAxesRefiner(_beta, _Rcutoff, nPass));
}

   
AxesRefiner* GeometricCutoffMeasure::createAxesRefiner(int nPass) const {
   return (new GeometricAxesRefiner(_jet_beta, _Rcutoff, nPass));
}
   
   
// Return all of the necessary TauComponents for specific input particles and axes
TauComponents MeasureDefinition::component_result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const {
   
   // first find partition
   TauPartition partition = get_partition(particles,axes);
   
   // then return result calculated from partition
   return component_result_from_partition(partition,axes);
}

TauPartition MeasureDefinition::get_partition(const std::vector<fastjet::PseudoJet>& particles,
                                            const std::vector<fastjet::PseudoJet>& axes) const {
   
   TauPartition myPartition(axes.size());
   
   // Figures out the partiting of the input particles into the various jet pieces
   // Based on which axis the parition is closest to
   for (unsigned i = 0; i < particles.size(); i++) {
      
      // find minimum distance; start with beam (-1) for reference
      int j_min = -1;
      double minRsq;
      if (_has_beam) minRsq = beam_distance_squared(particles[i]);
      else minRsq = std::numeric_limits<double>::max(); // make it large value
      
      // check to see which axis the particle is closest to
      for (unsigned j = 0; j < axes.size(); j++) {
         double tempRsq = jet_distance_squared(particles[i],axes[j]); // delta R distance
         if (tempRsq < minRsq) {
            minRsq = tempRsq;
            j_min = j;
         }
      }
      
      if (j_min == -1) {
         if (_has_beam) myPartition.push_back_beam(particles[i],i);
         else assert(_has_beam);  // this should never happen.
      } else {
         myPartition.push_back_jet(j_min,particles[i],i);
      }
   }
   
   return myPartition;
}

// Uses existing partition and calculates result
// TODO:  Can we cache this for speed up when doing area subtraction?
TauComponents MeasureDefinition::component_result_from_partition(const TauPartition& partition,
                                                     const std::vector<fastjet::PseudoJet>& axes) const {
   
   std::vector<double> jetPieces(axes.size(), 0.0);
   double beamPiece = 0.0;
   
   double tauDen = 0.0;
   if (!_has_denominator) tauDen = 1.0;  // if no denominator, then 1.0 for no normalization factor
   
   // first find jet pieces
   for (unsigned j = 0; j < axes.size(); j++) {
      std::vector<PseudoJet> thisPartition = partition.jet(j).constituents();
      for (unsigned i = 0; i < thisPartition.size(); i++) {
         jetPieces[j] += jet_numerator(thisPartition[i],axes[j]); //numerator jet piece
         if (_has_denominator) tauDen += denominator(thisPartition[i]); // denominator
      }
   }
   
   // then find beam piece
   if (_has_beam) {
      std::vector<PseudoJet> beamPartition = partition.beam().constituents();

      for (unsigned i = 0; i < beamPartition.size(); i++) {
         beamPiece += beam_numerator(beamPartition[i]); //numerator beam piece
         if (_has_denominator) tauDen += denominator(beamPartition[i]); // denominator
      }
   }
   return TauComponents(jetPieces, beamPiece, tauDen, _has_denominator, _has_beam);
}

   
   
   
   
} //namespace contrib

FASTJET_END_NAMESPACE
