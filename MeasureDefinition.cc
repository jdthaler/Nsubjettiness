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


//descriptions updated to include measure type -- TJW   
std::string DefaultMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Default Measure (should not be used directly)";
   return stream.str();
};

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

std::string DeprecatedGeometricMeasure::description() const {
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

std::string DeprecatedGeometricCutoffMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Geometric Cutoff Measure (beta = " << _jet_beta << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
};

// Added by TJW

std::string OriginalGeometricMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Original Geometric Cutoff Measure (Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
}; 

std::string ModifiedGeometricMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Modified Geometric Cutoff Measure (Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
}; 

std::string ConicalGeometricCutoffMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Conical Geometric Cutoff Measure (beta = " << _jet_beta << ", gamma = " << _jet_gamma << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
}; 
   

std::string XConeCutoffMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "XCone Cutoff Measure (beta = " << _jet_beta << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
}; 

AxesRefiner* MeasureDefinition::createAxesRefiner(int nPass) const {
   return (new GeneralAxesRefiner(this->create(), nPass));
}

AxesRefiner* DefaultMeasure::createAxesRefiner(int nPass) const {
   return (new ConicalAxesRefiner(_beta, _Rcutoff, nPass));
}

AxesRefiner* DeprecatedGeometricCutoffMeasure::createAxesRefiner(int nPass) const {
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
      if (has_beam()) minRsq = beam_distance_squared(particles[i]);
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
         assert(has_beam());  // should have beam for this to make sense.
         myPartition.push_back_beam(particles[i],i);
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
   if (!has_denominator()) tauDen = 1.0;  // if no denominator, then 1.0 for no normalization factor
   
   // first find jet pieces
   for (unsigned j = 0; j < axes.size(); j++) {
      std::vector<PseudoJet> thisPartition = partition.jet(j).constituents();
      for (unsigned i = 0; i < thisPartition.size(); i++) {
         jetPieces[j] += jet_numerator(thisPartition[i],axes[j]); //numerator jet piece
         if (has_denominator()) tauDen += denominator(thisPartition[i]); // denominator
      }
   }
   
   // then find beam piece
   if (has_beam()) {
      std::vector<PseudoJet> beamPartition = partition.beam().constituents();

      for (unsigned i = 0; i < beamPartition.size(); i++) {
         beamPiece += beam_numerator(beamPartition[i]); //numerator beam piece
         if (has_denominator()) tauDen += denominator(beamPartition[i]); // denominator
      }
   }
   
   // create jets for storage in TauComponents
   std::vector<PseudoJet> jets = partition.jets();
   
   return TauComponents(_tau_mode, jetPieces, beamPiece, tauDen, jets, axes);
}

// new methods added to generalize energy and angle squared for different measure types -- TJW
double DefaultMeasure::energy(const PseudoJet& jet) const {
   if (_measure_type == pt_R) {
      return jet.perp();
   }  else if (_measure_type == E_theta || _measure_type == lorentz_dot) {
      return jet.e();
   } else {
      assert(_measure_type == pt_R || _measure_type == E_theta || _measure_type == lorentz_dot);
      return std::numeric_limits<double>::quiet_NaN();
   }
}
   
double DefaultMeasure::angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const {
   if (_measure_type == pt_R) {
      return jet1.squared_distance(jet2);
   } else if (_measure_type == E_theta) {
      // doesn't seem to be a fastjet built in for this
      double dot = jet1.px()*jet2.px() + jet1.py()*jet2.py() + jet1.pz()*jet2.pz();
      double norm1 = sqrt(jet1.px()*jet1.px() + jet1.py()*jet1.py() + jet1.pz()*jet1.pz());
      double norm2 = sqrt(jet2.px()*jet2.px() + jet2.py()*jet2.py() + jet2.pz()*jet2.pz());
        
      double costheta = dot/(norm1 * norm2);
      if (costheta > 1.0) costheta = 1.0; // Need to handle case of numerical overflow
      double theta = acos(costheta);
      return theta*theta;   
   } else if (_measure_type == lorentz_dot) {
      double dotproduct = dot_product(jet1,jet2);
      return 2.0 * dotproduct / (jet1.e() * jet2.e());
   } else {
      assert(_measure_type == pt_R || _measure_type == E_theta);
      return std::numeric_limits<double>::quiet_NaN();
   }
}

} //namespace contrib

FASTJET_END_NAMESPACE
