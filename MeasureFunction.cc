//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
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


#include "MeasureFunction.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Measure Function
//
///////

// Returns the unnormalized sub-tau values, i.e. a std::vector of the contributions to tau_N of each Voronoi region (or region within R_0)
std::vector<double> MeasureFunction::subtaus_numerator(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes) {
      
   std::vector<double> tauNum(axes.size(), 0.0);
   
   for (unsigned i = 0; i < particles.size(); i++) {
      // find minimum distance; start with 0'th axis for reference
      int j_min = 0;
      double minR = distance(particles[i],axes[0]);
      for (unsigned j = 1; j < axes.size(); j++) {
         double tempR = distance(particles[i],axes[j]); // delta R distance
         if (tempR < minR) {minR = tempR; j_min = j;}
      }
      tauNum[j_min] += numerator(particles[i],axes[j_min]);
   }

   return tauNum;
}

//Calculates normalization for tau and subTau if _has_denominator is true, otherwise returns 1.0 (i.e. no normalization)
double MeasureFunction::tau_denominator(const std::vector <fastjet::PseudoJet>& particles) {
   if (_has_denominator) {
      double tauDen = 0.0;
      for (unsigned i = 0; i < particles.size(); i++) {
         tauDen += denominator(particles[i]);
      }
      return tauDen;
   }
   else return 1.0;
}

// Return all of the necessary TauComponents for specific input particles and axes
TauComponents MeasureFunction::result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) {
   std::vector<double> _subtaus_numerator = subtaus_numerator(particles, axes);
   double _tau_denominator = tau_denominator(particles);
   return TauComponents(_subtaus_numerator, _tau_denominator, _has_denominator);
}
   
} //namespace contrib

FASTJET_END_NAMESPACE
