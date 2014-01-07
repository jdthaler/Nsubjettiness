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

//NEW FILE CREATED BY TJW 12/25
//Update to move MeasureFunction class and definitions into separate .cc/.hh files

#include "MeasureFunction.hh" //new .hh file added by TJW 12/25

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Measure Function
//
///////

// Name of MeasureFunctor changed to MeasureFunction -- TJW 12/25
// all base MeasureFunction class function definitions moved from Nsubjettiness.cc -- TJW 12/25

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

// Calculates unnormalized tau_N by summing over subtaus
double MeasureFunction::tau_numerator(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) {
   std::vector<double> tau_vec = subtaus_numerator(particles, axes);
   double tauNum = 0.0;
   for (unsigned i = 0; i < tau_vec.size(); i++) {
      tauNum += tau_vec[i];
   }
   return tauNum;
}

//Calculates normalization for tau and subTau  
double MeasureFunction::tau_denominator(const std::vector <fastjet::PseudoJet>& particles) {
   double tauDen = 0.0;
   for (unsigned i = 0; i < particles.size(); i++) {
      tauDen += denominator(particles[i]);
   }
   return tauDen;
}

//normalizes subTaus according to the denominator calculated above
std::vector<double> MeasureFunction::subtaus_normalized(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet> &axes){
   std::vector<double> tauNum = subtaus_numerator(particles, axes);
   std::vector<double> tau_normalized(axes.size(), 0.0);
   double tauDen = tau_denominator(particles);
   for (unsigned i = 0; i < axes.size(); i++) {
      tau_normalized[i] = tauNum[i]/tauDen;
   }
   return tau_normalized;
}

//normalizes tau_N according to denominator calculated above.
double MeasureFunction::tau_normalized(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet> &axes){
   return tau_numerator(particles, axes)/tau_denominator(particles);
}

} //namespace contrib

FASTJET_END_NAMESPACE
