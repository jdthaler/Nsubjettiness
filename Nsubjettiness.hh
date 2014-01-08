// Nsubjettiness Package
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

#ifndef __FASTJET_CONTRIB_NSUBJETTINESS_HH__
#define __FASTJET_CONTRIB_NSUBJETTINESS_HH__

#include <fastjet/internal/base.hh>

#include "Njettiness.hh"

#include "fastjet/FunctionOfPseudoJet.hh"
#include <string>
#include <climits>

#ifndef G__DICTIONARY
typedef double Double32_t; // ROOT will store as 32-bit, but in code is double
#endif


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

//------------------------------------------------------------------------
/// \class Nsubjettiness
/// Nsubjettiness extends the concept of Njettiness to a jet shape, but other
/// than the set of particles considered, they are identical.  This class
/// wraps the core Njettiness code to provide the fastjet::FunctionOfPseudoJet
/// interface for convenience in larger analyses.  See NjettinessPlugin.hh for
/// definitions of tau_N and the constructor options.

class Nsubjettiness : public FunctionOfPseudoJet<Double32_t> {
public:

   //moved constructor definition to here to clean up code -- TJW 12/25
   Nsubjettiness(int N, 
      Njettiness::AxesMode mode, 
      double beta, 
      double R0, 
      double Rcutoff=std::numeric_limits<double>::max(), 
      bool normalized = true)
   : _njettinessFinder(mode, NsubParameters(beta, R0, Rcutoff)), _N(N), _normalized(normalized) {}

   //new constructor definition added by TJW 1/7 (normalization done by checking if the measure_mode is set to unnormalized_measure or not, seems a bit messy)
   Nsubjettiness(int N, 
      Njettiness::AxesMode axes_mode, 
      Njettiness::MeasureMode measure_mode, 
      double para1 = NAN, double para2 = NAN, double para3 = NAN) 
   : _njettinessFinder(axes_mode, measure_mode, para1, para2, para3), _N(N), _normalized(measure_mode != Njettiness::unnormalized_measure) {}

   /// returns tau_N, measured on the constituents of this jet 
   Double32_t result(const PseudoJet& jet) const;

   //To set axes for manual use 
   void setAxes(std::vector<fastjet::PseudoJet> myAxes) {
      // TODO:  Have this test that manual axes are being used
   	_njettinessFinder.setAxes(myAxes);
   }

private:

   mutable Njettiness _njettinessFinder; // should muck with this so result can be const without this mutable

   int _N;
   bool _normalized;

};

//result definition moved to Nsubjettiness.cc -- TJW 12/22

//------------------------------------------------------------------------
/// \class NsubjettinessRatios
// NsubjettinessRatios Calculates uses the results from Nsubjettiness to calculate the ratio
// tau_N/tau_M, where N and M are specified by the user. The ratio of different tau values
// proves to be an interesting value for analysis, and this class allows these values
// to be calculated with ease. -- comment added by TJW

class NsubjettinessRatio : public FunctionOfPseudoJet<Double32_t> {
public:

   NsubjettinessRatio(int N, int M, 
      Njettiness::AxesMode mode, 
      double beta, 
      double R0 = 1.0, //added R0 to constructor arguments and default it to 1.0 since it is unimportant for ratio -- TJW 12/28
      double Rcutoff=std::numeric_limits<double>::max())
   : _nsub_numerator(N, mode, beta, R0, Rcutoff), _nsub_denominator(M, mode, beta, R0, Rcutoff) {};
   // TODO:  Set constant for R0 = 1.0 so that we don't have magic numbers

   //changed return value from double to Double32_t to match Nsubjettiness class -- TJW 12/22
   //returns tau_N/tau_M based off the input jet using result function from Nsubjettiness 
   Double32_t result(const PseudoJet& jet) const;

private: 
   Nsubjettiness _nsub_numerator;
   Nsubjettiness _nsub_denominator;
};

//ratio result definition moved to Nsubjettiness.cc --TJW 12/22

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NSUBJETTINESS_HH__
