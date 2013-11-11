// Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-13
//  Jesse Thaler, Ken Van Tilburg, and Christopher K. Vermilion
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

   Nsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff=std::numeric_limits<double>::max(), bool normalized = true);
   
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

inline Nsubjettiness::Nsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff, bool normalized)
  : _njettinessFinder(mode, NsubParameters(beta, R0, Rcutoff)), _N(N), _normalized(normalized) {}


//result will allow user to choose whether they want numerator or normalized version
inline Double32_t Nsubjettiness::result(const PseudoJet& jet) const
{
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   if (_normalized) return _njettinessFinder.getTau(_N, particles); 
   else return _njettinessFinder.getTauNumerator(_N, particles);
} 

//Class NsubjettinessRatios
//Used to Calculate tau_N/tau_M based off results from class Nsubjettiness
//Requires two integers in constructor (N, M)
class NsubjettinessRatio : public FunctionOfPseudoJet<Double32_t> {
public:
   NsubjettinessRatio(int N, 
      int M, 
      Njettiness::AxesMode mode, 
      double beta, 
      double Rcutoff=std::numeric_limits<double>::max())
   : _nsub_numerator(N, mode, beta, 1.0, Rcutoff), _nsub_denominator(M, mode, beta, 1.0, Rcutoff) {};
   // TODO:  Set constant for R0 = 1.0 so that we don't have magic numbers

   //returns tau_N/tau_M based off the input jet, uses result function from Nsubjettiness 
   double result(const PseudoJet& jet) const;

private: 
   Nsubjettiness _nsub_numerator;
   Nsubjettiness _nsub_denominator;
};

inline double NsubjettinessRatio::result(const PseudoJet& jet) const {
   double numerator = _nsub_numerator.result(jet);
   double denominator = _nsub_denominator.result(jet);
   return numerator/denominator;
}

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NSUBJETTINESS_HH__
