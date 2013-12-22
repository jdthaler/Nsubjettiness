//  Nsubjettiness Package
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

//UPDATES FROM TJW:

#include "Nsubjettiness.hh"
#include "Njettiness.hh"
#include "NjettinessPlugin.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//constructor definition moved from Nsubjettiness.hh file, but maybe easier to define it in class definition? -- TJW
Nsubjettiness::Nsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff, bool normalized)
  : _njettinessFinder(mode, NsubParameters(beta, R0, Rcutoff)), _N(N), _normalized(normalized) {}

//result moved from Nsubjettiness.hh file --TJW
//result allows user to choose whether they want tau numerator or tau normalized
Double32_t Nsubjettiness::result(const PseudoJet& jet) const {
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   if (_normalized) return _njettinessFinder.getTau(_N, particles); 
   else return _njettinessFinder.getTauNumerator(_N, particles);
} 

//ratio moved from Nsubjettiness.hh file --TJW
//ratio result allows user to choose two tau values to create ratio of tau_M/tau_N
//changed return value from double to Double32_t to match Nsubjettiness class -- TJW
Double32_t NsubjettinessRatio::result(const PseudoJet& jet) const {
   double numerator = _nsub_numerator.result(jet);
   double denominator = _nsub_denominator.result(jet);
   return numerator/denominator;
}

} // namespace contrib

FASTJET_END_NAMESPACE
