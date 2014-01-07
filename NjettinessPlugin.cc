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

//NEW FILE CREATED BY TJW 12/22
//Update to NsubjettinessPlugin so that all inline functions are declared explicitly in this file

#include "NjettinessPlugin.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//NjettinessPlugin constructors (moved from NjettinessPlugin.hh -- TJW 12/22)
NjettinessPlugin::NjettinessPlugin(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff)
  : _N(N), _njettinessFinder(mode, NsubParameters(beta, R0, Rcutoff))
{}

NjettinessPlugin::NjettinessPlugin(int N, NsubGeometricParameters paraGeo)
  : _N(N), _njettinessFinder(paraGeo)
{}

//NjettinessPlugin functions (moved from NjettinessPlugin.hh -- TJW 12/22)
std::string NjettinessPlugin::description() const {return "NJettiness";}

//clusters the particles according to the Njettiness jet algorithm
void NjettinessPlugin::run_clustering(ClusterSequence& cs) const
{
   std::vector<fastjet::PseudoJet> particles = cs.jets();
   _njettinessFinder.getTau(_N, particles);
   std::vector<std::list<int> > partition = _njettinessFinder.getPartition(particles);

   std::vector<fastjet::PseudoJet> jet_indices_for_extras;

   // output clusterings for each jet
   for (size_t i = 0; i < partition.size(); ++i) {
      std::list<int>& indices = partition[i];
      if (indices.size() == 0) continue;
      std::list<int>::const_iterator it = indices.begin();
      while (indices.size() > 1) {
         int merge_i = indices.back(); indices.pop_back();
         int merge_j = indices.back(); indices.pop_back();
         int newIndex;
         double fakeDij = -1.0;
      
         cs.plugin_record_ij_recombination(merge_i, merge_j, fakeDij, newIndex);

         indices.push_back(newIndex);
      }
      double fakeDib = -1.0;
      
      int finalJet = indices.back();
      cs.plugin_record_iB_recombination(finalJet, fakeDib);
      jet_indices_for_extras.push_back(cs.jets()[finalJet]);  // Get the four vector for the final jets to compare later.
   }

   NjettinessExtras * extras = new NjettinessExtras(_njettinessFinder.currentTau(),jet_indices_for_extras,_njettinessFinder.currentTaus(),_njettinessFinder.currentAxes());
   cs.plugin_associate_extras(std::auto_ptr<ClusterSequence::Extras>(extras));
   
}


} // namespace contrib

FASTJET_END_NAMESPACE
