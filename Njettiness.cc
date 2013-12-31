//NEW FILE CREATED BY TJW 12/22
//Update to nsubjettiness so that all inline functions are declared explicitly in this file

//legal info below copied directly from Njettiness.hh

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

#include "Njettiness.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

// all base MeasureFunction function definitions moved to MeasureFunction.cc -- TJW 12/25

// WinnerTakeAllRecombiner function definitions moved into WinnerTakeAllRecombiner.cc -- TJW 12/28

// Minimization function definitions moved to AxesFinder.cc -- TJW 12/28

///////
//
// Main Njettiness Class
//
///////

//all following Njettiness functions moved from Njettiness.hh -- TJW 12/22

void Njettiness::establishTaus(const std::vector <fastjet::PseudoJet> & inputs) {

   // These are the basic function calls
   _current_subtaus_numerator = _functor->subtaus_numerator(inputs, _currentAxes); //Set numerator of subTaus from functor
   _current_tau_denominator = _functor->tau_denominator(inputs); //Set denominator from functor
   
   // These are derived quantities
   // (Save some computational time by not recalculating in the _functor)
   _current_tau_normalized = 0.0;
   _current_tau_numerator = 0.0;
   _current_subtaus_normalized.resize(_current_subtaus_numerator.size(),0.0);
   for (unsigned j = 0; j < _current_subtaus_numerator.size(); j++) {
      _current_subtaus_normalized[j] = _current_subtaus_numerator[j]/_current_tau_denominator;
      _current_tau_numerator += _current_subtaus_numerator[j];
      _current_tau_normalized += _current_subtaus_normalized[j];
   }

}

//Use NsubAxesMode to pick which type of axes to use
void Njettiness::establishAxes(unsigned int n_jets, const std::vector <fastjet::PseudoJet> & inputs) {
   _currentAxes = _axesFinder->getAxes(n_jets,inputs,_currentAxes);   
}


Njettiness::Njettiness(NsubGeometricParameters paraGeo) {
   double Rcutoff = paraGeo.Rcutoff();
   _functor = new GeometricMeasure(Rcutoff);
   _axesFinder = new AxesFinderFromGeometricMinimization(new AxesFinderFromKT(),Rcutoff);
}

//Constructor sets KmeansParameters from NsubAxesMode input
Njettiness::Njettiness(AxesMode axes, NsubParameters paraNsub) {

   _functor = new DefaultMeasure(paraNsub);  //Is there a way to do this without pointers?

   // memory management note, AxesFinderFromKmeansMinimization is responsible for deleting its subpointer.
   // TODO: convert to smart pointers

   switch (axes) {
      case WTA_kt_axes:
         _axesFinder = new AxesFinderFromWTA_KT(); 
         break;
      case WTA_ca_axes:
         _axesFinder = new AxesFinderFromWTA_CA(); 
         break;
      case WTA_onepass_kt_axes:
         _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromWTA_KT(), KmeansParameters(1,0.0001,1000,0.8), paraNsub); 
         break;
      case WTA_onepass_ca_axes:
         _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromWTA_CA(), KmeansParameters(1,0.0001,1000,0.8), paraNsub); 
         break;
      case kt_axes:
         _axesFinder = new AxesFinderFromKT();
         break;
      case ca_axes:
         _axesFinder = new AxesFinderFromCA();
         break;
      case antikt_0p2_axes:
         _axesFinder = new AxesFinderFromAntiKT(0.2);      
         break;
      case onepass_kt_axes:
         _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromKT(),KmeansParameters(1,0.0001,1000,0.8), paraNsub);      
         break;
      case onepass_ca_axes:
         _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromCA(),KmeansParameters(1,0.0001,1000,0.8), paraNsub);
         break;
      case onepass_antikt_0p2_axes:
         _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromAntiKT(0.2),KmeansParameters(1,0.0001,1000,0.8), paraNsub);
         break;
      case onepass_manual_axes:
         _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromUserInput(),KmeansParameters(1,0.0001,1000,0.8), paraNsub);
         break;
      case min_axes:
         _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromKT(),KmeansParameters(100,0.0001,1000,0.8), paraNsub);
         break;
      case manual_axes:
         _axesFinder = new AxesFinderFromUserInput();
         break;
      default:
         assert(false);
         break;
   }

}

Njettiness::~Njettiness() {
   delete _functor;
   delete _axesFinder;
}

// Partition a list of particles according to which N-jettiness axis they are closest to.
// Return a vector of length _currentAxes.size() (which should be N).
// Each vector element is a list of ints corresponding to the indices in
// particles of the particles belonging to that jet.
std::vector<std::list<int> > Njettiness::getPartition(const std::vector<fastjet::PseudoJet> & particles) {
   std::vector<std::list<int> > partitions(_currentAxes.size());

   int j_min = -1;
   for (unsigned i = 0; i < particles.size(); i++) {
      // find minimum distance
      double minR = 10000.0; // large number
      for (unsigned j = 0; j < _currentAxes.size(); j++) {
         double tempR = _functor->distance(particles[i],_currentAxes[j]); // delta R distance
         if (tempR < minR) {
            minR = tempR;
            j_min = j;
         }
      }
      if (_functor->do_cluster(particles[i],_currentAxes[j_min])) partitions[j_min].push_back(i);
   }
   return partitions;
}

// Having found axes, assign each particle in particles to an axis, and return a set of jets.
// Each jet is the sum of particles closest to an axis (Njet = Naxes).
std::vector<fastjet::PseudoJet> Njettiness::getJets(const std::vector<fastjet::PseudoJet> & particles) {
   
   std::vector<fastjet::PseudoJet> jets(_currentAxes.size());

   std::vector<std::list<int> > partition = getPartition(particles);
   for (unsigned j = 0; j < partition.size(); ++j) {
      std::list<int>::const_iterator it, itE;
      for (it = partition[j].begin(), itE = partition[j].end(); it != itE; ++it) {
         jets[j] += particles[*it];
      }
   }
   return jets;
}

} // namespace contrib

FASTJET_END_NAMESPACE
