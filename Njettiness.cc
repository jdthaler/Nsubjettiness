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
//Update to nsubjettiness so that all inline functions are declared explicitly in this file

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

//created separate function to set MeasureFunction and AxesFinder from modes in order to reduce redundant code in Njettiness constructors -- TJW 1/11
void Njettiness::setMeasureFunctionandAxesFinder(AxesMode axes_mode, MeasureMode measure_mode, double para1, double para2, double para3) {

   // definition of maximum Rcutoff for non-cutoff measures -- TJW 1/11
   double max_Rcutoff = std::numeric_limits<double>::max();

   //choose which MeasureFunction to use (added separate DefaultNormalizatedMeasure and DefaultUnnormalizedMeasure classes in MeasureFunction -- TJW 1/8)
   // added normalized_cutoff_measure, unnormalized_cutoff_measure, and geometric_cutoff_measure to use an explicit Rcutoff (as opposed to the largest possible number) -- TJW 1/10
   // value of Rcutoff set differently depending on measure chosen to reduce the number of measure checks necessary in the AxesFinder section -- TJW 1/11
   double Rcutoff;

   switch (measure_mode) {
      case normalized_measure:
         Rcutoff = max_Rcutoff; //no explicit Rcutoff, so set to maximum
         if(correctParameterCount(2, para1, para2, para3)) 
            _functor = new DefaultNormalizedMeasure(para1, para2, Rcutoff); //normalized_measure requires 2 parameters, beta and R0
         else { 
            std::cerr << "normalized_measure needs 2 parameters (beta and R0)" << std::endl;
            exit(1); }
         break;
      case unnormalized_measure:
         Rcutoff = max_Rcutoff; //no explicit Rcutoff, so set to maximum
         if(correctParameterCount(1, para1, para2, para3)) 
            _functor = new DefaultUnnormalizedMeasure(para1, Rcutoff); //unnormalized_measure requires 1 parameter, beta
         else {
            std::cerr << "unnormalized_measure needs 1 parameter (beta)" << std::endl;
            exit(1); }
         break;
      case geometric_measure:
         Rcutoff = max_Rcutoff; //no explicit Rcutoff, so set to maximum
         if(correctParameterCount(0, para1, para2, para3)) 
            _functor = new GeometricMeasure(Rcutoff); //geometric_measure requires 0 parameters
         else {
            std::cerr << "geometric_measure needs 0 parameters" << std::endl;
            exit(1); }
         break;
      case normalized_cutoff_measure:
         Rcutoff = para3; //Rcutoff parameter is 3rd parameter in normalized_cutoff_measure 
         if(correctParameterCount(3, para1, para2, para3)) 
            _functor = new DefaultNormalizedMeasure(para1, para2, Rcutoff); //normalized_cutoff_measure requires 3 parameters, beta, R0, and Rcutoff
         else { 
            std::cerr << "normalized_cutoff_measure has 3 parameters (beta, R0, Rcutoff)" << std::endl;
            exit(1); }
         break;
      case unnormalized_cutoff_measure:
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in normalized_cutoff_measure
         if(correctParameterCount(2, para1, para2, para3)) 
            _functor = new DefaultUnnormalizedMeasure(para1, para2); //unnormalized_cutoff_measure requires 2 parameters, beta and Rcutoff
         else {
            std::cerr << "unnormalized_cutoff_measure has 2 parameters (beta, Rcutoff)" << std::endl;
            exit(1); }
         break;
      case geometric_cutoff_measure:
         Rcutoff = para1; //Rcutoff parameter is 1st parameter in normalized_cutoff_measure
         if(correctParameterCount(1, para1, para2, para3)) 
            _functor = new GeometricMeasure(para1); //geometric_cutoff_measure only requires 1 parameter, Rcutoff
         else {
            std::cerr << "geometric_cutoff_measure has 1 parameter (Rcutoff)" << std::endl;
            exit(1); }
         break;
      default:
         assert(false);
         break;
   }   

   //choose which AxesFinder to use
   // cerr outputs are added to make sure minimization axes finders are only used with normalized_measure, unnormalized_measure, normalized_cutoff_measure, and unnormalized_cutoff_measure -- TJW 1/11
   // Rcutoff is set differently based on measure, so Rcutoff here is already set appropriately; para1 is always beta for non-geometric measures -- TJW 1/11
   switch (axes_mode) {
      case wta_kt_axes:
         _axesFinder = new AxesFinderFromWTA_KT(); 
         break;
      case wta_ca_axes:
         _axesFinder = new AxesFinderFromWTA_CA(); 
         break;
      case onepass_wta_kt_axes:
         if (measure_mode != geometric_measure && measure_mode != geometric_cutoff_measure) _axesFinder = new AxesFinderFromOnePassMinimization(new AxesFinderFromWTA_KT(), para1, Rcutoff);
         else {
            std::cerr << "minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure" << std::endl;
            exit(1); }
         break;
      case onepass_wta_ca_axes:
         if (measure_mode != geometric_measure && measure_mode != geometric_cutoff_measure) _axesFinder = new AxesFinderFromOnePassMinimization(new AxesFinderFromWTA_CA(), para1, Rcutoff);
         else {
            std::cerr << "minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure" << std::endl;
            exit(1); }
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
         if (measure_mode != geometric_measure && measure_mode != geometric_cutoff_measure) _axesFinder = new AxesFinderFromOnePassMinimization(new AxesFinderFromKT(), para1, Rcutoff);
         else {
            std::cerr << "minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure" << std::endl;
            exit(1); }
         break;
      case onepass_ca_axes:
         if (measure_mode != geometric_measure && measure_mode != geometric_cutoff_measure) _axesFinder = new AxesFinderFromOnePassMinimization(new AxesFinderFromCA(), para1, Rcutoff);
         else {
            std::cerr << "minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure" << std::endl;
            exit(1); }
         break;
      case onepass_antikt_0p2_axes:
         if (measure_mode != geometric_measure && measure_mode != geometric_cutoff_measure) _axesFinder = new AxesFinderFromOnePassMinimization(new AxesFinderFromAntiKT(0.2), para1, Rcutoff);
         else {
            std::cerr << "minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure" << std::endl;
            exit(1); }
         break;
      case onepass_manual_axes:
         if (measure_mode != geometric_measure && measure_mode != geometric_cutoff_measure) _axesFinder = new AxesFinderFromOnePassMinimization(new AxesFinderFromUserInput(), para1, Rcutoff);
         else {
            std::cerr << "minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure" << std::endl;
            exit(1); }
         break;
      case min_axes:
         if (measure_mode != geometric_measure && measure_mode != geometric_cutoff_measure) _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromKT(),KmeansParameters(100,0.0001,1000,0.8), para1, Rcutoff);
         else {
            std::cerr << "minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure" << std::endl;
            exit(1); }
         break;
      case manual_axes:
         _axesFinder = new AxesFinderFromUserInput();
         break;
      default:
         assert(false);
         break;
      }   
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
