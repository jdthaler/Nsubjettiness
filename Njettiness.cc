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

#include "Njettiness.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {


///////
//
// Main Njettiness Class
//
///////

Njettiness::Njettiness(AxesMode axes_mode, const MeasureDefinition & measure_spec)
: _axes_mode(axes_mode), _measure_spec(measure_spec.copy()) {
   setMeasureFunctionAndAxesFinder();  // call helper function to do the hard work
}
   
Njettiness::~Njettiness() {
   // clean house
   delete _measure_spec;
   delete _measureFunction;
   delete _axesFinder;
}

MeasureDefinition* Njettiness::createMeasureDef(MeasureMode measure_mode, int num_para, double para1, double para2, double para3) const {

   // definition of maximum Rcutoff for non-cutoff measures, changed later by other measures
   double Rcutoff = std::numeric_limits<double>::max();  //large number
   // Most (but not all) measures have some kind of beta value
   double beta = std::numeric_limits<double>::quiet_NaN();
   // The normalized measures have an R0 value.
   double R0 = std::numeric_limits<double>::quiet_NaN();
   
   // Find the MeasureFunction and set the parameters.
   switch (measure_mode) {
      case normalized_measure:
         beta = para1;
         R0 = para2;
         if(num_para == 2) {
            return new NormalizedMeasure(beta,R0);
         } else {
            throw Error("normalized_measure needs 2 parameters (beta and R0)");
         }
         break;
      case unnormalized_measure:
         beta = para1;
         if(num_para == 1) {
            return new UnnormalizedMeasure(beta);
         } else {
            throw Error("unnormalized_measure needs 1 parameter (beta)");
         }
         break;
      case geometric_measure:
         beta = para1;
         if (num_para == 1) {
            return new GeometricMeasure(beta);
         } else {
            throw Error("geometric_measure needs 1 parameter (beta)");
         }
         break;
      case normalized_cutoff_measure:
         beta = para1;
         R0 = para2;
         Rcutoff = para3; //Rcutoff parameter is 3rd parameter in normalized_cutoff_measure
         if (num_para == 3) {
            return new NormalizedCutoffMeasure(beta,R0,Rcutoff);
         } else {
            throw Error("normalized_cutoff_measure has 3 parameters (beta, R0, Rcutoff)");
            exit(1); }
         break;
      case unnormalized_cutoff_measure:
         beta = para1;
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in normalized_cutoff_measure
         if (num_para == 2) {
            return new UnnormalizedCutoffMeasure(beta,Rcutoff);
         } else {
            throw Error("unnormalized_cutoff_measure has 2 parameters (beta, Rcutoff)");
            exit(1); }
         break;
      case geometric_cutoff_measure:
         beta = para1;
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in geometric_cutoff_measure
         if(num_para == 2) {
           return new GeometricCutoffMeasure(beta,Rcutoff);
         } else {
            throw Error("geometric_cutoff_measure has 2 parameters (beta, Rcutoff)");
            exit(1); }
         break;
      default:
         assert(false);
         break;
   }
   return NULL;
}
   
// Parsing needed for constructor to set AxesFinder and MeasureFunction
// All of the parameter handling is here, and checking that number of parameters is correct.
void Njettiness::setMeasureFunctionAndAxesFinder() {

   // Get the correct MeasureFunction
   _measureFunction = _measure_spec->associatedMeasureFunction();
   
   // Choose which AxesFinder from user input.
   // Uses setOnePassAxesFinder helpful function to use beta and Rcutoff values about (if needed)
   switch (_axes_mode) {
      case wta_kt_axes:
         _axesFinder = new AxesFinderFromWTA_KT(); 
         break;
      case wta_ca_axes:
         _axesFinder = new AxesFinderFromWTA_CA(); 
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
      case onepass_wta_kt_axes:
         _axesFinder = _measure_spec->associatedOnePassAxesFinder(new AxesFinderFromWTA_KT());
         break;
      case onepass_wta_ca_axes:
         _axesFinder = _measure_spec->associatedOnePassAxesFinder(new AxesFinderFromWTA_CA());
         break;
      case onepass_kt_axes:
         _axesFinder = _measure_spec->associatedOnePassAxesFinder(new AxesFinderFromKT());
         break;
      case onepass_ca_axes:
         _axesFinder = _measure_spec->associatedOnePassAxesFinder(new AxesFinderFromCA());
         break;
      case onepass_antikt_0p2_axes:
         _axesFinder = _measure_spec->associatedOnePassAxesFinder(new AxesFinderFromAntiKT(0.2));
         break;
      case onepass_manual_axes:
         _axesFinder = _measure_spec->associatedOnePassAxesFinder(new AxesFinderFromUserInput());
         break;
      case min_axes: //full minimization is not defined for geometric_measure.
         //Defaults to 100 iteration to find minimum
         _axesFinder = _measure_spec->associatedMultiPassAxesFinder(new AxesFinderFromKT());
         break;
      case manual_axes:
         _axesFinder = new AxesFinderFromUserInput();
         break;
      default:
         assert(false);
         break;
      }   

}

// setAxes for Manual mode
void Njettiness::setAxes(const std::vector<fastjet::PseudoJet> & myAxes) {
   if (_axes_mode == manual_axes || _axes_mode == onepass_manual_axes) {
      _currentAxes = myAxes;
   }
   else {
      std::cerr << "You can only use setAxes if using manual_axes or onepass_manual_axes measure mode" << std::endl;
      exit(1);
   }
}
   
// Calculates and returns all TauComponents that user would want.
// This information is stored in _current_tau_components for later access as well.
TauComponents Njettiness::getTauComponents(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) {
   if (inputJets.size() <= n_jets) {  //if not enough particles, return zero
      _currentAxes = inputJets;
      _currentAxes.resize(n_jets,fastjet::PseudoJet(0.0,0.0,0.0,0.0));
      _current_tau_components = TauComponents();
      _seedAxes = _currentAxes;
      _currentJets = _currentAxes;
      _currentBeam = PseudoJet(0.0,0.0,0.0,0.0);
   } else {
      
      // Find axes and store information
      _currentAxes = _axesFinder->getAxes(n_jets,inputJets,_currentAxes); // sets current Axes
      _seedAxes = _axesFinder->seedAxes(); // sets seed Axes (if one pass minimization was used)

      
      // Find partition and store information
      // (jet information in _currentJets, beam in _currentBeam)
      _currentJets = _measureFunction->get_partition(inputJets,_currentAxes,&_currentBeam);
      
      // Find tau value and store information
      _current_tau_components = _measureFunction->result_from_partition(_currentJets, _currentAxes,&_currentBeam);  // sets current Tau Values
   }
   return _current_tau_components;
}
   
   
// Partition a list of particles according to which N-jettiness axis they are closest to.
// Return a vector of length _currentAxes.size() (which should be N).
// Each vector element is a list of ints corresponding to the indices in
// particles of the particles belonging to that jet.
std::vector<std::list<int> > Njettiness::getPartitionList(const std::vector<fastjet::PseudoJet> & particles) {
   // core code is in MeasureFunction
   return _measureFunction->get_partition_list(particles,_currentAxes);
}

   
} // namespace contrib

FASTJET_END_NAMESPACE
