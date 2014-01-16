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

#ifndef __FASTJET_CONTRIB_NJETTINESS_HH__
#define __FASTJET_CONTRIB_NJETTINESS_HH__

#include "MeasureFunction.hh" //new file added by TJW 12/25
#include "AxesFinder.hh" //new file added by TJW 12/25

#include "fastjet/PseudoJet.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/JetDefinition.hh"
#include <cmath>
#include <vector>
#include <list>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

//parameter classes moved to MeasureFunction.hh to avoid cross-referencing -- TJW 12/28

//base MeasureFunction class and function definitions moved to MeasureFunction.hh/.cc -- TJW 12/25
//Derived MeasureFunction classes moved to MeasureFunction.hh -- TJW/28

//WinnerTakeAllRecombiner class definition moved into WinnerTakeAllRecombiner.hh -- TJW 12/28

//AxesFinder base class definitions moved to AxesFinder.hh -- TJW 12/25
//AxesFinder from exlusive jets and derived class definitions moved to AxesFinder.hh -- TJW 12/28
//AxesFinder from hardest jets and derived class definitions moved to AxesFinder.hh -- TJW 12/28
//Manual AxesFinder class moved to AxesFinder.hh -- TJW 12/28
//Minimum axes classes moved to AxesFinder.hh -- TJW 12/28

// functions for minimization all moved to AxesFinder.cc -- TJW 12/22

///////
//
// Main Njettiness Class
//
///////

//------------------------------------------------------------------------
/// \class Njettiness
// Njettiness uses AxesFinder and MeasureFunction together in order to find tau_N for the event. The user specifies
// which AxesFinder and which MeasureFunction to use in the calculation, and then Njettiness returns tau_N for the event.
// It also can return information about the axes and jets it used in the calculation, as well as information about 
// how the event was partitioned. -- comment added by TJW
class Njettiness {
public:
   enum AxesMode {
      wta_kt_axes, //Winner Take All axes with kt
      wta_ca_axes, // Winner Take All axes with CA
      onepass_wta_kt_axes, //one-pass minimization of WTA axes with kt
      onepass_wta_ca_axes, //one-pass minimization of WTA axes with ca
      kt_axes,  // exclusive kt axes
      ca_axes,  // exclusive ca axes
      antikt_0p2_axes,  // inclusive hardest axes with antikt-0.2
      min_axes, // axes that minimize N-subjettiness (100 passes by default)
      manual_axes, // set your own axes with setAxes()
      onepass_kt_axes, // one-pass minimization from kt starting point
      onepass_ca_axes, // one-pass minimization from ca starting point
      onepass_antikt_0p2_axes,  // one-pass minimization from antikt-0.2 starting point 
      onepass_manual_axes  // one-pass minimization from manual starting point
   };

   //new MeasureMode enum added by TJW 1/7
   enum MeasureMode {
      normalized_measure, //default normalized measure
      unnormalized_measure, //default unnormalized measure
      geometric_measure, //geometric measure
      normalized_cutoff_measure, //default normalized measure with explicit Rcutoff
      unnormalized_cutoff_measure, //default unnormalized measure with explicit Rcutoff
      geometric_cutoff_measure //geometric measure with explicit Rcutoff
   };

private:
   MeasureFunction* _function;
   AxesFinder* _axesFinder;
   // new class added to contain all various components of tau -- TJW 1/14
   TauComponents _current_tau_components; //automatically set to have components of 0; these values will be set by the getTau function call -- TJW 1/15

   // added enum information so functions can specify output based on specific options, primarily for setAxes -- TJW 1/15
   AxesMode _current_axes_mode;
   MeasureMode _current_measure_mode;

   std::vector<fastjet::PseudoJet> _currentAxes;

   // these values are no longer necessary -- TJW 1/15
   // double _current_tau_normalized;
   // double _current_tau_numerator; //To return unnormalized values if wanted
   // double _current_tau_denominator; //To return normalization factor if wanted
   // std::vector<double> _current_subtaus_normalized; 
   // std::vector<double> _current_subtaus_numerator; //To return unnormalized values if wanted
   
   //Use NsubAxesMode to pick which type of axes to use
   // function definition moved from Njettiness.cc -- TJW 1/15
   void establishAxes(unsigned n_jets, const std::vector <fastjet::PseudoJet> & inputs) {
      _currentAxes = _axesFinder->getAxes(n_jets,inputs,_currentAxes);   
   }

   // function removed since it is no longer necessary -- TJW 1/15
   // void establishTaus(const std::vector <fastjet::PseudoJet> & inputs);

   // added for compilation of non C++11 users -- TJW 1/15
   bool isnan(double para) { return para != para; }

   //created new function to check to make sure input has correct number of parameters -- TJW 1/10
   bool correctParameterCount(int n, double para1, double para2, double para3, double para4){
      int numpara;
      if (!isnan(para1) && !isnan(para2) && !isnan(para3) && !isnan(para4)) numpara = 4;
      else if (!isnan(para1) && !isnan(para2) && !isnan(para3) && isnan(para4)) numpara = 3;
      else if (!isnan(para1) && !isnan(para2) && isnan(para3) && isnan(para4)) numpara = 2;
      else if (!isnan(para1) && isnan(para2) && isnan(para3) && isnan(para4)) numpara = 1;
      else numpara = 0;
      return n == numpara;
   }

   // created new function to set onepass_axes depending on input measure_mode and startingFinder-- TJW 1/13
   // made void so that it just sets _axesFinder instead of returning AxesFinder -- TJW 1/15
   void setOnePassAxesFinder(MeasureMode measure_mode, AxesFinder* startingFinder, double para1, double Rcutoff) {
      if (measure_mode == normalized_measure || measure_mode == unnormalized_measure || measure_mode == normalized_cutoff_measure || measure_mode == unnormalized_cutoff_measure) {
         _axesFinder = new AxesFinderFromOnePassMinimization(startingFinder, para1, Rcutoff);
      }
      else if (measure_mode == geometric_measure || measure_mode == geometric_cutoff_measure) {
         _axesFinder = new AxesFinderFromGeometricMinimization(startingFinder, Rcutoff);
      }
      else {
         std::cerr << "minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure, geometric_measure, geometric_cutoff_measure" << std::endl;
         exit(1); }
   }
 
   //created separate function to set MeasureFunction and AxesFinder in order to reduce redundant code in Njettiness constructors -- TJW 1/11
   void setMeasureFunctionandAxesFinder(AxesMode axes_mode, MeasureMode measure_mode, double para1, double para2, double para3, double para4);

public:
   Njettiness(AxesFinder* axesFinder, MeasureFunction* function) : _function(function), _axesFinder(axesFinder) {}

   // updated constructor to use three separate parameters instead of NsubParameters -- TJW 1/9
   // updated to use new private function defined above -- TJW 1/11
   // constructor removed since it is never explicitly used by the user -- TJW 1/13
   // Njettiness(AxesMode axes_mode, double beta, double R0, double Rcutoff) {
   //    setMeasureFunctionandAxesFinder(axes_mode, normalized_cutoff_measure, beta, R0, Rcutoff);
   // }

   //new constructor to include both AxesMode and MeasureMode enums, and parameters for them -- TJW 1/7
   // updated to use new private function defined above -- TJW 1/11
   // updated to include 4th parameter (if necessary) -- TJW 
   Njettiness(AxesMode axes_mode, MeasureMode measure_mode, double para1 = NAN, double para2 = NAN, double para3 = NAN, double para4 = NAN) : _current_axes_mode(axes_mode), _current_measure_mode(measure_mode) {
      setMeasureFunctionandAxesFinder(axes_mode, measure_mode, para1, para2, para3, para4);      
   }
   
   // updated constructor to use separate Rcutoff parameter instead of NsubGeometricParameters for initialization of geometric measure-- TJW 1/10
   // constructor definition moved from Njettiness.cc -- TJW 1/11   
   // this constructor should no longer exist because the user should specify geometric_measure with above constructor. -- TJW 1/14
   // Njettiness(double Rcutoff) {
   //    _function = new GeometricMeasure(Rcutoff);
   //    _axesFinder = new AxesFinderFromGeometricMinimization(new AxesFinderFromKT(),Rcutoff);
   // }

   // constructor definition moved from Njettiness.cc -- TJW 1/11
   ~Njettiness() {
      delete _function;
      delete _axesFinder;
   }

   void setMeasureFunction(MeasureFunction* newFunction) {_function = newFunction;}
   void setAxesFinder(AxesFinder* newAxesFinder) {_axesFinder = newAxesFinder;}
   
   // setAxes for Manual mode
   void setAxes(std::vector<fastjet::PseudoJet> myAxes) {
      if (_current_axes_mode == manual_axes || _current_axes_mode == onepass_manual_axes) {
         _currentAxes = myAxes;
      }
      else {
         std::cerr << "You can only use setAxes if using manual_axes or onepass_manual_axes measure mode" << std::endl;
         exit(1);
      }
   }
   
   // The value of N-subjettiness
   double getTau(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) {
      if (inputJets.size() <= n_jets) {
         _currentAxes = inputJets;
         _currentAxes.resize(n_jets,fastjet::PseudoJet(0.0,0.0,0.0,0.0));
         return 0.0;
      }
      establishAxes(n_jets, inputJets);  // sets current Axes
      _current_tau_components = _function->result(inputJets, _currentAxes);  // sets current Tau Values
      //establishTaus(inputJets); //no longer necessary -- TJW 1/15
      
      // return _current_tau_normalized;
      return _current_tau_components.tau_normalized();
   }

   // new function to return all TauComponents that user would want -- TJW 1/15
   TauComponents getTauComponents(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) {
      getTau(n_jets, inputJets);
      return _current_tau_components;
   }

   // Alternative function call to return just numerator information
   // Function for retrieving the unnormalized tau_N
   // removed since nothing uses it -- TJW 1/15
   // double getTauNumerator(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) { 
   //    getTau(n_jets,inputJets);      
   //    return _current_tau_numerator;
   // }

   // get axes used by getTau.
   std::vector<fastjet::PseudoJet> currentAxes() { return _currentAxes;}
   
   // get subTau values calculated in getTau.
   // std::vector<double> currentTaus() { return _current_subtaus_normalize; }
   std::vector<double> currentTaus() { return _current_tau_components.subtaus_normalized(); }

   // get total Tau value calculated in getTau.
   // double currentTau() { return _current_tau_normalized; }
   double currentTau() { return _current_tau_components.tau_normalized(); }

//   TauComponents currentTauComponents() {return _current_tau_components; }

   // these functions aren't used by anyone and can be replaced by NjettinessComponents function -- TJW 1/15
   // double currentTauNormalized() { return _current_tau_normalized; }
   // double currentTauNumerator() { return _current_tau_numerator; }
   // double currentTauDenominator() { return _current_tau_denominator; }
   // std::vector<double> currentSubTausNumerator() { return _current_subtaus_numerator; }
   // std::vector<double> currentSubTausNormalized() { return _current_subtaus_normalized; }


   // partition inputs by Voronoi (each vector stores indices corresponding to inputJets)
   std::vector<std::list<int> > getPartition(const std::vector<fastjet::PseudoJet> & inputJets);

   // partition inputs by Voronoi
   std::vector<fastjet::PseudoJet> getJets(const std::vector<fastjet::PseudoJet> & inputJets);

};

//all Njettiness function definitions moved to Njettiness.cc -- TJW 12/22

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

