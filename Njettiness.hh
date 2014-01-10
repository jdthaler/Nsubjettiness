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
      wta_onepass_kt_axes, //one-pass minimization of WTA axes with kt
      wta_onepass_ca_axes, //one-pass minimization of WTA axes with ca
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
      geometric_measure //geometric measure
   };

private:
   MeasureFunction* _functor;
   AxesFinder* _axesFinder;

   std::vector<fastjet::PseudoJet> _currentAxes;

   double _current_tau_normalized;
   double _current_tau_numerator; //To return unnormalized values if wanted
   double _current_tau_denominator; //To return normalization factor if wanted

   std::vector<double> _current_subtaus_normalized; 
   std::vector<double> _current_subtaus_numerator; //To return unnormalized values if wanted
   
   void establishAxes(unsigned n_jets, const std::vector <fastjet::PseudoJet> & inputs);
   void establishTaus(const std::vector <fastjet::PseudoJet> & inputs);
   
public:
   Njettiness(MeasureFunction* functor, AxesFinder* axesFinder) : _functor(functor), _axesFinder(axesFinder) {}

   //updated constructor to use three separate parameters instead of NsubParameters -- TJW 1/9
   Njettiness(AxesMode axes, double beta, double R0, double Rcutoff);
   Njettiness(NsubGeometricParameters paraGeo);

   //new constructor to include both AxesMode and MeasureMode enums, and parameters for them -- TJW 1/7
   Njettiness(AxesMode axes_mode, MeasureMode measure_mode, double para1 = NAN, double para2 = NAN, double para3 = NAN);

   ~Njettiness();
   
   void setMeasureFunction(MeasureFunction* newFunctor) {_functor = newFunctor;}
   void setAxesFinder(AxesFinder* newAxesFinder) {_axesFinder = newAxesFinder;}
   
   // setAxes for Manual mode
   void setAxes(std::vector<fastjet::PseudoJet> myAxes) {
      _currentAxes = myAxes;
   }
   
   // The value of N-subjettiness
   double getTau(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) {
      if (inputJets.size() <= n_jets) {
         _currentAxes = inputJets;
         _currentAxes.resize(n_jets,fastjet::PseudoJet(0.0,0.0,0.0,0.0));
         return 0.0;
      }
      establishAxes(n_jets, inputJets);  // sets current Axes
      establishTaus(inputJets); // sets current Tau Values
      
      return _current_tau_normalized;
   }

   // Alternative function call to return just numerator information
   // Function for retrieving the unnormalized tau_N
   double getTauNumerator(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) { 
      getTau(n_jets,inputJets);      
      return _current_tau_numerator;
   }

   // get axes used by getTau.
   std::vector<fastjet::PseudoJet> currentAxes() { return _currentAxes;}
   
   // get subTau values calculated in getTau.
   std::vector<double> currentTaus() { return _current_subtaus_normalized; }

   // get total Tau value calculated in getTau.
   double currentTau() { return _current_tau_normalized; }

   double currentTauNormalized() { return _current_tau_normalized; }
   double currentTauNumerator() { return _current_tau_numerator; }
   double currentTauDenominator() { return _current_tau_denominator; }
   std::vector<double> currentSubTausNumerator() { return _current_subtaus_numerator; }
   std::vector<double> currentSubTausNormalized() { return _current_subtaus_normalized; }


   // partition inputs by Voronoi (each vector stores indices corresponding to inputJets)
   std::vector<std::list<int> > getPartition(const std::vector<fastjet::PseudoJet> & inputJets);

   // partition inputs by Voronoi
   std::vector<fastjet::PseudoJet> getJets(const std::vector<fastjet::PseudoJet> & inputJets);

};

//all Njettiness function definitions moved to Njettiness.cc -- TJW 12/22

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

