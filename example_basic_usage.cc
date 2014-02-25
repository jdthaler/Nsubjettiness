//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-13
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


#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <sstream>
#include "Nsubjettiness.hh" // In external code, this should be fastjet/contrib/Nsubjettiness.hh
#include "Njettiness.hh"
#include "NjettinessPlugin.hh"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
void PrintJets(const vector <PseudoJet>& jets, bool min_axes);
void analyze(const vector<PseudoJet> & input_particles);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  //----------------------------------------------------------
  // illustrate how Nsubjettiness contrib works

  analyze(event);

  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {return;}
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    event.push_back(particle);
  }
}

// Helper Function for Printing out Jet Information
void PrintJets(const vector <PseudoJet>& jets, bool commentOut = false);

////////
//
//  Main Routine for Analysis 
//
///////

void analyze(const vector<PseudoJet> & input_particles) {
   
   ////////
   //
   //  Start of analysis.  First find anti-kT jets, then find N-subjettiness values of those jets
   //
   ///////
   
   // Initial clustering with anti-kt algorithm
   JetAlgorithm algorithm = antikt_algorithm;
   double jet_rad = 1.00; // jet radius for anti-kt algorithm
   JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
   ClusterSequence clust_seq(input_particles,jetDef);
   vector<PseudoJet> antikt_jets  = sorted_by_pt(clust_seq.inclusive_jets());
   
   for (int j = 0; j < 2; j++) { // Two hardest jets per event
      if (antikt_jets[j].perp() < 200) continue;
      
      vector<PseudoJet> jet_constituents = clust_seq.constituents(antikt_jets[j]);
      
      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << "Analyzing Jet " << j + 1 << ":" << endl;
      cout << "-------------------------------------------------------------------------------------" << endl;
      
      
      ////////
      //
      //  Basic checks of tau values first
      //
      //  If you don't want to know the directions of the subjets,
      //  then you can use the simple function Nsubjettiness.
      //
      //  Recommended usage for Nsubjettiness:
      //  AxesMode:  kt_axes, wta_kt_axes, onepass_kt_axes, or onepass_wta_kt_axes
      //  MeasureMode:  unnormalized_measure
      //  beta with kt_axes: 2.0
      //  beta with wta_kt_axes: anything greater than 0.0 (particularly good for 1.0)
      //  beta with onepass_kt_axes or onepass_wta_kt_axes:  between 1.0 and 3.0
      //
      ///////
      
      
      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << "N-subjettiness with Unnormalized Measure (in GeV)" << endl;
      cout << "beta = 1.0:  One-pass Winner-Take-All kT Axes" << endl;
      cout << "beta = 2.0:  One-pass E-Scheme kT Axes" << endl;
      cout << "-------------------------------------------------------------------------------------" << endl;
      
      
      // Now loop through all options
      cout << setprecision(6) << right << fixed;
      
      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << setw(14) << "beta"
         << setw(14) << "tau1"
         << setw(14) << "tau2"
         << setw(14) << "tau3"
         << setw(14) << "tau2/tau1"
         << setw(14) << "tau3/tau2"
         << endl;
      
      
      Njettiness::AxesMode axisMode1 = Njettiness::onepass_wta_kt_axes;
      Njettiness::AxesMode axisMode2 = Njettiness::onepass_kt_axes;
      Njettiness::MeasureMode measureMode = Njettiness::unnormalized_measure;
      double beta1 = 1.0;
      double beta2 = 2.0;
      
      // define Nsubjettiness functions (beta = 1.0)
      Nsubjettiness         nSub1_beta1(1,  axisMode1,measureMode,beta1);
      Nsubjettiness         nSub2_beta1(2,  axisMode1,measureMode,beta1);
      Nsubjettiness         nSub3_beta1(3,  axisMode1,measureMode,beta1);
      NsubjettinessRatio   nSub21_beta1(2,1,axisMode1,measureMode,beta1);
      NsubjettinessRatio   nSub32_beta1(3,2,axisMode1,measureMode,beta1);
      
      // define Nsubjettiness functions (beta = 2.0)
      Nsubjettiness         nSub1_beta2(1,  axisMode2,measureMode,beta2);
      Nsubjettiness         nSub2_beta2(2,  axisMode2,measureMode,beta2);
      Nsubjettiness         nSub3_beta2(3,  axisMode2,measureMode,beta2);
      NsubjettinessRatio   nSub21_beta2(2,1,axisMode2,measureMode,beta2);
      NsubjettinessRatio   nSub32_beta2(3,2,axisMode2,measureMode,beta2);
      
      
      // calculate Nsubjettiness values (beta = 1.0)
      double tau1_beta1 = nSub1_beta1(antikt_jets[j]);
      double tau2_beta1 = nSub2_beta1(antikt_jets[j]);
      double tau3_beta1 = nSub3_beta1(antikt_jets[j]);
      double tau21_beta1 = nSub21_beta1(antikt_jets[j]);
      double tau32_beta1 = nSub32_beta1(antikt_jets[j]);

      // calculate Nsubjettiness values (beta = 2.0)
      double tau1_beta2 = nSub1_beta2(antikt_jets[j]);
      double tau2_beta2 = nSub2_beta2(antikt_jets[j]);
      double tau3_beta2 = nSub3_beta2(antikt_jets[j]);
      double tau21_beta2 = nSub21_beta2(antikt_jets[j]);
      double tau32_beta2 = nSub32_beta2(antikt_jets[j]);
      
      // Output results (beta = 1.0)
      cout << setw(14) << beta1
         << setw(14) << tau1_beta1
         << setw(14) << tau2_beta1
         << setw(14) << tau3_beta1
         << setw(14) << tau21_beta1
         << setw(14) << tau32_beta1
         << endl;

      // Output results (beta = 2.0)
      cout << setw(14) << beta2
         << setw(14) << tau1_beta2
         << setw(14) << tau2_beta2
         << setw(14) << tau3_beta2
         << setw(14) << tau21_beta2
         << setw(14) << tau32_beta2
         << endl;
      
      cout << "-------------------------------------------------------------------------------------" << endl;

   }

}
