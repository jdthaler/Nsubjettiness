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
   //  This code will check multiple axes/measure modes
   //  First thing we do is establish the various modes we will check
   //
   ///////
   
   // Getting all of the axes modes to test
   vector<Njettiness::AxesMode> _testAxesModes;
   vector<string> _testAxesNames;
   _testAxesModes.push_back(Njettiness::kt_axes);                  _testAxesNames.push_back("            KT:");
   _testAxesModes.push_back(Njettiness::ca_axes);                  _testAxesNames.push_back("            CA:");
   _testAxesModes.push_back(Njettiness::antikt_0p2_axes);          _testAxesNames.push_back("        AKT0.2:");
   _testAxesModes.push_back(Njettiness::wta_kt_axes);              _testAxesNames.push_back("        WTA KT:");
   _testAxesModes.push_back(Njettiness::wta_ca_axes);              _testAxesNames.push_back("        WTA CA:");
   _testAxesModes.push_back(Njettiness::onepass_kt_axes);          _testAxesNames.push_back("    OnePass KT:");
   _testAxesModes.push_back(Njettiness::onepass_ca_axes);          _testAxesNames.push_back("    OnePass CA:");
   _testAxesModes.push_back(Njettiness::onepass_antikt_0p2_axes);  _testAxesNames.push_back("OnePass AKT0.2:");
   _testAxesModes.push_back(Njettiness::onepass_wta_kt_axes);      _testAxesNames.push_back("OnePass WTA KT:");
   _testAxesModes.push_back(Njettiness::onepass_wta_ca_axes);      _testAxesNames.push_back("OnePass WTA CA:");
   _testAxesModes.push_back(Njettiness::min_axes);                 _testAxesNames.push_back("#     Min Axes:");  // Putting in # because min_axes uses random number seed

   //
   // Note:  Njettiness::min_axes is not guarenteed to give a global
   // minimum, only a local minimum, and different choices of the random
   // number seed can give different results.  For that reason,
   // the one-pass minimization are recommended over min_axes.
   //
   
   // Getting a smaller list of recommended modes
   // These are the ones that are more likely to give sensible results (and are all IRC safe)
   vector<Njettiness::AxesMode> _testKeyAxesModes;
   vector<string> _testKeyAxesNames;
   vector<string> _testKeyAxesNamesLong;
   _testKeyAxesModes.push_back(Njettiness::kt_axes);                    _testKeyAxesNamesLong.push_back("KT Axes:");
   _testKeyAxesModes.push_back(Njettiness::wta_kt_axes);                _testKeyAxesNamesLong.push_back("Winner-Take-All KT Axes:");
   _testKeyAxesModes.push_back(Njettiness::onepass_kt_axes);            _testKeyAxesNamesLong.push_back("One-Pass Minimization starting from KT:");
   _testKeyAxesModes.push_back(Njettiness::onepass_wta_kt_axes);        _testKeyAxesNamesLong.push_back("One-Pass Minimization starting from WTA KT:");


   // Getting some of the measure modes to test
   // (When applied to a single jet we won't test the cutoff measures,
   // since cutoffs aren't typically helpful when applied to single jets)
   vector<Njettiness::MeasureMode> _testMeasureModes;
   vector<double> _testBeta;
   vector<double> _testR0;
   vector<string> _testMeasureNames;
   _testMeasureModes.push_back(Njettiness::normalized_measure);   _testBeta.push_back(1.0); _testR0.push_back(1.0); _testMeasureNames.push_back("Normalized Measure (beta = 1.0, R0 = 1.0):");
   _testMeasureModes.push_back(Njettiness::unnormalized_measure); _testBeta.push_back(1.0); _testR0.push_back(NAN); _testMeasureNames.push_back("Unnormalized Measure (beta = 1.0, in GeV):");
   _testMeasureModes.push_back(Njettiness::normalized_measure);   _testBeta.push_back(2.0); _testR0.push_back(1.0); _testMeasureNames.push_back("Normalized Measure (beta = 2.0, R0 = 1.0):");
   _testMeasureModes.push_back(Njettiness::unnormalized_measure); _testBeta.push_back(2.0); _testR0.push_back(NAN); _testMeasureNames.push_back("Unnormalized Measure (beta = 2.0, in GeV):");
   _testMeasureModes.push_back(Njettiness::geometric_measure);    _testBeta.push_back(NAN); _testR0.push_back(NAN); _testMeasureNames.push_back("Geometric Measure  (in GeV):");

   
   // When doing Njettiness as a jet algorithm, want to test the cutoff measures.
   // (Since they are not senisible without a cutoff)
   vector<Njettiness::MeasureMode> _testCutoffMeasureModes;
   vector<double> _testCutoffBeta;
   vector<double> _testRcutoff;
   vector<string> _testCutoffMeasureNames;
   _testCutoffMeasureModes.push_back(Njettiness::unnormalized_cutoff_measure); _testCutoffBeta.push_back(1.0); _testRcutoff.push_back(0.8); _testCutoffMeasureNames.push_back("Unnormalized Measure (beta = 1.0, Rcut = 0.8):");
   _testCutoffMeasureModes.push_back(Njettiness::unnormalized_cutoff_measure); _testCutoffBeta.push_back(2.0); _testRcutoff.push_back(0.8); _testCutoffMeasureNames.push_back("Unnormalized Measure (beta = 1.0, Rcut = 0.8):");
   // TODO:  Figure out what to do with Geometric Measure, since order of arguments makes this not work here.
//   _testCutoffMeasureModes.push_back(Njettiness::geometric_cutoff_measure);    _testCutoffBeta.push_back(NAN); _testRcutoff.push_back(0.8); _testMeasureNames.push_back("Geometric Measure (Rcut = 0.8):");

   
   
   
   /////// N-subjettiness /////////////////////////////

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
   
   // small number to show equivalence of doubles
   double epsilon = 0.0001;
   
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
      cout << "Outputting N-subjettiness Values" << endl;
      cout << "-------------------------------------------------------------------------------------" << endl;

      
      // Now loop through all options
      cout << setprecision(6) << right << fixed;
      for (unsigned iM = 0; iM < _testMeasureModes.size(); iM++) {
         
         cout << "-------------------------------------------------------------------------------------" << endl;
         cout << _testMeasureNames[iM] << endl;
         cout << "       AxisMode"
            << setw(14) << "tau1"
            << setw(14) << "tau2"
            << setw(14) << "tau3"
            << setw(14) << "tau2/tau1"
            << setw(14) << "tau3/tau2"
            << endl;
         
         for (unsigned iA = 0; iA < _testAxesModes.size(); iA++) {
            
            // This case doesn't work, so skip it.
            if (_testAxesModes[iA] == Njettiness::min_axes && _testMeasureModes[iM] == Njettiness::geometric_measure) continue;
            
            // define Nsubjettiness functions
            Nsubjettiness        nSub1(1,    _testAxesModes[iA], _testMeasureModes[iM], _testBeta[iM], _testR0[iM]);
            Nsubjettiness        nSub2(2,    _testAxesModes[iA], _testMeasureModes[iM], _testBeta[iM], _testR0[iM]);
            Nsubjettiness        nSub3(3,    _testAxesModes[iA], _testMeasureModes[iM], _testBeta[iM], _testR0[iM]);
            NsubjettinessRatio   nSub21(2,1, _testAxesModes[iA], _testMeasureModes[iM], _testBeta[iM], _testR0[iM]);
            NsubjettinessRatio   nSub32(3,2, _testAxesModes[iA], _testMeasureModes[iM], _testBeta[iM], _testR0[iM]);
            
            // calculate Nsubjettiness values
            double tau1 = nSub1(antikt_jets[j]);
            double tau2 = nSub2(antikt_jets[j]);
            double tau3 = nSub3(antikt_jets[j]);
            double tau21 = nSub21(antikt_jets[j]);
            double tau32 = nSub32(antikt_jets[j]);
            
            // Make sure calculations are consistent
            assert(abs(tau21 - tau2/tau1) < epsilon);
            assert(abs(tau32 - tau3/tau2) < epsilon);
            
            // Output results:
            cout << _testAxesNames[iA]
               << setw(14) << tau1
               << setw(14) << tau2
               << setw(14) << tau3
               << setw(14) << tau21
               << setw(14) << tau32
               << endl;
         }
      }

      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << "Done Outputting N-subjettiness Values" << endl;
      cout << "-------------------------------------------------------------------------------------" << endl;

      
      ////////
      //
      //  Finding axes/jets found by N-subjettiness partitioning
      //
      //  This uses the NjettinessPlugin as a jet finder (in this case, acting on an individual jet)
      //  Recommended usage is same as above
      //
      ///////

      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << "Outputting N-subjettiness Subjets" << endl;
      cout << "-------------------------------------------------------------------------------------" << endl;

      
      // Loop through all options, this time setting up jet finding
      cout << setprecision(6) << left << fixed;
      for (unsigned iM = 0; iM < _testMeasureModes.size(); iM++) {
         
         for (unsigned iA = 0; iA < _testKeyAxesModes.size(); iA++) {

            // This case doesn't work, so skip it.
            if (_testAxesModes[iA] == Njettiness::min_axes && _testMeasureModes[iM] == Njettiness::geometric_measure) continue;
            
            // define the plugins
            NjettinessPlugin nsub_plugin1(1, _testKeyAxesModes[iA], _testMeasureModes[iM], _testBeta[iM], _testR0[iM]);
            NjettinessPlugin nsub_plugin2(2, _testKeyAxesModes[iA], _testMeasureModes[iM], _testBeta[iM], _testR0[iM]);
            NjettinessPlugin nsub_plugin3(3, _testKeyAxesModes[iA], _testMeasureModes[iM], _testBeta[iM], _testR0[iM]);

            // define the corresponding jet definitions
            JetDefinition nsub_jetDef1(&nsub_plugin1);
            JetDefinition nsub_jetDef2(&nsub_plugin2);
            JetDefinition nsub_jetDef3(&nsub_plugin3);

            // and the corresponding cluster sequences
            ClusterSequence nsub_seq1(antikt_jets[j].constituents(), nsub_jetDef1);
            ClusterSequence nsub_seq2(antikt_jets[j].constituents(), nsub_jetDef2);
            ClusterSequence nsub_seq3(antikt_jets[j].constituents(), nsub_jetDef3);

            vector<PseudoJet> jets1 = nsub_seq1.inclusive_jets();
            vector<PseudoJet> jets2 = nsub_seq2.inclusive_jets();
            vector<PseudoJet> jets3 = nsub_seq3.inclusive_jets();

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << _testMeasureNames[iM] << endl;
            cout << _testKeyAxesNamesLong[iA] << endl;
            
            bool commentOut = false;
            if (_testAxesModes[iA] == Njettiness::min_axes) commentOut = true;  // have to comment out min_axes, because it has random values
            
            // This helper function tries to find out if the jets have tau information for printing
            PrintJets(jets1,commentOut);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJets(jets2,commentOut);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJets(jets3,commentOut);
            
            // Also try to find axes location (if njettiness_extras works)
            vector<PseudoJet> axes1;
            vector<PseudoJet> axes2;
            vector<PseudoJet> axes3;
            const NjettinessExtras * extras1 = njettiness_extras(nsub_seq1);
            const NjettinessExtras * extras2 = njettiness_extras(nsub_seq2);
            const NjettinessExtras * extras3 = njettiness_extras(nsub_seq3);

            axes1 = extras1->axes();
            axes2 = extras2->axes();
            axes3 = extras3->axes();
            
            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
            cout << "Axes Used for Above Subjets" << endl;

            PrintJets(axes1,commentOut);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJets(axes2,commentOut);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJets(axes3,commentOut);
            
         }
      }
      
      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << "Done Outputting N-subjettiness Subjets" << endl;
      cout << "-------------------------------------------------------------------------------------" << endl;

   }
   
   
   ////////// N-jettiness as a jet algorithm ///////////////////////////

   ////////
   //
   //  You can also find jets event-wide with Njettiness.
   //  In this case, Winner-Take-All axes are a must, since the other axes get trapped in local minima
   //
   //  Recommended usage of NjettinessPlugin (event-wide)
   //  AxesMode:  wta_kt_axes or onepass_wta_kt_axes
   //  MeasureMode:  unnormalized_measure
   //  beta with wta_kt_axes: anything greater than 0.0 (particularly good for 1.0)
   //  beta with onepass_wta_kt_axes:  between 1.0 and 3.0
   //
   ///////
   
   cout << "-------------------------------------------------------------------------------------" << endl;
   cout << "Using N-jettiness as a Jet Algorithm" << endl;
   cout << "-------------------------------------------------------------------------------------" << endl;

   
   for (unsigned iM = 0; iM < _testCutoffMeasureModes.size(); iM++) {
      
      for (unsigned iA = 0; iA < _testKeyAxesModes.size(); iA++) {
         
         // define the plugins
         NjettinessPlugin njet_plugin2(2, _testKeyAxesModes[iA], _testCutoffMeasureModes[iM], _testBeta[iM], _testRcutoff[iM]);
         NjettinessPlugin njet_plugin3(3, _testKeyAxesModes[iA], _testCutoffMeasureModes[iM], _testBeta[iM], _testRcutoff[iM]);
         NjettinessPlugin njet_plugin4(4, _testKeyAxesModes[iA], _testCutoffMeasureModes[iM], _testBeta[iM], _testRcutoff[iM]);
   
         // and the jet definitions
         JetDefinition njet_jetDef2(&njet_plugin2);
         JetDefinition njet_jetDef3(&njet_plugin3);
         JetDefinition njet_jetDef4(&njet_plugin4);

         // and the cluster sequences
         ClusterSequence njet_seq2(input_particles, njet_jetDef2);
         ClusterSequence njet_seq3(input_particles, njet_jetDef3);
         ClusterSequence njet_seq4(input_particles, njet_jetDef4);

         // and find the jets
         vector<PseudoJet> njet_jets2 = njet_seq2.inclusive_jets();
         vector<PseudoJet> njet_jets3 = njet_seq3.inclusive_jets();
         vector<PseudoJet> njet_jets4 = njet_seq4.inclusive_jets();

         cout << "-------------------------------------------------------------------------------------" << endl;
         cout << _testCutoffMeasureNames[iM] << endl;
         cout << _testKeyAxesNamesLong[iA] << endl;
         
         PrintJets(njet_jets2);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintJets(njet_jets3);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintJets(njet_jets4);

         // The axes might point in a different direction than the jets
         // Using the NjettinessExtras pointer (ClusterSequence::Extras) to access that information
         vector<PseudoJet> njet_axes2;
         vector<PseudoJet> njet_axes3;
         vector<PseudoJet> njet_axes4;
         const NjettinessExtras * extras2 = njettiness_extras(njet_seq2);
         const NjettinessExtras * extras3 = njettiness_extras(njet_seq3);
         const NjettinessExtras * extras4 = njettiness_extras(njet_seq4);
         
         njet_axes2 = extras2->axes();
         njet_axes3 = extras3->axes();
         njet_axes4 = extras4->axes();
         
         cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
         cout << "Axes Used for Above Jets" << endl;
         
         PrintJets(njet_axes2);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintJets(njet_axes3);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintJets(njet_axes4);
         
         bool calculateArea = false;
         if (calculateArea) {
            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
            cout << "Adding Area Information (quite slow)" << endl;
            
            double ghost_maxrap = 5.0; // e.g. if particles go up to y=5
            AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap));
            
            // Defining cluster sequences with area
            ClusterSequenceArea njet_seq_area2(input_particles, njet_jetDef2, area_def);
            ClusterSequenceArea njet_seq_area3(input_particles, njet_jetDef3, area_def);
            ClusterSequenceArea njet_seq_area4(input_particles, njet_jetDef4, area_def);
            
            vector<PseudoJet> njet_jets_area2 = njet_seq_area2.inclusive_jets();
            vector<PseudoJet> njet_jets_area3 = njet_seq_area3.inclusive_jets();
            vector<PseudoJet> njet_jets_area4 = njet_seq_area4.inclusive_jets();
            
            PrintJets(njet_jets_area2);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJets(njet_jets_area3);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJets(njet_jets_area4);
         }
         
      }
   }
   
   cout << "-------------------------------------------------------------------------------------" << endl;
   cout << "Done Using N-jettiness as a Jet Algorithm" << endl;
   cout << "-------------------------------------------------------------------------------------" << endl;

}


void PrintJets(const vector <PseudoJet>& jets, bool commentOut) {
   
   string commentStr = "";
   if (commentOut) commentStr = "#";
   
   // gets extras information
   if (jets.size() == 0) return;
   const NjettinessExtras * extras = njettiness_extras(jets[0]);
   
   bool useExtras = (extras != NULL);
   bool useArea = jets[0].has_area();
   
   // define nice tauN header
   int N = jets.size();
   stringstream ss(""); ss << "tau" << N; string tauName = ss.str();
   
   cout << fixed << right;
   
   cout << commentStr << setw(5) << "jet #" << "   "
   <<  setw(10) << "rap"
   <<  setw(10) << "phi"
   <<  setw(11) << "pt"
   <<  setw(11) << "m"
   <<  setw(11) << "e";
   if (useExtras) cout << setw(14) << tauName;
   if (useArea) cout << setw(10) << "area";
   cout << endl;
   
   fastjet::PseudoJet total(0,0,0,0);
   
   // print out individual jet information
   for (unsigned i = 0; i < jets.size(); i++) {
      cout << commentStr << setw(5) << i+1  << "   "
      << setprecision(4) <<  setw(10) << jets[i].rap()
      << setprecision(4) <<  setw(10) << jets[i].phi()
      << setprecision(4) <<  setw(11) << jets[i].perp()
      << setprecision(4) <<  setw(11) << max(jets[i].m(),0.0) // needed to fix -0.0 issue on some compilers.
      << setprecision(4) <<  setw(11) << jets[i].e();
      if (useExtras) cout << setprecision(6) <<  setw(14) << max(extras->subTau(jets[i]),0.0);
      if (useArea) cout << setprecision(4) << setw(10) << (jets[i].has_area() ? jets[i].area() : 0.0 );
      cout << endl;
      total += jets[i];
   }
   
   // print out total jet
   if (useExtras) {
      double beamTau = extras->beamTau();
      
      if (beamTau > 0.0) {
         cout << commentStr << setw(5) << " beam" << "   "
         <<  setw(10) << ""
         <<  setw(10) << ""
         <<  setw(11) << ""
         <<  setw(11) << ""
         <<  setw(11) << ""
         <<  setw(14) << setprecision(6) << beamTau
         << endl;
      }
      
      cout << commentStr << setw(5) << "total" << "   "
      <<  setprecision(4) << setw(10) << total.rap()
      <<  setprecision(4) << setw(10) << total.phi()
      <<  setprecision(4) << setw(11) << total.perp()
      <<  setprecision(4) << setw(11) << max(total.m(),0.0) // needed to fix -0.0 issue on some compilers.
      <<  setprecision(4) <<  setw(11) << total.e()
      <<  setprecision(6) << setw(14) << extras->totalTau();
      if (useArea) cout << setprecision(4) << setw(10) << (total.has_area() ? total.area() : 0.0);
      cout << endl;
   }
   
}


