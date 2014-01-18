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

/// Helper function for output (bool used to specify if axes should be commented in the final result)
void PrintJets(const vector <PseudoJet>& jets, bool min_axes = false) {

   if (jets.size() == 0) return;
   const NjettinessExtras * extras = njettiness_extras(jets[0]);

   if (jets[0].has_area()) {
      if (extras == NULL) {
         if (min_axes) printf("#%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity", "phi", "pt","m","e","area"); // label the columns
         else printf("%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity", "phi", "pt","m","e","area"); // label the columns
         for (unsigned int i = 0; i < jets.size(); i++) {
            if (min_axes) printf("#%5u %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),max(jets[i].m(),0.0),jets[i].e(),(jets[i].has_area() ? jets[i].area() : 0.0 ));
            else printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),max(jets[i].m(),0.0),jets[i].e(),(jets[i].has_area() ? jets[i].area() : 0.0 ));
         }
      }
      else {
         fastjet::PseudoJet total(0,0,0,0);
         if (min_axes) printf("#%5s %10s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity", "phi", "pt","m","e","subTau","area"); // label the columns
         else printf("%5s %10s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity", "phi", "pt","m","e","subTau","area"); // label the columns
         for (unsigned int i = 0; i < jets.size(); i++) {
            if (min_axes) printf("#%5u %10.3f %10.3f %10.3f %10.3f %10.3f %10.6f %10.3f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),max(jets[i].m(),0.0),jets[i].e(),extras->subTau(jets[i]),(jets[i].has_area() ? jets[i].area() : 0.0 ));
            else printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f %10.6f %10.3f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),max(jets[i].m(),0.0),jets[i].e(),extras->subTau(jets[i]),(jets[i].has_area() ? jets[i].area() : 0.0 ));
            total += jets[i];
         }   
         if (min_axes) printf("%5s %10.3f %10.3f %10.3f %10.3f %10.3f %10.6f %10.3f\n","total", total.rap(), total.phi(), total.perp(),total.m(),total.e(),extras->totalTau(),(total.has_area() ? total.area() : 0.0 ));
         else printf("%5s %10.3f %10.3f %10.3f %10.3f %10.3f %10.6f %10.3f\n","total", total.rap(), total.phi(), total.perp(),total.m(),total.e(),extras->totalTau(),(total.has_area() ? total.area() : 0.0 ));
      }
   } else {
      if (extras == NULL) {
         if (min_axes) printf("#%5s %10s %10s %10s %10s %10s\n","jet #", "rapidity", "phi", "pt","m","e"); // label the columns
         else printf("%5s %10s %10s %10s %10s %10s\n","jet #", "rapidity", "phi", "pt","m","e"); // label the columns
         for (unsigned int i = 0; i < jets.size(); i++) {
            if (min_axes) printf("#%5u %10.3f %10.3f %10.3f %10.3f %10.3f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),max(jets[i].m(),0.0),jets[i].e());
            else printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),max(jets[i].m(),0.0),jets[i].e());
         }
      }
      else {
         fastjet::PseudoJet total(0,0,0,0);
         if (min_axes) printf("#%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity", "phi", "pt","m","e","subTau"); // label the columns
         else printf("%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity", "phi", "pt","m","e","subTau"); // label the columns
         for (unsigned int i = 0; i < jets.size(); i++) {
            if (min_axes) printf("#%5u %10.3f %10.3f %10.3f %10.3f %10.3f %10.6f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),max(jets[i].m(),0.0),jets[i].e(),extras->subTau(jets[i]));
            else printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f %10.6f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),max(jets[i].m(),0.0),jets[i].e(),extras->subTau(jets[i]));
            total += jets[i];
         }   
         if (min_axes) printf("%5s %10.3f %10.3f %10.3f %10.3f %10.3f %10.6f\n","total", total.rap(), total.phi(), total.perp(),total.m(),total.e(),extras->totalTau());
         else printf("%5s %10.3f %10.3f %10.3f %10.3f %10.3f %10.6f\n","total", total.rap(), total.phi(), total.perp(),total.m(),total.e(),extras->totalTau());
      }
   }

}

////////
//
//  Main Routine for Analysis 
//
///////

void analyze(const vector<PseudoJet> & input_particles) {

   /////// N-subjettiness /////////////////////////////
   
   // Initial clustering with anti-kt algorithm
   JetAlgorithm algorithm = antikt_algorithm; 
   double jet_rad = 1.00; // jet radius for anti-kt algorithm
   JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
   ClusterSequence clust_seq(input_particles,jetDef);
   vector<PseudoJet> antikt_jets  = sorted_by_pt(clust_seq.inclusive_jets());
   
   // Defining Nsubjettiness parameters
   double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
   double R0 = 1.0; // Characteristic jet radius for normalization	      
   //double Rcut = 1.0; // maximum R particles can be from axis to be included in jet
   
   // used to choose whether or not min_axes are used
   bool test_min_axes = true;
   // small number to show equivalence of doubles
   double epsilon = 0.0001;

   for (int j = 0; j < 2; j++) { // Two hardest jets per event
      if (antikt_jets[j].perp() > 200) {
         vector<PseudoJet> jet_constituents = clust_seq.constituents(antikt_jets[j]);
         
         //
         // If you don't want subjets, you can use the simple functor Nsubjettiness:
         // Recommended usage is Njettiness::onepass_kt_axes mode.
         //
         
         //
         // Note:  all instances of axes mode Njettiness::min_axes will only be used 
         // with the bool test_min_axes. Otherwise the min values are set to 0.  
         // This method is not guarenteed to give a global
         // minimum, only a local minimum, and different choices of the random
         // number seed can give different results.  For that reason,
         // Njettiness::onepass_kt_axes is the recommended usage.
         //

         //
         // Normalized Measure Values (kT, min, onepass kT, Winner take all kT)
         //

         // 1-subjettiness
         Nsubjettiness nSub1KT(1, Njettiness::kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau1 = nSub1KT(antikt_jets[j]);
         Nsubjettiness nSub1OnePass(1, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau1onepass = nSub1OnePass(antikt_jets[j]);
         Nsubjettiness nSub1WTAKT(1, Njettiness::wta_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau1WTA = nSub1WTAKT(antikt_jets[j]);

         // 2-subjettiness
         Nsubjettiness nSub2KT(2, Njettiness::kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau2 = nSub2KT(antikt_jets[j]);
         Nsubjettiness nSub2OnePass(2, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau2onepass = nSub2OnePass(antikt_jets[j]);
         Nsubjettiness nSub2WTAKT(2, Njettiness::wta_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau2WTA = nSub2WTAKT(antikt_jets[j]);

         // 3-subjettiness
         Nsubjettiness nSub3KT(3, Njettiness::kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau3 = nSub3KT(antikt_jets[j]);
         Nsubjettiness nSub3OnePass(3, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau3onepass = nSub3OnePass(antikt_jets[j]);
         Nsubjettiness nSub3WTAKT(3, Njettiness::wta_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau3WTA = nSub3WTAKT(antikt_jets[j]);

         //
         // User can also calculate ratios using new NsubjettinessRatio class
         //

         NsubjettinessRatio nSub3over2KT(3, 2, Njettiness::kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau3over2KT = nSub3over2KT(antikt_jets[j]);
         NsubjettinessRatio nSub2over1KT(2, 1, Njettiness::kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau2over1KT = nSub2over1KT(antikt_jets[j]);

         NsubjettinessRatio nSub3over2OnePass(3, 2, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau3over2onepass = nSub3over2OnePass(antikt_jets[j]);
         NsubjettinessRatio nSub2over1OnePass(2, 1, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau2over1onepass = nSub2over1OnePass(antikt_jets[j]);

         NsubjettinessRatio nSub3over2WTAKT(3, 2, Njettiness::wta_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau3over2WTA = nSub3over2WTAKT(antikt_jets[j]);
         NsubjettinessRatio nSub2over1WTAKT(2, 1, Njettiness::wta_kt_axes, Njettiness::normalized_measure, beta, R0);
         double tau2over1WTA = nSub2over1WTAKT(antikt_jets[j]);

         //
         //  Unnormalized Measure Values (kT, min, onepass kT, Winner take all kT)
         //

         // 1-subjettiness
         Nsubjettiness nSub1KT_unnorm(1, Njettiness::kt_axes, Njettiness::unnormalized_measure, beta);
         double tau1_unnorm = nSub1KT_unnorm(antikt_jets[j]);
         Nsubjettiness nSub1OnePass_unnorm(1, Njettiness::onepass_kt_axes, Njettiness::unnormalized_measure, beta);
         double tau1onepass_unnorm = nSub1OnePass_unnorm(antikt_jets[j]);
         Nsubjettiness nSub1WTAKT_unnorm(1, Njettiness::wta_kt_axes, Njettiness::unnormalized_measure, beta);
         double tau1WTA_unnorm = nSub1WTAKT_unnorm(antikt_jets[j]);

         // 2-subjettiness
         Nsubjettiness nSub2KT_unnorm(2, Njettiness::kt_axes, Njettiness::unnormalized_measure, beta);
         double tau2_unnorm = nSub2KT_unnorm(antikt_jets[j]);
         Nsubjettiness nSub2OnePass_unnorm(2, Njettiness::onepass_kt_axes, Njettiness::unnormalized_measure, beta);
         double tau2onepass_unnorm = nSub2OnePass_unnorm(antikt_jets[j]);
         Nsubjettiness nSub2WTAKT_unnorm(2, Njettiness::wta_kt_axes, Njettiness::unnormalized_measure, beta);
         double tau2WTA_unnorm = nSub2WTAKT_unnorm(antikt_jets[j]);

         // 3-subjettiness
         Nsubjettiness nSub3KT_unnorm(3, Njettiness::kt_axes, Njettiness::unnormalized_measure, beta);
         double tau3_unnorm = nSub3KT_unnorm(antikt_jets[j]);
         Nsubjettiness nSub3OnePass_unnorm(3, Njettiness::onepass_kt_axes, Njettiness::unnormalized_measure, beta);
         double tau3onepass_unnorm = nSub3OnePass_unnorm(antikt_jets[j]);
         Nsubjettiness nSub3WTAKT_unnorm(3, Njettiness::wta_kt_axes, Njettiness::unnormalized_measure, beta);
         double tau3WTA_unnorm = nSub3WTAKT_unnorm(antikt_jets[j]);
         
         // min_axes currently in testing, all information defined in here
         double tau1min, tau2min, tau3min;
         double tau1min_unnorm, tau2min_unnorm, tau3min_unnorm;
         double tau3over2min, tau2over1min;
         if (test_min_axes) {  
            Nsubjettiness nSub1Min(1, Njettiness::min_axes, Njettiness::normalized_measure, beta, R0);
            tau1min = nSub1Min(antikt_jets[j]);

            Nsubjettiness nSub2Min(2, Njettiness::min_axes, Njettiness::normalized_measure, beta, R0);
            tau2min = nSub2Min(antikt_jets[j]);

            Nsubjettiness nSub3Min(3, Njettiness::min_axes, Njettiness::normalized_measure, beta, R0);
            tau3min = nSub3Min(antikt_jets[j]);

            Nsubjettiness nSub1Min_unnorm(1, Njettiness::min_axes, Njettiness::unnormalized_measure, beta);
            tau1min_unnorm = nSub1Min_unnorm(antikt_jets[j]);

            Nsubjettiness nSub2Min_unnorm(2, Njettiness::min_axes, Njettiness::unnormalized_measure, beta);
            tau2min_unnorm = nSub2Min_unnorm(antikt_jets[j]);

            Nsubjettiness nSub3Min_unnorm(3, Njettiness::min_axes, Njettiness::unnormalized_measure, beta);
            tau3min_unnorm = nSub3Min_unnorm(antikt_jets[j]);

            NsubjettinessRatio nSub3over2min(3, 2, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, beta, R0);
            tau3over2min = nSub3over2min(antikt_jets[j]);

            NsubjettinessRatio nSub2over1min(2, 1, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, beta, R0);
            tau2over1min = nSub2over1min(antikt_jets[j]);
         }
         
         // Show the ratio values are the same (use epsilon to properly check double equality): 

         assert(tau3over2KT - tau3/tau2 < epsilon);
         assert(tau2over1KT - tau2/tau1 < epsilon);
         assert(tau3over2onepass - tau3onepass/tau2onepass < epsilon);
         assert(tau2over1onepass - tau2onepass/tau1onepass < epsilon);
         assert(tau3over2WTA - tau3WTA/tau2WTA < epsilon);
         assert(tau2over1WTA - tau2WTA/tau1WTA < epsilon);

         //
         // Or, if you want subjets, use the FastJet plugin on a jet's constituents
         // Recommended usage is onepass_kt_axes or wta_kt_axes mode.
         //
         
         NjettinessPlugin nsub_plugin1(1, Njettiness::kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsub_jetDef1(&nsub_plugin1);
         ClusterSequence nsub_seq1(antikt_jets[j].constituents(), nsub_jetDef1);
         vector<PseudoJet> kt1jets = nsub_seq1.inclusive_jets();
         
         NjettinessPlugin nsub_plugin2(2, Njettiness::kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsub_jetDef2(&nsub_plugin2);
         ClusterSequence nsub_seq2(antikt_jets[j].constituents(), nsub_jetDef2);
         vector<PseudoJet> kt2jets = nsub_seq2.inclusive_jets();

         NjettinessPlugin nsub_plugin3(3, Njettiness::kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsub_jetDef3(&nsub_plugin3);
         ClusterSequence nsub_seq3(antikt_jets[j].constituents(), nsub_jetDef3);
         vector<PseudoJet> kt3jets = nsub_seq3.inclusive_jets();

         NjettinessPlugin nsubOnePass_plugin1(1, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsubOnePass_jetDef1(&nsubOnePass_plugin1);
         ClusterSequence nsubOnePass_seq1(antikt_jets[j].constituents(), nsubOnePass_jetDef1);
         vector<PseudoJet> onepass1jets = nsubOnePass_seq1.inclusive_jets();

         NjettinessPlugin nsubOnePass_plugin2(2, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsubOnePass_jetDef2(&nsubOnePass_plugin2);
         ClusterSequence nsubOnePass_seq2(antikt_jets[j].constituents(), nsubOnePass_jetDef2);
         vector<PseudoJet> onepass2jets = nsubOnePass_seq2.inclusive_jets();
         
         NjettinessPlugin nsubOnePass_plugin3(3, Njettiness::onepass_kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsubOnePass_jetDef3(&nsubOnePass_plugin3);
         ClusterSequence nsubOnePass_seq3(antikt_jets[j].constituents(), nsubOnePass_jetDef3);
         vector<PseudoJet> onepass3jets = nsubOnePass_seq3.inclusive_jets();
                 
         NjettinessPlugin nsubWTA_plugin1(1, Njettiness::wta_kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsubWTA_jetDef1(&nsubWTA_plugin1);
         ClusterSequence nsubWTA_seq1(antikt_jets[j].constituents(), nsubWTA_jetDef1);
         vector<PseudoJet> wta1jets = nsubWTA_seq1.inclusive_jets();

         NjettinessPlugin nsubWTA_plugin2(2, Njettiness::wta_kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsubWTA_jetDef2(&nsubWTA_plugin2);
         ClusterSequence nsubWTA_seq2(antikt_jets[j].constituents(), nsubWTA_jetDef2);
         vector<PseudoJet> wta2jets = nsubWTA_seq2.inclusive_jets();
         
         NjettinessPlugin nsubWTA_plugin3(3, Njettiness::wta_kt_axes, Njettiness::normalized_measure, 1.0, 1.0);
         JetDefinition nsubWTA_jetDef3(&nsubWTA_plugin3);
         ClusterSequence nsubWTA_seq3(antikt_jets[j].constituents(), nsubWTA_jetDef3);
         vector<PseudoJet> wta3jets = nsubWTA_seq3.inclusive_jets();

         // min_axes currently in testing
         vector<PseudoJet> min1jets, min2jets, min3jets;
         if (test_min_axes) {        
            NjettinessPlugin nsubMin_plugin1(1, Njettiness::min_axes, 1.0, 1.0, 1.0);
            JetDefinition nsubMin_jetDef1(&nsubMin_plugin1);
            ClusterSequence nsubMin_seq1(antikt_jets[j].constituents(), nsubMin_jetDef1);
            min1jets = nsubMin_seq1.inclusive_jets();
         
            NjettinessPlugin nsubMin_plugin2(2, Njettiness::min_axes, 1.0, 1.0, 1.0);
            JetDefinition nsubMin_jetDef2(&nsubMin_plugin2);
            ClusterSequence nsubMin_seq2(antikt_jets[j].constituents(), nsubMin_jetDef2);
            min2jets = nsubMin_seq2.inclusive_jets();

            NjettinessPlugin nsubMin_plugin3(3, Njettiness::min_axes, 1.0, 1.0, 1.0);
            JetDefinition nsubMin_jetDef3(&nsubMin_plugin3);
            ClusterSequence nsubMin_seq3(antikt_jets[j].constituents(), nsubMin_jetDef3);
            min3jets = nsubMin_seq3.inclusive_jets();
         }
         else {
            min1jets.resize(0);
            min2jets.resize(0);
            min3jets.resize(0);
         }

         cout << "Jet " << j + 1 << ":" << endl;                 
         printf("-------------------------------------------------------------------------------------"); printf("\n");
         printf("-------------------------------------------------------------------------------------"); printf("\n");
         cout << "Beta = " << beta << endl;
         cout << "kT Axes:" << endl;
         PrintJets(kt1jets);
         PrintJets(kt2jets);
         PrintJets(kt3jets);
         cout << "One Pass Minimization Axes from kT" << endl;
         PrintJets(onepass1jets);
         PrintJets(onepass2jets);
         PrintJets(onepass3jets);            
         cout << "Winner Take All Axes with kT" << endl;
         PrintJets(wta1jets);
         PrintJets(wta2jets);
         PrintJets(wta3jets);            
         // currently in testing (bool set so that it is commented in output file)
         if (test_min_axes) {
         cout << "#Multi-Pass Minimization Axes:" << endl;
            PrintJets(min1jets, true);
            PrintJets(min2jets, true);
            PrintJets(min3jets, true);
         }

         printf("-------------------------------------------------------------------------------------"); printf("\n");
         cout << "Beta = " << beta << setprecision(6) << endl;
         cout << "Normalized values:" << endl;
         cout << "     kT: " << "tau1: " << tau1 << "  tau2: " << tau2 << "  tau3: " << tau3 << "  tau2/tau1: " << tau2over1KT << "  tau3/tau2: " << tau3over2KT << endl;
         cout << "     OnePass: " << "tau1: " << tau1onepass << "  tau2: " << tau2onepass << "  tau3: " << tau3onepass << "  tau2/tau1: " << tau2over1onepass << "  tau3/tau2: " << tau3over2onepass << endl;
         cout << "     WTA kT: " << "tau1: " << tau1WTA << "  tau2: " << tau2WTA << "  tau3: " << tau3WTA << "  tau2/tau1: " << tau2over1WTA << "  tau3/tau2: " << tau3over2WTA << endl;
         if (test_min_axes) cout << "#    Min: " << "tau1: " << tau1min << "  tau2: " << tau2min << "  tau3: " << tau3min << "  tau2/tau1: " << tau2over1min << "  tau3/tau2: " << tau3over2min << endl;
         cout << endl;
         cout << "Unnormalized values:" << endl;
         cout << "     kT: " << "tau1: " << tau1_unnorm << "  tau2: " << tau2_unnorm << "  tau3: " << tau3_unnorm << "  tau2/tau1: " << tau2over1KT << "  tau3/tau2: " << tau3over2KT << endl;
         cout << "     OnePass: " << "tau1: " << tau1onepass_unnorm << "  tau2: " << tau2onepass_unnorm << "  tau3: " << tau3onepass_unnorm << "  tau2/tau1: " << tau2over1onepass << "  tau3/tau2: " << tau3over2onepass << endl;
         cout << "     WTA kT: " << "tau1: " << tau1WTA_unnorm << "  tau2: " << tau2WTA_unnorm << "  tau3: " << tau3WTA_unnorm << "  tau2/tau1: " << tau2over1WTA << "  tau3/tau2: " << tau3over2WTA << endl;
         if (test_min_axes) cout << "#    Min: " << "tau1: " << tau1min_unnorm << "  tau2: " << tau2min_unnorm << "  tau3: " << tau3min_unnorm << "  tau2/tau1: " << tau2over1min << "  tau3/tau2: " << tau3over2min << endl;
         cout << endl;
         printf("-------------------------------------------------------------------------------------"); printf("\n");
         printf("-------------------------------------------------------------------------------------"); printf("\n");
      }
   }
   
   
   ////////// N-jettiness as a jet algorithm ///////////////////////////

   // You can also find jets with Njettiness:
   
   //  Using Winner-Take-All axes and standard measure
   NjettinessPlugin njet_plugin(3, Njettiness::wta_kt_axes, Njettiness::unnormalized_cutoff_measure, 1.0, 1.0);
   JetDefinition njet_jetDef(&njet_plugin);
   ClusterSequence njet_seq(input_particles, njet_jetDef);
   vector<PseudoJet> njet_jets = njet_seq.inclusive_jets();

   // Using WTA but geometric measure
   NjettinessPlugin geo_plugin(3, Njettiness::wta_kt_axes, Njettiness::geometric_cutoff_measure, 1.0);
   JetDefinition geo_jetDef(&geo_plugin);
   ClusterSequence geo_seq(input_particles, geo_jetDef);
   vector<PseudoJet> geo_jets = geo_seq.inclusive_jets();
      
   // The axes might point in a different direction than the jets
   // Using the NjettinessExtras pointer (ClusterSequence::Extras) to access that information
   vector<PseudoJet> njet_axes;
   const NjettinessExtras * extras = njettiness_extras(njet_seq);
   if (extras != NULL) {
      njet_axes = extras->axes();
   }
   
   printf("-------------------------------------------------------------------------------------"); printf("\n");
   // cout << "Event-wide Jets from One-Pass Minimization (beta = 1.0)" << endl;
   cout << "Event-wide Jets from Winner-take-all kT Axes (beta = 1.0)" << endl;
   PrintJets(njet_jets);
   cout << "Event-wide Axis Location for Above Jets" << endl;
   PrintJets(njet_axes);
   cout << "Event-wide Jets from Geometric Measure" << endl;
   PrintJets(geo_jets);
   printf("-------------------------------------------------------------------------------------"); printf("\n");

   // You can also find jet areas using this method (quite slow, though)

   double ghost_maxrap = 5.0; // e.g. if particles go up to y=5
   AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap));
   
   ClusterSequenceArea njet_seq_area(input_particles, njet_jetDef,area_def);
   vector<PseudoJet> njet_jets_area = njet_seq_area.inclusive_jets();

   ClusterSequenceArea geo_seq_area(input_particles, geo_jetDef,area_def);
   vector<PseudoJet> geo_jets_area = geo_seq_area.inclusive_jets();
      
   // The axes might point in a different direction than the jets
   // Using the NjettinessExtras pointer (ClusterSequence::Extras) to access that information
   vector<PseudoJet> njet_axes_area;
   const NjettinessExtras * extras_area = njettiness_extras(njet_seq_area);
   if (extras_area != NULL) {
      njet_axes_area = extras_area->axes();
   }
   
   printf("-------------------------------------------------------------------------------------"); printf("\n");
   // cout << "Event-wide Jets from One-Pass Minimization (beta = 1.0) (with area information)" << endl;
   cout << "Event-wide Jets from Winner-take-all axes (beta = 1.0) (with area information)" << endl;
   PrintJets(njet_jets_area);
   cout << "Event-wide Axis Location for Above Jets (with area information)" << endl;
   PrintJets(njet_axes_area);
   cout << "Event-wide Jets from Geometric Measure (with area information)" << endl;
   PrintJets(geo_jets_area);
   printf("-------------------------------------------------------------------------------------"); printf("\n");

}


