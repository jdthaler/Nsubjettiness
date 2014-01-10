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

//NEW FILE CREATED BY TJW 12/25
//Update to move AxesFinder class and definitions into separate .cc/.hh files

#include "AxesFinder.hh" //new .hh file added by TJW 12/25

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Functions for minimization.
//
///////

// Minimization function definitions moved from Njettiness.hh -- TJW 12/22

// Given starting axes, update to find better axes by using Kmeans clustering around the old axes --comment added by TJW
//updated function arguments to use three separate parameters instead of NsubParameters-- TJW 1/9                                  
template <int N>
std::vector<LightLikeAxis> AxesFinderFromKmeansMinimization::UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes, 
                                  const std::vector <fastjet::PseudoJet> & inputJets,
                                  double beta, double R0, double Rcutoff, double precision) {
   assert(old_axes.size() == N);
   
   // some storage, declared static to save allocation/re-allocation costs
   static LightLikeAxis new_axes[N];
   static fastjet::PseudoJet new_jets[N];
   for (int n = 0; n < N; ++n) {
      new_axes[n].reset(0.0,0.0,0.0,0.0);
#ifdef FASTJET2
      new_jets[n].reset(0.0,0.0,0.0,0.0);
#else
      // use cheaper reset if available
      new_jets[n].reset_momentum(0.0,0.0,0.0,0.0);
#endif
   }

   // no longer necessary since beta and Rcutoff are defined in arguments -- TJW 1/9
   // double beta = paraNsub.beta();
   // double Rcutoff = paraNsub.Rcutoff();
   
   /////////////// Assignment Step //////////////////////////////////////////////////////////
   std::vector<int> assignment_index(inputJets.size()); 
   int k_assign = -1;
   
   for (unsigned i = 0; i < inputJets.size(); i++){
      double smallestDist = 1000000.0;
      for (int k = 0; k < N; k++) {
         double thisDist = old_axes[k].DistanceSq(inputJets[i]);
         if (thisDist < smallestDist) {
            smallestDist = thisDist;
            k_assign = k;
         }
      }
      if (smallestDist > sq(Rcutoff)) {k_assign = -1;}
      assignment_index[i] = k_assign;
   }
   
   //////////////// Update Step /////////////////////////////////////////////////////////////
   double distPhi, old_dist;
   for (unsigned i = 0; i < inputJets.size(); i++) {
      int old_jet_i = assignment_index[i];
      if (old_jet_i == -1) {continue;}

      const fastjet::PseudoJet& inputJet_i = inputJets[i];
      LightLikeAxis& new_axis_i = new_axes[old_jet_i];
      double inputPhi_i = inputJet_i.phi();
      double inputRap_i = inputJet_i.rap();
            
      // optimize pow() call
      // add noise (the precision term) to make sure we don't divide by zero
      if (beta == 1.0) {
         double DR = std::sqrt(sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i));
         old_dist = 1.0/DR;
      } else if (beta == 2.0) {
         old_dist = 1.0;
      } else if (beta == 0.0) {
         double DRSq = sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i);
         old_dist = 1.0/DRSq;
      } else {
         old_dist = sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i);
         old_dist = std::pow(old_dist, (0.5*beta-1.0));
      }
      
      // TODO:  Put some of these addition functions into light-like axes
      // rapidity sum
      new_axis_i.set_rap(new_axis_i.rap() + inputJet_i.perp() * inputRap_i * old_dist);
      // phi sum
      distPhi = inputPhi_i - old_axes[old_jet_i].phi();
      if (fabs(distPhi) <= M_PI){
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * inputPhi_i * old_dist );
      } else if (distPhi > M_PI) {
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * (-2*M_PI + inputPhi_i) * old_dist );
      } else if (distPhi < -M_PI) {
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * (+2*M_PI + inputPhi_i) * old_dist );
      }
      // weights sum
      new_axis_i.set_weight( new_axis_i.weight() + inputJet_i.perp() * old_dist );
      // momentum magnitude sum
      new_jets[old_jet_i] += inputJet_i;
   }
   // normalize sums
   for (int k = 0; k < N; k++) {
      if (new_axes[k].weight() == 0) {
         // no particles were closest to this axis!  Return to old axis instead of (0,0,0,0)
         new_axes[k] = old_axes[k];
      } else {
         new_axes[k].set_rap( new_axes[k].rap() / new_axes[k].weight() );
         new_axes[k].set_phi( new_axes[k].phi() / new_axes[k].weight() );
         new_axes[k].set_phi( std::fmod(new_axes[k].phi() + 2*M_PI, 2*M_PI) );
         new_axes[k].set_mom( std::sqrt(new_jets[k].modp2()) );
      }
   }
   std::vector<LightLikeAxis> new_axes_vec(N);
   for (unsigned k = 0; k < N; ++k) new_axes_vec[k] = new_axes[k];
   return new_axes_vec;
}

// Given N starting axes, this function updates all axes to find N better axes. 
// (This is just a wrapper for the templated version above.)
//updated function arguments to use three separate parameters instead of NsubParameters-- TJW 1/9                                  
std::vector<LightLikeAxis> AxesFinderFromKmeansMinimization::UpdateAxes(const std::vector <LightLikeAxis> & old_axes, 
                                      const std::vector <fastjet::PseudoJet> & inputJets, double beta, double R0, double Rcutoff, double precision) {
   int N = old_axes.size();
   switch (N) {
      case 1: return UpdateAxesFast<1>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 2: return UpdateAxesFast<2>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 3: return UpdateAxesFast<3>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 4: return UpdateAxesFast<4>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 5: return UpdateAxesFast<5>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 6: return UpdateAxesFast<6>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 7: return UpdateAxesFast<7>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 8: return UpdateAxesFast<8>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 9: return UpdateAxesFast<9>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 10: return UpdateAxesFast<10>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 11: return UpdateAxesFast<11>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 12: return UpdateAxesFast<12>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 13: return UpdateAxesFast<13>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 14: return UpdateAxesFast<14>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 15: return UpdateAxesFast<15>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 16: return UpdateAxesFast<16>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 17: return UpdateAxesFast<17>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 18: return UpdateAxesFast<18>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 19: return UpdateAxesFast<19>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      case 20: return UpdateAxesFast<20>(old_axes, inputJets, beta, R0, Rcutoff, precision);
      default: std::cout << "N-jettiness is hard-coded to only allow up to 20 jets!" << std::endl;
         return std::vector<LightLikeAxis>();
   }

}

// Go from internal LightLikeAxis to PseudoJet
// TODO:  Make part of LightLikeAxis class.
std::vector<fastjet::PseudoJet> ConvertToPseudoJet(const std::vector <LightLikeAxis>& axes) {
   
   int n_jets = axes.size();
   
   double px, py, pz, E;
   std::vector<fastjet::PseudoJet> FourVecJets;
   for (int k = 0; k < n_jets; k++) {
      E = axes[k].mom();
      pz = (std::exp(2.0*axes[k].rap()) - 1.0) / (std::exp(2.0*axes[k].rap()) + 1.0) * E;
      px = std::cos(axes[k].phi()) * std::sqrt( std::pow(E,2) - std::pow(pz,2) );
      py = std::sin(axes[k].phi()) * std::sqrt( std::pow(E,2) - std::pow(pz,2) );
      fastjet::PseudoJet temp = fastjet::PseudoJet(px,py,pz,E);
      FourVecJets.push_back(temp);
   }
   return FourVecJets;
}

// GetMinimumAxes replaced with getAxes, first 4 lines declare instances of parameters to match previous definition -- TJW 12/28  

// uses minimization of N-jettiness to continually update axes for as many times as the user specifies, 
// or until convergence. The function returns the axes found at the minimum -- comment added by TJW
std::vector<fastjet::PseudoJet> AxesFinderFromKmeansMinimization::getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& currentAxes) {
	
	std::vector<fastjet::PseudoJet> seedAxes = _startingFinder->getAxes(n_jets, inputJets, currentAxes);   
	KmeansParameters para = _paraKmeans;
  // parameter definitions added and paraNsub definition removed -- TJW 1/9
	double beta = _beta;
  double R0 = _R0;
  double Rcutoff = _Rcutoff;
  MeasureFunction* functor = _functor;
    double noise = 0, tau = 10000.0, tau_tmp, cmp;
    std::vector< LightLikeAxis > new_axes(n_jets, LightLikeAxis(0,0,0,0)), old_axes(n_jets, LightLikeAxis(0,0,0,0));
    std::vector<fastjet::PseudoJet> tmp_min_axes, min_axes;
    for (int l = 0; l < para.n_iterations(); l++) { // Do minimization procedure multiple times
      // Add noise to guess for the axes
      for (int k = 0; k < n_jets; k++) {
         if ( 0 == l ) { // Don't add noise on first pass
            old_axes[k].set_rap( seedAxes[k].rap() );
            old_axes[k].set_phi( seedAxes[k].phi() );
         } else {
            noise = ((double)rand()/(double)RAND_MAX) * para.noise_range() * 2 - para.noise_range();
            old_axes[k].set_rap( seedAxes[k].rap() + noise );
            noise = ((double)rand()/(double)RAND_MAX) * para.noise_range() * 2 - para.noise_range();
            old_axes[k].set_phi( seedAxes[k].phi() + noise );
         }
      }
      cmp = 100.0; int h = 0;
      while (cmp > para.precision() && h < para.halt()) { // Keep updating axes until near-convergence or too many update steps
         cmp = 0.0; h++;
         new_axes = UpdateAxes(old_axes, inputJets, beta, R0, Rcutoff, para.precision()); // Update axes //updated to use separate parameters instead of NsubParameters -- TJW 1/9
         for (int k = 0; k < n_jets; k++) {
            cmp += old_axes[k].Distance(new_axes[k]);
         }
         cmp = cmp / ((double) n_jets);
         old_axes = new_axes;
      }
      tmp_min_axes = ConvertToPseudoJet(old_axes); // Convert axes directions into four-std::vectors

      tau_tmp = functor->tau(inputJets, tmp_min_axes); 
      if (tau_tmp < tau) {tau = tau_tmp; min_axes = tmp_min_axes;} // Keep axes and tau only if they are best so far
   }	
   return min_axes;
}

//definition of getAxes moved from AxesFinder.hh -- TJW 12/28

// uses minimization of the geometric distance in order to find the minimum axes. It continually updates until it reaches convergence
// or it reaches the maximum number of attempts. -- comment added by TJW 
std::vector<fastjet::PseudoJet> AxesFinderFromGeometricMinimization::getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& currentAxes) {

    std::vector<fastjet::PseudoJet> seedAxes = _startingFinder->getAxes(n_jets, particles, currentAxes);
    double seedTau = _functor->tau(particles,seedAxes);
         
    for (int i = 0; i < _nAttempts; i++) {
            
        std::vector<fastjet::PseudoJet> newAxes(seedAxes.size(),fastjet::PseudoJet(0,0,0,0));

        for (unsigned int i = 0; i < particles.size(); i++) {
            double minDist = 100000000.0; //large number    
            int minJ = -1; //bad ref
            for (unsigned int j = 0; j < seedAxes.size(); j++) {
                double tempDist = _functor->distance(particles[i],seedAxes[j]);
                if (tempDist < minDist) {
                    minDist = tempDist;
                    minJ = j;
                }
            }
               
            if (_functor->do_cluster(particles[i],seedAxes[minJ])) {
                newAxes[minJ] += particles[i];
            }
        }

        seedAxes = newAxes;
            
        double tempTau = _functor->tau(particles,newAxes);

        if (fabs(tempTau - seedTau) < _accuracy) break;
        seedTau = tempTau;
    }
        
    return seedAxes;
}      

} //namespace contrib

FASTJET_END_NAMESPACE
