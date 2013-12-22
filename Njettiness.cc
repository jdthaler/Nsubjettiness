//NEW FILE CREATED BY TJW
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


#include "Nsubjettiness.hh"
#include "Njettiness.hh"
#include "NjettinessPlugin.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Measure Functor
//
///////

// N-subjettiness pieces in MeasureFunctor moved from Njettiness.hh -- TJW

std::vector<double> MeasureFunctor::subtaus_numerator(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& axes) {
// Returns the unnormalized sub-tau values, i.e. a std::vector of the contributions to tau_N of each Voronoi region (or region within R_0)
      
   std::vector<double> tauNum(axes.size(), 0.0);
   
   for (unsigned i = 0; i < particles.size(); i++) {
      // find minimum distance; start with 0'th axis for reference
      int j_min = 0;
      double minR = distance(particles[i],axes[0]);
      for (unsigned j = 1; j < axes.size(); j++) {
         double tempR = distance(particles[i],axes[j]); // delta R distance
         if (tempR < minR) {minR = tempR; j_min = j;}
      }
      tauNum[j_min] += numerator(particles[i],axes[j_min]);
   }

   return tauNum;
}

// Calculates unnormalized tau_N
double MeasureFunctor::tau_numerator(const std::vector <fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) {
   std::vector<double> tau_vec = subtaus_numerator(particles, axes);
   double tauNum = 0.0;
   for (unsigned i = 0; i < tau_vec.size(); i++) {
      tauNum += tau_vec[i];
   }
   return tauNum;
}

//Calculates normalization for tau and sub Tau 
double MeasureFunctor::tau_denominator(const std::vector <fastjet::PseudoJet>& particles) {
   double tauDen = 0.0;
   for (unsigned i = 0; i < particles.size(); i++) {
      tauDen += denominator(particles[i]);
   }
   return tauDen;
}

//calculates normalized subTaus
std::vector<double> MeasureFunctor::subtaus_normalized(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet> &axes){
   std::vector<double> tauNum = subtaus_numerator(particles, axes);
   std::vector<double> tau_normalized(axes.size(), 0.0);
   double tauDen = tau_denominator(particles);
   for (unsigned i = 0; i < axes.size(); i++) {
      tau_normalized[i] = tauNum[i]/tauDen;
   }
   return tau_normalized;
}

//calculates normalized tau_N
double MeasureFunctor::tau_normalized(const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet> &axes){
   return tau_numerator(particles, axes)/tau_denominator(particles);
}

///////
//
// Axes Finder Options
//
///////

//functions for WTA Recombiner moved from Njettiness.hh -- TJW

void WTARecombiner::recombine(const fastjet::PseudoJet &pa, const fastjet::PseudoJet &pb, fastjet::PseudoJet &pab) const{
  if (pa.perp() > pb.perp()) {
    pab.reset_PtYPhiM(pa.perp() + pb.perp(), pa.rap(), pa.phi());
  }
  else if (pb.perp() > pa.perp()) {
    pab.reset_PtYPhiM(pa.perp() + pb.perp(), pb.rap(), pb.phi());
  }
}

std::string WTARecombiner::description() const{
  return "WTA scheme recombination";
}

///////
//
// Functions for minimization.
//
///////

// Minimization function definitions moved from Njettiness.hh -- TJW

// Given starting axes, update to find better axes
template <int N>
std::vector<LightLikeAxis> AxesFinderFromKmeansMinimization::UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes, 
                                  const std::vector <fastjet::PseudoJet> & inputJets,
                                  NsubParameters paraNsub, double precision) {
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

   
   double beta = paraNsub.beta();
   double Rcutoff = paraNsub.Rcutoff();
   
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

// Given starting axes, update to find better axes
// (This is just a wrapper for the templated version above.)
std::vector<LightLikeAxis> AxesFinderFromKmeansMinimization::UpdateAxes(const std::vector <LightLikeAxis> & old_axes, 
                                      const std::vector <fastjet::PseudoJet> & inputJets, NsubParameters paraNsub, double precision) {
   int N = old_axes.size();
   switch (N) {
      case 1: return UpdateAxesFast<1>(old_axes, inputJets, paraNsub, precision);
      case 2: return UpdateAxesFast<2>(old_axes, inputJets, paraNsub, precision);
      case 3: return UpdateAxesFast<3>(old_axes, inputJets, paraNsub, precision);
      case 4: return UpdateAxesFast<4>(old_axes, inputJets, paraNsub, precision);
      case 5: return UpdateAxesFast<5>(old_axes, inputJets, paraNsub, precision);
      case 6: return UpdateAxesFast<6>(old_axes, inputJets, paraNsub, precision);
      case 7: return UpdateAxesFast<7>(old_axes, inputJets, paraNsub, precision);
      case 8: return UpdateAxesFast<8>(old_axes, inputJets, paraNsub, precision);
      case 9: return UpdateAxesFast<9>(old_axes, inputJets, paraNsub, precision);
      case 10: return UpdateAxesFast<10>(old_axes, inputJets, paraNsub, precision);
      case 11: return UpdateAxesFast<11>(old_axes, inputJets, paraNsub, precision);
      case 12: return UpdateAxesFast<12>(old_axes, inputJets, paraNsub, precision);
      case 13: return UpdateAxesFast<13>(old_axes, inputJets, paraNsub, precision);
      case 14: return UpdateAxesFast<14>(old_axes, inputJets, paraNsub, precision);
      case 15: return UpdateAxesFast<15>(old_axes, inputJets, paraNsub, precision);
      case 16: return UpdateAxesFast<16>(old_axes, inputJets, paraNsub, precision);
      case 17: return UpdateAxesFast<17>(old_axes, inputJets, paraNsub, precision);
      case 18: return UpdateAxesFast<18>(old_axes, inputJets, paraNsub, precision);
      case 19: return UpdateAxesFast<19>(old_axes, inputJets, paraNsub, precision);
      case 20: return UpdateAxesFast<20>(old_axes, inputJets, paraNsub, precision);
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

// Get minimization axes
std::vector<fastjet::PseudoJet> AxesFinderFromKmeansMinimization::GetMinimumAxes(const std::vector <fastjet::PseudoJet> & seedAxes, const std::vector <fastjet::PseudoJet> & inputJets, KmeansParameters para, 
                                          NsubParameters paraNsub,MeasureFunctor* functor) {
   int n_jets = seedAxes.size();
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
         new_axes = UpdateAxes(old_axes, inputJets, paraNsub, para.precision()); // Update axes
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

///////
//
// Main Njettiness Class
//
///////

//all following Njettiness functions moved from Njettiness.hh -- TJW

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
