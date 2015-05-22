//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: XConePlugin.hh 748 2014-10-02 06:13:28Z tjwilk $
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

#ifndef __FASTJET_CONTRIB_XConePlugin_HH__
#define __FASTJET_CONTRIB_XConePlugin_HH__

#include <fastjet/config.h>

// #include "Njettiness.hh"
// #include "MeasureDefinition.hh"
// #include "AxesDefinition.hh"
// #include "AxesRefiner.hh"
// #include "TauComponents.hh"
#include "NjettinessPlugin.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib {

//TJW-- UPDATE THIS DOCUMENTATION
/// The XCone jet algorithm
/**
 * An exclusive jet finder that identifies N jets; first N axes are found, then
 * particles are assigned to the nearest (DeltaR) axis and for each axis the
 * corresponding jet is simply the four-momentum sum of these particles.
 *
 * Axes can be found in several ways, specified by the AxesMode argument.  The
 * recommended choices are
 *
 * KT_Axes              : exclusive kT
 * WTA_KT_axes          : exclusive kT with winner-take-all-recombination
 * OnePass_KT_Axes      : one-pass minimization seeded by kt (pretty good)
 * OnePass_WTA_KT_Axes  : one-pass minimization seeded by wta_kt
 *
 * For the UnnormalizedMeasure(beta), N-jettiness is defined as:
 *
 * tau_N = Sum_{all particles i} p_T^i min((DR_i1)^beta, (DR_i2)^beta, ...)
 *
 *   DR_ij is the distance sqrt(Delta_phi^2 + Delta_rap^2) between particle i
 *   and jet j.
 * 
 * The NormalizedMeausure include an extra parameter R0, and the various cutoff
 * measures include an Rcutoff, which effectively defines an angular cutoff
 * similar in effect to a cone-jet radius.
 *
 */

// class XConePlugin : public JetDefinition::Plugin {
class XConePlugin : public NjettinessPlugin {
public:

   // Constructor with N, R, and beta as the options
   XConePlugin(int N, double R0, double beta)
   : NjettinessPlugin(N, OnePass_GenET_GenKT_Axes(calc_delta(beta), calc_power(beta), R0), XConeMeasure(beta, R0))
    , _N(N), _R0(R0), _beta(beta) {}
   
   // The things that are required by base class.
   virtual std::string description () const;
   virtual double R() const {return _R0;} // TODO: make this not stupid
   // virtual void run_clustering(ClusterSequence&) const;

   virtual ~XConePlugin() {}

private:

   double calc_delta(double beta) {
    double delta;
    if (beta > 1) delta = 1/(beta - 1);
    else delta = std::numeric_limits<int>::max(); 
    return delta;
   }

   double calc_power(double beta) {
    return (double)1.0/beta;
   }

   int _N;
   double _R0;
   double _beta;

public:
   
   // // Alternative constructors that define the measure via enums and parameters
   // // These constructors are likely be deprecated in a future version.
   // XConePlugin(int N,
   //                  Njettiness::AxesMode axes_mode,
   //                  Njettiness::MeasureMode measure_mode)
   // : _njettinessFinder(axes_mode, measure_mode, 0), _N(N) {}
   
   
   // XConePlugin(int N,
   //                  Njettiness::AxesMode axes_mode,
   //                  Njettiness::MeasureMode measure_mode,
   //                  double para1)
   // : _njettinessFinder(axes_mode, measure_mode, 1, para1), _N(N) {}
   
   
   // XConePlugin(int N,
   //                  Njettiness::AxesMode axes_mode,
   //                  Njettiness::MeasureMode measure_mode,
   //                  double para1,
   //                  double para2)
   // : _njettinessFinder(axes_mode, measure_mode, 2, para1, para2), _N(N) {}
   
   
   // XConePlugin(int N,
   //                  Njettiness::AxesMode axes_mode,
   //                  Njettiness::MeasureMode measure_mode,
   //                  double para1,
   //                  double para2,
   //                  double para3)
   // : _njettinessFinder(axes_mode, measure_mode, 3, para1, para2, para3), _N(N) {}
   
   
   // // Old constructor for backwards compatibility with v1.0,
   // // where NormalizedCutoffMeasure was the only option
   // XConePlugin(int N,
   //                  Njettiness::AxesMode mode,
   //                  double beta,
   //                  double R0,
   //                  double Rcutoff=std::numeric_limits<double>::max())
   // : _njettinessFinder(mode, NormalizedCutoffMeasure(beta, R0, Rcutoff)), _N(N) {}


};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_XConePlugin_HH__