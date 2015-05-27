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

#ifndef __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__
#define __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__

#include <fastjet/config.h>

#include "Njettiness.hh"
#include "MeasureDefinition.hh"
#include "AxesDefinition.hh"
#include "TauComponents.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib {


/// The Njettiness jet algorithm
/**
 * An exclusive jet finder that identifies N jets; first N axes are found, then
 * particles are assigned to the nearest (DeltaR) axis and for each axis the
 * corresponding jet is simply the four-momentum sum of these particles.
 *
 * As of version 2.2, it is recommended to use the XConePlugin, which has
 * sensible default values for jet finding.
 *
 * Axes can be found in several ways, specified by the AxesDefinition argument.
 * For recommendations on which axes to use, please see the README file.
 * 
 * Jet regions are determined by the MeasureDefinition. For example, 
 * for the UnnormalizedMeasure(beta), N-jettiness is defined as:
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
 * Other measures introduced in version 2.2 include OriginalGeometricMeasure,
 * ModifiedGeometricMeasure, and ConicalGeometricMeasure, which define N-jettiness
 * through dot products of particle momenta with light-like axes. OriginalGeometricMeasure
 * produces football-shaped jets due to its central weighting of the beam measure,
 * but ModifiedGeometric and ConicalGeometric both deform the original geometric measure
 * to allow for cone-shaped jets. The size of these cones can be controlled through Rcutoff
 * just as in the other measures. See the README file or MeasureDefinition.hh for information
 * on how to call these measures.

 */

class NjettinessPlugin : public JetDefinition::Plugin {
public:

   // Constructor with same arguments as Nsubjettiness.
   NjettinessPlugin(int N,
                    const AxesDefinition & axes_def,
                    const MeasureDefinition & measure_def)
   : _njettinessFinder(axes_def, measure_def), _N(N) {}
   

   // The things that are required by base class.
   virtual std::string description () const;
   virtual double R() const {return -1.0;} // TODO: make this not stupid
   virtual void run_clustering(ClusterSequence&) const;

   // Using manual axes with Njettiness Plugin
   void setAxes(const std::vector<fastjet::PseudoJet> & myAxes) {
      // Cross check that manual axes are being used is in Njettiness
      _njettinessFinder.setAxes(myAxes);
   }

   virtual ~NjettinessPlugin() {}

private:

   Njettiness _njettinessFinder;
   int _N;

public:
   
   // Alternative constructors that define the measure via enums and parameters
   // These constructors are likely be deprecated in a future version.
   NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode)
   : _njettinessFinder(axes_mode, measure_mode, 0), _N(N) {}
   
   
   NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode,
                    double para1)
   : _njettinessFinder(axes_mode, measure_mode, 1, para1), _N(N) {}
   
   
   NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode,
                    double para1,
                    double para2)
   : _njettinessFinder(axes_mode, measure_mode, 2, para1, para2), _N(N) {}
   
   
   NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode,
                    double para1,
                    double para2,
                    double para3)
   : _njettinessFinder(axes_mode, measure_mode, 3, para1, para2, para3), _N(N) {}
   
   
   // Old constructor for backwards compatibility with v1.0,
   // where NormalizedCutoffMeasure was the only option
   NjettinessPlugin(int N,
                    Njettiness::AxesMode mode,
                    double beta,
                    double R0,
                    double Rcutoff=std::numeric_limits<double>::max())
   : _njettinessFinder(mode, NormalizedCutoffMeasure(beta, R0, Rcutoff)), _N(N) {}


};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__