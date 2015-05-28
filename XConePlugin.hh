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

#ifndef __FASTJET_CONTRIB_XCONEPLUGIN_HH__
#define __FASTJET_CONTRIB_XCONEPLUGIN_HH__

#include <fastjet/config.h>

#include "NjettinessPlugin.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib {

/// The XCone jet algorithm
/**
 * An exclusive jet finder that identifies N jets; first N axes are found, then
 * particles are assigned to the nearest (DeltaR) axis and for each axis the
 * corresponding jet is simply the four-momentum sum of these particles.
 *
 * The XConePlugin is based on NjettinessPlugin, but with sensible default
 * values for the AxesDefinition and MeasureDefinition.  There are three arguments
 *
 *   int N:       number of exclusive jets to be found
 *   double R0:   approximate jet radius
 *   double beta: determines style of jet finding with the recommended values being:
 *                beta = 2:  standard "mean" jets where jet momentum/axis align approximately.
 *                beta = 1:  recoil-free "median" variant where jet axis points at hardest cluster.
 *
 * The AxesDefinition is OnePass_GenET_GenKT_Axes, which uses a generalized kT
 * clustering algorithm matched to the beta value.
 *
 * The MeasureDefinition is the XConeMeasure, which is based on the
 * ConicalGeometric measure.
 *
 */

// class XConePlugin
class XConePlugin : public NjettinessPlugin {
public:

   // Constructor with N, R0, and beta as the options.  beta = 2.0 is the default
   // All this does is use the NjettinessPlugin with OnePass_GenET_GenKT_Axes and the XConeMeasure.
   // For more advanced usage, call NjettinessPlugin directly
   // Note that the order of the R0 and beta values is reversed from the XConeMeasure to
   // standard usage for Plugins.
   XConePlugin(int N, double R0, double beta = 2.0)
   : NjettinessPlugin(N,
                      OnePass_GenET_GenKT_Axes(calc_delta(beta), calc_power(beta), R0), // use recommended axes method only
                      XConeMeasure(beta, R0)  // use recommended XCone measure.
                      ),
   _N(N), _R0(R0), _beta(beta)
   {}
   
   // The things that are required by base class.
   virtual std::string description () const;
   virtual double R() const {return _R0;}
   
   // run_clustering is done by NjettinessPlugin

   virtual ~XConePlugin() {}

private:

   // static calls to use with in the constructor, sets the recommended delta value
   static double calc_delta(double beta) {
      double delta;
      if (beta > 1) delta = 1/(beta - 1);
      else delta = std::numeric_limits<int>::max(); // use winner take all
      return delta;
   }

   // static calls to use with in the constructor, sets the recommended p value
   static double calc_power(double beta) {
      return (double) 1.0/beta;
   }

   double _N;
   double _R0;
   double _beta;

public:

};

   
// class PseudoXConePlugin
// A "poor man's" version of XCone with no minimization step
// Right now, just for testing purposes by the developers
class PseudoXConePlugin : public NjettinessPlugin {
public:
   
   // Constructor with N, R0, and beta as the options.  beta = 2.0 is the default
   // All this does is use the NjettinessPlugin with GenET_GenKT_Axes and the XConeMeasure.
   PseudoXConePlugin(int N, double R0, double beta = 2.0)
   : NjettinessPlugin(N,
                      GenET_GenKT_Axes(calc_delta(beta), calc_power(beta), R0), // poor man's axes
                      XConeMeasure(beta, R0)  // use recommended XCone measure.
                      ),
   _N(N), _R0(R0), _beta(beta)
   {}
   
   // The things that are required by base class.
   virtual std::string description () const;
   virtual double R() const {return _R0;}
   
   // run_clustering is done by NjettinessPlugin
   
   virtual ~PseudoXConePlugin() {}
   
private:
   
   // static calls to use with in the constructor, sets the recommended delta value
   static double calc_delta(double beta) {
      double delta;
      if (beta > 1) delta = 1/(beta - 1);
      else delta = std::numeric_limits<int>::max(); // use winner take all
      return delta;
   }
   
   // static calls to use with in the constructor, sets the recommended p value
   static double calc_power(double beta) {
      return (double) 1.0/beta;
   }
   
   double _N;
   double _R0;
   double _beta;
   
public:
   
};

   
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_XConePlugin_HH__