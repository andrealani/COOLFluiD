// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_Equidistant_hh
#define COOLFluiD_FluxReconstructionMethod_Equidistant_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/BasePointDistribution.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a Equidistant distribution of points
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
class Equidistant : public BasePointDistribution {
public:  // methods

  /// Constructor
  Equidistant(const std::string& name);

  /// Destructor
  ~Equidistant();
  
  /// Get the 1D coordinates of the point distribution
  std::vector<CFreal> getLocalCoords1D(CFPolyOrder::Type solOrder);
  
  /// Get the maximum distance in 1D between two subsequent points
  CFreal getSubcellResolution(CFPolyOrder::Type solOrder);
  
  /// Set up private data and data
  void setup();
  
  /// Unsetup up private data and data
  void unsetup();
  
  /// Gets the polymorphic type name
  std::string getPolymorphicTypeName() {return getClassName();}
  
}; // class Equidistant

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_Equidistant_hh

