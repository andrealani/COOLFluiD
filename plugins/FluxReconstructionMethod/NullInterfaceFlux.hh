// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NullInterfaceFlux_hh
#define COOLFluiD_FluxReconstructionMethod_NullInterfaceFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/BaseInterfaceFlux.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represent the base interface flux computer
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class NullInterfaceFlux : public BaseInterfaceFlux {
public:  // methods

  /// Constructor
  NullInterfaceFlux(const std::string& name);

  /// Destructor
  ~NullInterfaceFlux();
  
  /// Add compute the term to add in the jacobian
  void compute();
  
  /// Set up private data and data
  void setup();
  
  /// Unsetup up private data and data
  void unsetup();
  
  /// Gets the polymorphic type name
  std::string getPolymorphicTypeName() {return getClassName();}
  
}; // class NullInterfaceFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_NullInterfaceFlux_hh

