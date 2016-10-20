// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_FluxReconstructionStrategy_hh
#define COOLFluiD_FluxReconstructionMethod_FluxReconstructionStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a flux reconstruction strategy
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class FluxReconstructionStrategy : public FluxReconstructionSolverStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      FluxReconstructionSolverData,FluxReconstructionStrategy > PROVIDER;

public:  // methods

  /// Constructor
  FluxReconstructionStrategy(const std::string& name);

  /// Destructor
  ~FluxReconstructionStrategy();

  /// Add compute the term to add in the jacobian
  void compute();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "FluxReconstructionStrategy";
  }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

  /// Set up private data and data
  virtual void setup();

  /// Set up private data and data
  virtual void unsetup();

private: // data

}; // class FluxReconstructionStrategy

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_FluxReconstructionStrategy_hh

