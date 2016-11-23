// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BaseFluxPntDistribution_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BaseFluxPntDistribution_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represent the base flux point distribution
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class BaseFluxPntDistribution : public FluxReconstructionSolverStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      FluxReconstructionSolverData,BaseFluxPntDistribution > PROVIDER;

public:  // methods

  /// Constructor
  BaseFluxPntDistribution(const std::string& name);

  /// Destructor
  virtual ~BaseFluxPntDistribution();

  /// Add compute the term to add in the jacobian
  virtual void compute() = 0;

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BaseFluxPntDistribution";
  }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

  /// Set up private data and data
  virtual void setup();
  
  /// Unsetup up private data and data
  virtual void unsetup();
  
private: // data

}; // class BaseFluxPntDistribution

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_BaseFluxPntDistribution_hh

