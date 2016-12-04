// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_BasePointDistribution_hh
#define COOLFluiD_FluxReconstructionMethod_BasePointDistribution_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represent the base 1D point distribution for either flux or solution points
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class BasePointDistribution : public FluxReconstructionSolverStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      FluxReconstructionSolverData,BasePointDistribution > PROVIDER;

public:  // methods

  /// Constructor
  BasePointDistribution(const std::string& name);

  /// Destructor
  virtual ~BasePointDistribution();

  /// Get the 1D coordinates of the point distribution 
  virtual std::vector<CFreal> getLocalCoords1D(CFPolyOrder::Type solOrder) = 0;

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BasePointDistribution";
  }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

  /// Set up private data and data
  virtual void setup();
  
  /// Unsetup up private data and data
  virtual void unsetup();
  
private: // data
  

}; // class BasePointDistribution

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_BasePointDistribution_hh

