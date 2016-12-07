// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BaseCorrectionFunction_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BaseCorrectionFunction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class represent the base correction function computer
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class BaseCorrectionFunction : public FluxReconstructionSolverStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      FluxReconstructionSolverData,BaseCorrectionFunction > PROVIDER;

public:  // methods

  /// Constructor
  BaseCorrectionFunction(const std::string& name);

  /// Destructor
  virtual ~BaseCorrectionFunction();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BaseCorrectionFunction";
  }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

  /// Set up private data and data
  virtual void setup();
  
  /// Unsetup up private data and data
  virtual void unsetup();
  
private: // data

}; // class BaseCorrectionFunction

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_BaseCorrectionFunction_hh

