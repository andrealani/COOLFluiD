// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_NullCorrectionFunction_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_NullCorrectionFunction_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent the null correction function computer
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
class NullCorrectionFunction : public BaseCorrectionFunction {
public:  // methods

  /// Constructor
  NullCorrectionFunction(const std::string& name);

  /// Destructor
  ~NullCorrectionFunction();

  /// Compute the Null correction function of an instance of FluxReconstructionElementData
  void computeCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< RealVector > >& corrcts);
    
  /// Compute the divergence of the Null correction function of an instance of FluxReconstructionElementData
  void computeDivCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< CFreal > >& corrcts);
  
  /// Set up private data and data
  void setup();
  
  /// Unsetup up private data and data
  void unsetup();
  
  /// Gets the polymorphic type name
  std::string getPolymorphicTypeName() {return getClassName();}
  
}; // class NullCorrectionFunction

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_NullCorrectionFunction_hh

