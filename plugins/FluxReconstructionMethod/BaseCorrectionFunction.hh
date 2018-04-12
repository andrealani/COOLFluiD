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
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/** 
 * This class represent the base correction function computer
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
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
    
  /**
   * Compute the correction function of an instance of FluxReconstructionElementData. 
   * corrfcts contains for each solution point the contrbution of each flux point (being a realvector with the dimensionality of the test case)
   */
  virtual void computeCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< RealVector > >& corrcts) = 0;
    
  /**
   * Compute the divergence of the correction function of an instance of FluxReconstructionElementData. 
   * corrfcts contains for each solution point the contrbution of each flux point
   */
  virtual void computeDivCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< CFreal > >& corrcts) = 0;

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

