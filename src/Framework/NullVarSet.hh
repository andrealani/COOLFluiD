// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullVarSet_hh
#define COOLFluiD_Framework_NullVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConvectiveVarSet.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a null variables transformer
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API NullVarSet : public ConvectiveVarSet {
public:

  /// Default constructor without arguments
  NullVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /// Default destructor
  ~NullVarSet();

  /// Set up the private data and give the maximum size of states physical
  /// data to store
  void setup()
  {
    CFLog(VERBOSE, "Calling NullVarSet::setup() => " <<
          "this is a Null Variable Set" << "\n");
  }

  /// Gets the block separator for this variable set
  CFuint getBlockSeparator() const
  {
    CFLog(VERBOSE, "Calling NullVarSet::getBlockSeparator() => " <<
          "this is a Null Variable Set" << "\n");
    return 0;
  }

  /// Set the jacobians
  virtual void computeJacobians();

  /// Set the matrix of the right eigenvectors and the matrix of the eigenvalues
  void computeEigenValuesVectors(RealMatrix& rightEv,
				 RealMatrix& leftEv,
				 RealVector& eValues,
				 const RealVector& normal);
  
  /// Set the PhysicalData corresponding to the given State
  void computePhysicalData(const State& state, RealVector& data)
  {
    CFLog(VERBOSE,"NullVarSet::computePhysicalData() => "
	  << "this is a Null VarSet" << "\n");
  }
  
  /// Set the PhysicalData corresponding to the given State
  void computeStateFromPhysicalData(const RealVector& data, State& state)
  {
    CFLog(VERBOSE,"NullVarSet::computeStateFromPhysicalData() => "
  << "this is a Null VarSet" << "\n");
  }

  /// Compute the convective flux
  void computeFlux(const RealVector& pdata, const RealVector& normals);
  
  /// Compute the physical convective flux
  void computeStateFlux(const RealVector& pdata);
 
  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar);
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
  }
  
}; // end of class NullVarSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullVarSet_hh
