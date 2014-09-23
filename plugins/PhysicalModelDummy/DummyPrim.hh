// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_PhysicalModelDummy_DummyPrim_hh
#define COOLFluiD_PhysicalModelDummy_DummyPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "PhysicalModelDummy/DummyVarSet.hh"

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Dummy physical model for conservative variables
class DummyPrim : public DummyVarSet {

public:

  /// Constructor
  DummyPrim(Common::SafePtr< Framework::BaseTerm > term);

  /// Default destructor
  ~DummyPrim();

  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
  }
  
protected:

  /// Compute the convective flux
  virtual void computeFlux(const Framework::State& vars,const RealVector& normals );

  /// Compute the physical convective flux
  virtual void computeFlux(const Framework::State& vars);

  /// Compute the convective flux, projected in normal
  virtual void computeFlux (const RealVector& pdata, const RealVector& normals);

  /// Compute the physical convective flux
  virtual void computeStateFlux (const RealVector& pdata);

  /// Set the PhysicalData corresponding to the given State
  virtual void computePhysicalData (const Framework::State& state, RealVector& pdata);

}; // end of class DummyPrim

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

#endif // COOLFluiD_PhysicalModelDummy_DummyPrim_hh

