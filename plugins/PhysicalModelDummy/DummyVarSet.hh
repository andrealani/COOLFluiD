// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_PhysicalModelDummy_DummyVarSet_hh
#define COOLFluiD_PhysicalModelDummy_DummyVarSet_hh

#include "Framework/ConvectiveVarSet.hh"
#include "PhysicalModelDummy/DummyTerm.hh"

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Dummy physical model for conservative variables
class DummyVarSet : public Framework::ConvectiveVarSet {

public: // classes

  /// Constructor
  DummyVarSet(Common::SafePtr< Framework::BaseTerm > term) :
    Framework::ConvectiveVarSet(term),
    _model( term.d_castTo< DummyTerm >() )
  {
  }

  /// Default destructor
  virtual ~DummyVarSet()
  {
  }


protected:

  /// Compute the convective flux
  virtual void computeFlux(
    const Framework::State& vars,
    const RealVector& normals )
  {
    CFLog(INFO,"DummyVarSet::computeFlux( . )\n");
  }

  /// Compute the physical convective flux
  virtual void computeFlux(
    const Framework::State& vars)
  {
    CFLog(INFO,"DummyVarSet::computeFlux( . )\n");
  }

private:

  /// Acquaintance of the term
  Common::SafePtr< DummyTerm > _model;

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

#endif // COOLFluiD_PhysicalModelDummy_DummyVarSet_hh

