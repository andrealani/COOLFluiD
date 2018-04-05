// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NullDiffusiveVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullDiffusiveVarSet,
               DiffusiveVarSet,
               FrameworkLib,
               2>
nullDiffusiveVarSetProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullDiffusiveVarSet::NullDiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
DiffusiveVarSet(name, model)
{
}

//////////////////////////////////////////////////////////////////////////////

NullDiffusiveVarSet::~NullDiffusiveVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NullDiffusiveVarSet::getFlux(const RealVector& values,
                                         const std::vector<RealVector*>& gradients,
                                         const RealVector& normal,
                                         const CFreal& radius)
{
  CFLog(VERBOSE,"NullDiffusiveVarSet::getFlux()" << "\n");
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& NullDiffusiveVarSet::getFlux(const RealVector& values,
                                         const std::vector<RealVector*>& gradients,
                                         const CFreal& radius)
{
  CFLog(VERBOSE,"NullDiffusiveVarSet::getFlux()" << "\n");
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
