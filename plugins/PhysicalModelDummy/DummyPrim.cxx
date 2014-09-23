// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "PhysicalModelDummy/Dummy.hh"
#include "PhysicalModelDummy/DummyPrim.hh"
#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<
    DummyPrim,
    ConvectiveVarSet,
    DummyModule,
    1 >
  dummyPrimProvider("PhysicalModelDummyPrim");

//////////////////////////////////////////////////////////////////////////////

DummyPrim::DummyPrim(Common::SafePtr< BaseTerm > term) :
  DummyVarSet(term)
{
  setVarNames( term.d_castTo< DummyTerm >()->getVarNames() );
}

//////////////////////////////////////////////////////////////////////////////

DummyPrim::~DummyPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void DummyPrim::computeFlux(const State& vars, const RealVector& normals)
{
  CFLog(INFO,"DummyPrim::computeFlux( . )\n");
}

//////////////////////////////////////////////////////////////////////////////

void DummyPrim::computeFlux(const State& vars)
{
  CFLog(INFO,"DummyPrim::computeFlux( )\n");
}

//////////////////////////////////////////////////////////////////////////////

void DummyPrim::computeFlux (const RealVector& pdata, const RealVector& normals)
{
  CFLog(INFO,"DummyPrim::computeFlux( )\n");
}

//////////////////////////////////////////////////////////////////////////////

void DummyPrim::computeStateFlux (const RealVector& pdata)
{
  CFLog(INFO,"DummyPrim::computeStateFlux( )\n");
}

//////////////////////////////////////////////////////////////////////////////

void DummyPrim::computePhysicalData (const State& state, RealVector& pdata)
{
  CFLog(INFO,"DummyPrim::computePhysicalData( )\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

