// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NonLinearAdv/NonLinearAdv.hh"
#include "NonLinearAdv2DPrim.hh"
#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"
#include "NonLinearAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {

//////////////////////////////////////////// //////////////////////////

Environment::ObjectProvider<NonLinearAdv2DPrim, ConvectiveVarSet, NonLinearAdvModule, 1>
NonLinearAdv2DPrimProvider("NonLinearAdv2DPrim");

//////////////////////////////////////////////////////////////////////////////

NonLinearAdv2DPrim::NonLinearAdv2DPrim(Common::SafePtr<BaseTerm> term) :
  NonLinearAdv2DVarSet(term)
{
  vector<std::string> names(1);
  names[0] = "u";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

NonLinearAdv2DPrim::~NonLinearAdv2DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::getAverageAdvectionVector(RealVector& vec, const CFuint iVar) const
{
  const RealVector& linearData = getModel()->getPhysicalData();

  cf_assert(vec.size() == DIM_2D);

  vec[XX] = linearData[NonLinearAdvTerm::VX];
  vec[YY] = linearData[NonLinearAdvTerm::VY];
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::computeJacobians()
{
  const RealVector& linearData = getModel()->getPhysicalData();

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[XX](0,0) = linearData[NonLinearAdvTerm::VX];
  (*jacobians)[YY](0,0) = linearData[NonLinearAdvTerm::VY];
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::setEigenValuesVectors(RealMatrix& rightEv,
            RealMatrix& leftEv,
            RealVector& eValues,
            const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU   = linearData[NonLinearAdvTerm::VX];
  const CFreal avV   = linearData[NonLinearAdvTerm::VY];
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  rightEv = 1.0;
  leftEv = 1.0;

  eValues[0] = avU * nx + avV * ny;
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::computeFlux(
  const State& vars,
	const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];

    _fluxArray[0] = exp(vars[0])*nx + vars[0]*ny;
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::computeFlux(const State& vars)
{
    _physFlux(0,XX) = exp(vars[0]);
    _physFlux(0,YY) = vars[0];
  
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::computeScalarJacobian(
  const RealVector& normal,
  RealVector& jacob)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  jacob[0] = linearData[NonLinearAdvTerm::VX]*normal[XX] +
             linearData[NonLinearAdvTerm::VY]*normal[YY];
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::computePhysicalData(const State& state,
					    RealVector& data)
{

  data[NonLinearAdvTerm::VX] = exp(state[0]);
  data[NonLinearAdvTerm::VY] = 1.0;
}


//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::computeFlux(const RealVector& pdata, const RealVector& normals)
{
 /// @todo broken after release 2009.3 (added method)
 throw Common::NotImplementedException(FromHere(), "NonLinearAdv2DPrim::computeFlux()");
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DPrim::computeStateFlux(const RealVector& vars)
{
 /// @todo broken after release 2009.3 (added method)
 throw Common::NotImplementedException(FromHere(), "NonLinearAdv2DPrim::computeStateFlux()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
