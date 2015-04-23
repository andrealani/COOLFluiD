// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Burgers/Burgers.hh"
#include <numeric>

#include "Burgers2DPrim.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Burgers {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Burgers2DPrim, ConvectiveVarSet, BurgersModule, 1>
burgers2DPrimProvider("Burgers2DPrim");

//////////////////////////////////////////////////////////////////////////////

Burgers2DPrim::Burgers2DPrim(Common::SafePtr<BaseTerm> term) :
  Burgers2DVarSet(term),
  _model(term.d_castTo<BurgersTerm>())
{
  vector<std::string> names(1);
  names[0] = "u";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Burgers2DPrim::~Burgers2DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Burgers2DPrim::computeJacobians()
{
  const RealVector& linearData = getModel()->getPhysicalData();

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[0](0,0) = linearData[BurgersTerm::VX];
  (*jacobians)[1](0,0) = linearData[BurgersTerm::VY];
}

//////////////////////////////////////////////////////////////////////////////

void Burgers2DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
                                          RealMatrix& leftEv,
                                          RealVector& eValues,
                                          const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU = linearData[BurgersTerm::VX];
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  rightEv = 1.0;
  leftEv = 1.0;

  eValues[0] = avU*nx + ny;
}

//////////////////////////////////////////////////////////////////////////////

void Burgers2DPrim::computeFlux (const RealVector& vars,
				 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];

  // fx = u^2/2, fy = u
  _fluxArray[0] = 0.5*vars[0]*vars[0]*nx + vars[0]*ny;
}

//////////////////////////////////////////////////////////////////////////////

void Burgers2DPrim::computeFlux (const RealVector& vars)
{
  // fx = u^2/2, fy = u
  _physFlux(0,XX) = 0.5*vars[0]*vars[0];
  _physFlux(0,YY) = vars[0];
}


//////////////////////////////////////////////////////////////////////////////

void Burgers2DPrim::computeStateFlux(const RealVector& vars)
{
 /// @todo broken after release 2009.3
 throw Common::NotImplementedException(FromHere(), "Burgers2DPrim::computeStateFlux()");
}


//////////////////////////////////////////////////////////////////////////////

void Burgers2DPrim::computeScalarJacobian(const RealVector& normal,
                                   RealVector& jacob)
{
  const RealVector& linearData = getModel()->getPhysicalData();
  jacob[0] = linearData[BurgersTerm::VX]*normal[XX] + normal[YY];
}

//////////////////////////////////////////////////////////////////////////////

void Burgers2DPrim::computePhysicalData(const State& state, RealVector& data)
{
  data[BurgersTerm::VX] = state[0];
  data[BurgersTerm::VY] = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Burgers

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
