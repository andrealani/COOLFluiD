// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdv/RotationAdv.hh"
#include "RotationAdv2DPrim.hh"
#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"
#include "RotationAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////// //////////////////////////

Environment::ObjectProvider<RotationAdv2DPrim, ConvectiveVarSet, RotationAdvModule, 1>
rotationAdv2DPrimProvider("RotationAdv2DPrim");

//////////////////////////////////////////////////////////////////////////////

RotationAdv2DPrim::RotationAdv2DPrim(Common::SafePtr<BaseTerm> term) :
  RotationAdv2DVarSet(term)
{
  vector<std::string> names(1);
  names[0] = "u";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

RotationAdv2DPrim::~RotationAdv2DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DPrim::getAverageAdvectionVector(RealVector& vec, const CFuint iVar) const
{
  const RealVector& linearData = getModel()->getPhysicalData();

  cf_assert(vec.size() == DIM_2D);

  vec[XX] = linearData[RotationAdvTerm::VX];
  vec[YY] = linearData[RotationAdvTerm::VY];
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DPrim::computeJacobians()
{
  const RealVector& linearData = getModel()->getPhysicalData();

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[XX](0,0) = linearData[RotationAdvTerm::VX];
  (*jacobians)[YY](0,0) = linearData[RotationAdvTerm::VY];
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
            RealMatrix& leftEv,
            RealVector& eValues,
            const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU   = linearData[RotationAdvTerm::VX];
  const CFreal avV   = linearData[RotationAdvTerm::VY];
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  rightEv = 1.0;
  leftEv = 1.0;

  eValues[0] = avU * nx + avV * ny;
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DPrim::computeScalarJacobian(
  const RealVector& normal,
  RealVector& jacob)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  jacob[0] = linearData[RotationAdvTerm::VX]*normal[XX] +
             linearData[RotationAdvTerm::VY]*normal[YY];
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DPrim::computePhysicalData(const State& state,
					    RealVector& data)
{
  RealVector& node = state.getCoordinates();

  const CFreal OX = getModel()->getOX();
  const CFreal OY = getModel()->getOY();

  data[RotationAdvTerm::u] = state[0]; 
  data[RotationAdvTerm::VX] = - (node[YY] - OY);
  data[RotationAdvTerm::VY] =   (node[XX] - OX);

  if (getModel()->isClockwise()) data *= -1.0;
}

//////////////////////////////////////////////////////////////////////////////
//
void RotationAdv2DPrim::computeStateFromPhysicalData(const RealVector& data,
                                              State& state)
{
  state[0] = data[0];
}
                                                


//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
