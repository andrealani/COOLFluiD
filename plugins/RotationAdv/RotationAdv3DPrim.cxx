// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdv/RotationAdv3DPrim.hh"
#include "RotationAdv/RotationAdv.hh"
#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////// //////////////////////////

Environment::ObjectProvider<RotationAdv3DPrim, ConvectiveVarSet, RotationAdvModule, 1>
rotationAdv3DPrimProvider("RotationAdv3DPrim");

//////////////////////////////////////////////////////////////////////////////

RotationAdv3DPrim::RotationAdv3DPrim(Common::SafePtr<BaseTerm> term) :
  RotationAdv3DVarSet(term)
{
  vector<std::string> names(1);
  names[0] = "u";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

RotationAdv3DPrim::~RotationAdv3DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv3DPrim::getAverageAdvectionVector(RealVector& vec, const CFuint iVar) const
{
  const RealVector& linearData = getModel()->getPhysicalData();

  cf_assert(vec.size() == DIM_3D);

  vec[XX] = linearData[RotationAdvTerm::VX];
  vec[YY] = linearData[RotationAdvTerm::VY];
  vec[ZZ] = linearData[RotationAdvTerm::VZ];
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv3DPrim::computeJacobians()
{
  const RealVector& linearData =
    getModel()->getPhysicalData();

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[XX](0,0) = linearData[RotationAdvTerm::VX];
  (*jacobians)[YY](0,0) = linearData[RotationAdvTerm::VY];
  (*jacobians)[ZZ](0,0) = linearData[RotationAdvTerm::VZ];
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv3DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
            RealMatrix& leftEv,
            RealVector& eValues,
            const RealVector& normal)
{
  const RealVector& linearData =
    getModel()->getPhysicalData();

  const CFreal avU   = linearData[RotationAdvTerm::VX];
  const CFreal avV   = linearData[RotationAdvTerm::VY];
  const CFreal avW   = linearData[RotationAdvTerm::VZ];
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];

  rightEv = 1.0;
  leftEv = 1.0;

  eValues[0] = avU * nx + avV * ny +  avW * nz;
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv3DPrim::computeScalarJacobian(const RealVector& normal,
               RealVector& jacob)
{
  const RealVector& linearData =
    getModel()->getPhysicalData();

  jacob[0] = linearData[RotationAdvTerm::VX]*normal[XX] +
             linearData[RotationAdvTerm::VY]*normal[YY] +
             linearData[RotationAdvTerm::VZ]*normal[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv3DPrim::computePhysicalData(const State& state,
					    RealVector& data)
{
  RealVector& node = state.getCoordinates();
  const CFreal OX = getModel()->getOX();
  const CFreal OZ = getModel()->getOZ();

  data[RotationAdvTerm::VX] = - (node[ZZ] - OZ);
  data[RotationAdvTerm::VY] =    0.0;
  data[RotationAdvTerm::VZ] =   (node[XX] - OX);

  if (getModel()->isClockwise()) data *= -1.0;
}

/////////////////////////////////////////////////////////////////////////////

void RotationAdv3DPrim::computeStateFromPhysicalData(const RealVector& data,
                                                     State& state)
{
  state[0] = data[0];
}
                                               
//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
