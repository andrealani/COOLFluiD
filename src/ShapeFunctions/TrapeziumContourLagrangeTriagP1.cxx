// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "TrapeziumContourLagrangeTriagP1.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

IntegratorImplProvider<ContourIntegratorImpl,
                       TrapeziumContourLagrangeTriagP1>
trapeziumContourLagrangeTriagP1Provider;

//////////////////////////////////////////////////////////////////////////////

void TrapeziumContourLagrangeTriagP1::setup()
{
  _coeff.resize(3);
  for (CFuint iFace = 0; iFace < _coeff.size(); ++iFace) {
    _coeff[iFace].resize(2);
    _coeff[iFace][0] = 0.5;
    _coeff[iFace][1] = 0.5;
  }

  _pattern.setNbShapeFunctions(3);
}

//////////////////////////////////////////////////////////////////////////////

void TrapeziumContourLagrangeTriagP1::computeSolutionAtQuadraturePoints(
  const std::vector<Framework::State*>& states,
        std::vector<Framework::State*>& values)
{
  CFLogDebugMin("TrapeziumContourLagrangeTriagP1::computeSolutionAtQuadraturePoints()\n");
  
  copy(*states[0], *values[0]);
  copy(*states[1], *values[1]);
  copy(*states[1], *values[2]);
  copy(*states[2], *values[3]);
  copy(*states[2], *values[4]);
  copy(*states[0], *values[5]);
}

//////////////////////////////////////////////////////////////////////////////

void TrapeziumContourLagrangeTriagP1::computeFaceJacobianDetAtQuadraturePoints
(const std::vector<Framework::Node*>& nodes,
 std::vector<RealVector>& faceJacobian)
{
  CFLogDebugMin("TrapeziumContourLagrangeTriagP1::computeFaceJacobianDetAtQuadraturePoints()\n");
  
  cf_assert(faceJacobian.size() >= 3);
  const CFreal face01 = MathFunctions::getDistance(*nodes[0], *nodes[1]);
  const CFreal face12 = MathFunctions::getDistance(*nodes[1], *nodes[2]);
  const CFreal face20 = MathFunctions::getDistance(*nodes[2], *nodes[0]);
  
  faceJacobian[0][0] = faceJacobian[0][1] = face01;
  faceJacobian[1][0] = faceJacobian[1][1] = face12;
  faceJacobian[2][0] = faceJacobian[2][1] = face20;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
