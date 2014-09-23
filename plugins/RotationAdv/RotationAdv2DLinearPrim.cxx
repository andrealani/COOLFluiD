// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdv/RotationAdv.hh"
#include "RotationAdv2DLinearPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "RotationAdvPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RotationAdv2DLinearPrim, JacobianLinearizer, RotationAdvModule, 1> rotationAdv2DLinearPrimProvider("RotationAdv2DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

RotationAdv2DLinearPrim::RotationAdv2DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<RotationAdvTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

RotationAdv2DLinearPrim::~RotationAdv2DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DLinearPrim::linearize(const vector<State*>& statesInCell)
{
  cf_assert(_model.isNotNull());

  RealVector& linearData = _model->getPhysicalData();

  const CFreal OX = _model->getOX();
  const CFreal OY = _model->getOY();

  CFreal sumX = 0.;
  CFreal sumY = 0.;

  const CFuint nbStates = statesInCell.size();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    RealVector& node = statesInCell[iState]->getCoordinates();
    sumX += (node[XX]-OX);
    sumY += (node[YY]-OY);
  }

  sumX /= nbStates;
  sumY /= nbStates;

  // this is anti-clockwise rotation
  // which is positive in right-hand convention
  linearData[RotationAdvTerm::VX] = -sumY;
  linearData[RotationAdvTerm::VY] =  sumX;

  if (_model->isClockwise()) linearData *= -1.0; 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
