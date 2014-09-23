// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdv/RotationAdv3DLinearPrim.hh"
#include "RotationAdv/RotationAdv.hh"
#include "RotationAdv/RotationAdvPhysicalModel.hh"

#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RotationAdv3DLinearPrim, JacobianLinearizer, RotationAdvModule, 1> rotationAdv3DLinearPrimProvider("RotationAdv3DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

RotationAdv3DLinearPrim::RotationAdv3DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<RotationAdvTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

RotationAdv3DLinearPrim::~RotationAdv3DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv3DLinearPrim::linearize(const vector<State*>& statesInCell)
{
  cf_assert(_model.isNotNull());
  RealVector& linearData = _model->getPhysicalData();

  const CFreal OX = _model->getOX();
   const CFreal OY = _model->getOY();
  const CFreal OZ = _model->getOZ();

  CFreal sumX = 0.;
  CFreal sumY = 0.;
  CFreal sumZ = 0.;
  const CFuint nbStates = statesInCell.size();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
   /// @todo broken after release 2009.3 (instruction "sumY += (node[YY]-OY);" 
   /// was already commented)
     RealVector& node = statesInCell[iState]->getCoordinates();
     sumX += (node[XX]-OX);
     sumY += (node[YY]-OY);
     sumZ += (node[ZZ]-OZ); 
  }

  sumX /= nbStates;
   sumY /= nbStates;
  sumZ /= nbStates;

  // this is anti-clockwise rotation
  // which is positive in right-hand convention
  linearData[RotationAdvTerm::VX] = -sumZ;
  linearData[RotationAdvTerm::VY] =  0.0;
  linearData[RotationAdvTerm::VZ] =  sumX;
  if (_model->isClockwise()) linearData *= -1.0; 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
