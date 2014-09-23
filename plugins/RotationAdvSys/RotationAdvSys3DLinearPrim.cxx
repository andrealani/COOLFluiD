// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdvSys/RotationAdvSys3DLinearPrim.hh"
#include "RotationAdvSys/RotationAdvSys.hh"
#include "RotationAdvSys/RotationAdvSysPhysicalModel.hh"

#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RotationAdvSys3DLinearPrim, JacobianLinearizer, RotationAdvSysModule, 1> rotationAdvsys3DLinearPrimProvider("RotationAdvSys3DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

RotationAdvSys3DLinearPrim::RotationAdvSys3DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().d_castTo<RotationAdvSysTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

RotationAdvSys3DLinearPrim::~RotationAdvSys3DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys3DLinearPrim::linearize(const vector<State*>& statesInCell)
{
  cf_assert(_model.isNotNull());

  RealVector& linearData = _model->getPhysicalData();


cf_assert(_model.isNotNull());


  const CFreal OX0 = _model->getOX0();
  const CFreal OY0 = _model->getOY0();
  const CFreal OZ0 = _model->getOZ0();


  const CFreal OX1 = _model->getOX1();
  const CFreal OY1 = _model->getOY1();
  const CFreal OZ1 = _model->getOZ1();


  const CFreal OX2 = _model->getOX2();
  const CFreal OY2 = _model->getOY2();
  const CFreal OZ2 = _model->getOZ2();


  const CFreal OX3 = _model->getOX3();
  const CFreal OY3 = _model->getOY3();
  const CFreal OZ3 = _model->getOZ3();



  CFreal sumX0 = 0.;
  CFreal sumY0 = 0.;
  CFreal sumZ0 = 0.;


  CFreal sumX1 = 0.;
  CFreal sumY1 = 0.;
  CFreal sumZ1 = 0.;


  CFreal sumX2 = 0.;
  CFreal sumY2 = 0.;
  CFreal sumZ2 = 0.;


  CFreal sumX3 = 0.;
  CFreal sumY3 = 0.;
  CFreal sumZ3 = 0.;


  const CFuint nbStates = statesInCell.size();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    RealVector& node = statesInCell[iState]->getCoordinates();
    sumX0 += (node[XX]-OX0);
    sumY0 += (node[YY]-OY0);
    sumZ0 += (node[ZZ]-OZ0);

    sumX1 += (node[XX]-OX1);
    sumY1 += (node[YY]-OY1);
    sumZ1 += (node[ZZ]-OZ1);

    sumX2 += (node[XX]-OX2);
    sumY2 += (node[YY]-OY2);
    sumZ2 += (node[ZZ]-OZ2);


    sumY3 += (node[YY]-OY3);
    sumX3 += (node[XX]-OX3);
    sumZ3 += (node[ZZ]-OZ3);


  }

  sumX0 /= nbStates;
  sumY0 /= nbStates;
  sumZ0 /= nbStates;


  sumX1 /= nbStates;
  sumY1 /= nbStates;
  sumZ1 /= nbStates;


  sumX2 /= nbStates;
  sumY2 /= nbStates;
  sumZ2 /= nbStates;


  sumX3 /= nbStates;
  sumY3 /= nbStates;
  sumZ3 /= nbStates;


  // this is anti-clockwise rotation
  // which is positive in right-hand convention
     linearData[RotationAdvSysTerm::C0X] = -sumZ0;
     linearData[RotationAdvSysTerm::C0Y] =  0.0;
     linearData[RotationAdvSysTerm::C0Z] =  sumX0;
 
 
     linearData[RotationAdvSysTerm::C1X] = -sumZ1;
     linearData[RotationAdvSysTerm::C1Y] =  0.0;
     linearData[RotationAdvSysTerm::C1Z] =  sumX1;

     linearData[RotationAdvSysTerm::C2X] = -sumZ2;
     linearData[RotationAdvSysTerm::C2Y] =  0.0;
     linearData[RotationAdvSysTerm::C2Z] =  sumX2;

  
     linearData[RotationAdvSysTerm::C3X] = -sumZ3;
     linearData[RotationAdvSysTerm::C3Y] =  0.0;
     linearData[RotationAdvSysTerm::C3Z] =  sumX3;

  
  
    if (_model->isClockwise()) linearData *= -1.0;
  


}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
