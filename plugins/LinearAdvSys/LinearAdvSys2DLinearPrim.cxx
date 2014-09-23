// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdvSys/LinearAdvSys.hh"
#include "LinearAdvSys/LinearAdvSys2DLinearPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"
#include "LinearAdvSysPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinearAdvSys2DLinearPrim, JacobianLinearizer, LinearAdvSysModule, 1> linearAdvSys2DLinearPrimProvider("LinearAdvSys2DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

LinearAdvSys2DLinearPrim::LinearAdvSys2DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
 _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<LinearAdvSysTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

LinearAdvSys2DLinearPrim::~LinearAdvSys2DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DLinearPrim::linearize(const vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();

  const CFreal c0x = _model->getc0x();
  const CFreal c1x = _model->getc1x();
  const CFreal c2x = _model->getc2x();
  const CFreal c3x = _model->getc3x();

  const CFreal c0y = _model->getc0y();
  const CFreal c1y = _model->getc1y();
  const CFreal c2y = _model->getc2y();
  const CFreal c3y = _model->getc3y();
    
  linearData[LinearAdvSysTerm::C1X] = c1x;
  linearData[LinearAdvSysTerm::C2X] = c2x;
  linearData[LinearAdvSysTerm::C3X] = c3x;
  linearData[LinearAdvSysTerm::C0Y] = c0y;
  linearData[LinearAdvSysTerm::C1Y] = c1y;
  linearData[LinearAdvSysTerm::C2Y] = c2y; 
  linearData[LinearAdvSysTerm::C3Y] = c3y;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdvSys
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
