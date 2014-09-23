// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdvSys/LinearAdvSys.hh"
#include "AdvectionDiffusionSys2DPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "LinearAdvSys/LinearAdvSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


  namespace Physics {
namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AdvectionDiffusionSys2DPrim, DiffusiveVarSet,
	       LinearAdvSysModule, 2>
adSys2DPrimProvider("AdvectionDiffusionSys2DPrim");

//////////////////////////////////////////////////////////////////////////////

AdvectionDiffusionSys2DPrim::AdvectionDiffusionSys2DPrim
(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  AdvectionDiffusionSys2DVarSet(name, model),
  _linearadvModel(model->getConvectiveTerm().d_castTo<LinearAdvSysTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

AdvectionDiffusionSys2DPrim::~AdvectionDiffusionSys2DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void AdvectionDiffusionSys2DPrim::setGradientVars(const vector<RealVector*>& states,
					 RealMatrix& values,
					 const CFuint stateSize)
{
  const CFuint nbStates = stateSize;

  for (CFuint i = 0; i < nbStates; ++i) {
    const RealVector& state = *states[i];
    values(0,i) = state[0];
    values(1,i) = state[1];
    values(2,i) = state[2];
    values(3,i) = state[3];
  }
}

//////////////////////////////////////////////////////////////////////////////

void AdvectionDiffusionSys2DPrim::setGradientState(const RealVector& state)
{
}
//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
