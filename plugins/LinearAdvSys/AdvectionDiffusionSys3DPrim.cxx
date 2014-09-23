// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdvSys/LinearAdvSys.hh"
#include "AdvectionDiffusionSys3DPrim.hh"
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

Environment::ObjectProvider<AdvectionDiffusionSys3DPrim, DiffusiveVarSet,
	       LinearAdvSysModule, 2>
adSys3DPrimProvider("AdvectionDiffusionSys3DPrim");

//////////////////////////////////////////////////////////////////////////////

AdvectionDiffusionSys3DPrim::AdvectionDiffusionSys3DPrim
(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  AdvectionDiffusionSys3DVarSet(name, model),
  _linearadvModel(model->getConvectiveTerm().d_castTo<LinearAdvSysTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

AdvectionDiffusionSys3DPrim::~AdvectionDiffusionSys3DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void AdvectionDiffusionSys3DPrim::setGradientVars(const vector<RealVector*>& states,
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

void AdvectionDiffusionSys3DPrim::setGradientState(const RealVector& state)
{
}
//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
