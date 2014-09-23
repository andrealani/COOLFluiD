// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdv/LinearAdv.hh"
#include "AdvectionDiffusion3DPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "LinearAdv/LinearAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AdvectionDiffusion3DPrim, DiffusiveVarSet,
	       LinearAdvModule, 2>
ad3DPrimProvider("AdvectionDiffusion3DPrim");

//////////////////////////////////////////////////////////////////////////////

AdvectionDiffusion3DPrim::AdvectionDiffusion3DPrim
(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  AdvectionDiffusion3DVarSet(name, model),
  _linearadvModel(model->getConvectiveTerm().d_castTo<LinearAdvTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

AdvectionDiffusion3DPrim::~AdvectionDiffusion3DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void AdvectionDiffusion3DPrim::setGradientVars(const vector<RealVector*>& states,
					 const vector<RealVector*>& values,
					 const CFuint stateSize)
{
  const CFuint nbStates = stateSize;

   RealVector& u = *values[0];
  for (CFuint i = 0; i < nbStates; ++i) {
    const RealVector& state = *states[i];
    u[i] = state[0];
  }
}

//////////////////////////////////////////////////////////////////////////////

void AdvectionDiffusion3DPrim::setGradientState(const RealVector& state)
{
}
//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
