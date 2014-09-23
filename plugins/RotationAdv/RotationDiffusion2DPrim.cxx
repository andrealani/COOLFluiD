// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdv/RotationAdv.hh"
#include "RotationDiffusion2DPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "RotationAdv/RotationAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RotationDiffusion2DPrim, DiffusiveVarSet,
	       RotationAdvModule, 2>
rd2DPrimProvider("RotationDiffusion2DPrim");

//////////////////////////////////////////////////////////////////////////////

RotationDiffusion2DPrim::RotationDiffusion2DPrim
(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  RotationDiffusion2DVarSet(name, model),
  _RotationAdvModel(model->getConvectiveTerm().d_castTo<RotationAdvTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

RotationDiffusion2DPrim::~RotationDiffusion2DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void RotationDiffusion2DPrim::setGradientVars(const vector<RealVector*>& states,
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

void RotationDiffusion2DPrim::setGradientState(const RealVector& state)
{
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
