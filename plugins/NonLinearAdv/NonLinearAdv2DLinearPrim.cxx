// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NonLinearAdv/NonLinearAdv.hh"
#include "NonLinearAdv2DLinearPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "NonLinearAdvPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NonLinearAdv2DLinearPrim, JacobianLinearizer, NonLinearAdvModule, 1> NonLinearAdv2DLinearPrimProvider("NonLinearAdv2DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

NonLinearAdv2DLinearPrim::NonLinearAdv2DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<NonLinearAdvTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

NonLinearAdv2DLinearPrim::~NonLinearAdv2DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DLinearPrim::linearize(const vector<State*>& statesInCell)
{
  cf_assert(_model.isNotNull());

  RealVector& linearData = _model->getPhysicalData();

  const CFuint nbStates = statesInCell.size();
  CFreal mean_u = 0.0;
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    cf_assert(statesInCell[iState]->size() == 1);
    mean_u += (*statesInCell[iState])[0];  
  }

  mean_u /= nbStates;

  linearData[NonLinearAdvTerm::VX] = exp(mean_u);
  linearData[NonLinearAdvTerm::VY] = 1.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
