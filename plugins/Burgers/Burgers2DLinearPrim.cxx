// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Burgers/Burgers.hh"
#include "Burgers2DLinearPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"
#include "BurgersPhysicalModel.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Burgers {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Burgers2DLinearPrim, JacobianLinearizer, BurgersModule, 1> burgers2DLinearPrimProvider("Burgers2DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

Burgers2DLinearPrim::Burgers2DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().d_castTo<BurgersTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Burgers2DLinearPrim::~Burgers2DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Burgers2DLinearPrim::linearize(const std::vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();

  const CFuint nbStates = statesInCell.size();
  CFreal sum = 0.;
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    cf_assert(statesInCell[iState]->size() == 1);
    sum += (*statesInCell[iState])[0];
  }
  sum *= 1./nbStates;

  linearData[BurgersTerm::VX] = sum;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Burgers

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
