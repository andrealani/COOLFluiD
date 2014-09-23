// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdv/LinearAdv.hh"
#include "LinearAdv2DLinearPrim.hh"
#include "LinearAdvPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinearAdv2DLinearPrim, JacobianLinearizer, LinearAdvModule, 1> LinearAdv2DLinearPrimProvider("LinearAdv2DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

LinearAdv2DLinearPrim::LinearAdv2DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model((model->getImplementor()->getConvectiveTerm().d_castTo<LinearAdvTerm>()))
{
}

//////////////////////////////////////////////////////////////////////////////

LinearAdv2DLinearPrim::~LinearAdv2DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdv2DLinearPrim::linearize(const std::vector<State*>& statesInCell)
{
  cf_assert(_model.isNotNull());

  RealVector& linearData = _model->getPhysicalData();

  const std::vector<RealMatrix> *const jacob = m_physmodel->getImplementor()->getJacobians();

  cf_assert(jacob != CFNULL);

  linearData[LinearAdvTerm::VX] = ((*jacob)[XX])(0,0);
  linearData[LinearAdvTerm::VY] = ((*jacob)[YY])(0,0);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
}
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
