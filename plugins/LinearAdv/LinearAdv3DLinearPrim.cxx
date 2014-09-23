// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdv/LinearAdv.hh"
#include "LinearAdv3DLinearPrim.hh"
#include "LinearAdvTerm.hh"
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

Environment::ObjectProvider<LinearAdv3DLinearPrim, JacobianLinearizer, LinearAdvModule, 1> linearAdv3DLinearPrimProvider("LinearAdv3DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

LinearAdv3DLinearPrim::LinearAdv3DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model((model->getImplementor()->getConvectiveTerm().d_castTo<LinearAdvTerm>()))
{
}

//////////////////////////////////////////////////////////////////////////////

LinearAdv3DLinearPrim::~LinearAdv3DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdv3DLinearPrim::linearize(const std::vector<State*>& statesInCell)
{
  cf_assert(_model.isNotNull());

  RealVector& linearData = _model->getPhysicalData();

  const std::vector<RealMatrix> *const jacob =
    m_physmodel->getImplementor()->getJacobians();

  cf_assert(jacob != CFNULL);

  linearData[LinearAdvTerm::VX] = ((*jacob)[XX])(0,0);
  linearData[LinearAdvTerm::VY] = ((*jacob)[YY])(0,0);
  linearData[LinearAdvTerm::VZ] = ((*jacob)[ZZ])(0,0);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
