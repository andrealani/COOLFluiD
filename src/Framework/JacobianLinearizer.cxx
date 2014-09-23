// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/JacobianLinearizer.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

JacobianLinearizer::JacobianLinearizer(Common::SafePtr<Framework::PhysicalModel> model) :
  Common::OwnedObject(),
  _sumZ(model->getNbEq()),
  _avZ(model->getNbEq()),
  _maxNbStates(0),
  _extraValuesDelta(CFNULL),
  _upStates(CFNULL),
  m_physmodel(model)
{
}

//////////////////////////////////////////////////////////////////////////////

JacobianLinearizer::~JacobianLinearizer()
{
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
