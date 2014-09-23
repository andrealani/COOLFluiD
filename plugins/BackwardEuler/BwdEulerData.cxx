// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "BackwardEuler/BackwardEuler.hh"

#include "BwdEulerData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<BwdEulerData>, BwdEulerData, BackwardEulerModule> nullBwdEulerComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void BwdEulerData::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

BwdEulerData::BwdEulerData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner),
    m_lss()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

BwdEulerData::~BwdEulerData()
{
}

//////////////////////////////////////////////////////////////////////////////

void BwdEulerData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

