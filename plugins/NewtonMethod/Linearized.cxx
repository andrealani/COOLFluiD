// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Linearized.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "Framework/CFL.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Linearized,
               ConvergenceMethod,
               NewtonMethodModule,
               1>
LinearizedMethodProvider("Linearized");

//////////////////////////////////////////////////////////////////////////////

Linearized::Linearized(const std::string& name)
  : CrankNicholson(name)
{
  // new default commands are set for the Linearized scheme
  m_setupStr   = "LinearizedSetup";
  m_unSetupStr = "LinearizedUnSetup";
  m_prepareStr = "LinearizedPrepare";
  m_intermediateStr = "LinearizedIntermediate";
}

//////////////////////////////////////////////////////////////////////////////

Linearized::~Linearized()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

