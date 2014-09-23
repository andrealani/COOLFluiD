// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearizedBDF2.hh"
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

Environment::ObjectProvider<LinearizedBDF2,
               ConvergenceMethod,
               NewtonMethodModule,
               1>
LinearizedBDF2MethodProvider("LinearizedBDF2");

//////////////////////////////////////////////////////////////////////////////

LinearizedBDF2::LinearizedBDF2(const std::string& name)
  : BDF2(name)
{
  // new default commands are set for the LinearizedBDF2 scheme
  m_setupStr   = "LinearizedBDF2Setup";
  m_unSetupStr = "LinearizedUnSetup";
  m_prepareStr = "LinearizedPrepare";
  m_intermediateStr = "BDF2Intermediate";

}

//////////////////////////////////////////////////////////////////////////////

LinearizedBDF2::~LinearizedBDF2()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

