// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewmarkImplicit.hh"
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

Environment::ObjectProvider<NewmarkImplicit,
                            ConvergenceMethod,
                            NewtonMethodModule,
                            1>
NewmarkImplicitMethodProvider("NewmarkImplicit");

//////////////////////////////////////////////////////////////////////////////

NewmarkImplicit::NewmarkImplicit(const std::string& name)
  : NewtonIterator(name)
{
  // new default commands are set for the NewmarkImplicit scheme
  m_initStr = "ResetSystem";
  m_setupStr   = "NewmarkSetup";
  m_unSetupStr = "NewmarkUnSetup";
  m_prepareStr = "NewmarkPrepare";
  m_updateSolStr = "NewmarkImplicitUpdateSol";
}

//////////////////////////////////////////////////////////////////////////////

NewmarkImplicit::~NewmarkImplicit()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

