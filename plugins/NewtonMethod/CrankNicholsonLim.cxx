// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CrankNicholsonLim.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/CFL.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CrankNicholsonLim,
               ConvergenceMethod,
               NewtonMethodModule,
               1>
crankNichsolsonLimMethodProvider("CrankNicholsonLim");

//////////////////////////////////////////////////////////////////////////////

CrankNicholsonLim::CrankNicholsonLim(const std::string& name)
  : CrankNicholson(name)
{
  // new default commands are set for the Limited Crank-Nicholson scheme
  m_setupStr   = "CrankNichLimSetup";
  m_unSetupStr = "CrankNichLimUnSetup";
  m_intermediateStr = "CrankNichLimIntermediate";
  m_prepareStr = "CrankNichLimPrepare";
  m_initStr = "CrankNichLimInit";
}

//////////////////////////////////////////////////////////////////////////////

CrankNicholsonLim::~CrankNicholsonLim()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

