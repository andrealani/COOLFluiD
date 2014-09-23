// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/FwdEulerData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<FwdEulerData>, FwdEulerData, ForwardEulerLib> nullFwdEulerComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void FwdEulerData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("VarID","Variable for which to compute the Norm of the dU");
   options.addConfigOption< bool >("TimeAccurate","True if time accurate");
   options.addConfigOption< bool >("PrintHistory","Print convergence history for each iteration");
}

//////////////////////////////////////////////////////////////////////////////

FwdEulerData::FwdEulerData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner),
    m_achieved(false)
{
  addConfigOptionsTo(this);

  m_varID = 0;
  setParameter("VarID",&m_varID);

  m_printHistory = false;
  setParameter("PrintHistory",&m_printHistory);

  m_isTimeAccurate = false;
  setParameter("TimeAccurate",&m_isTimeAccurate);
}

//////////////////////////////////////////////////////////////////////////////

FwdEulerData::~FwdEulerData()
{
}

//////////////////////////////////////////////////////////////////////////////

void FwdEulerData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

