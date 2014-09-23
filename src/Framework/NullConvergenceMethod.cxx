// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullConvergenceMethod.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullConvergenceMethod,
               ConvergenceMethod,
               FrameworkLib,
               1>
nullConvergenceMethodProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullConvergenceMethod::NullConvergenceMethod(const std::string& name)
  : ConvergenceMethod(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullConvergenceMethod::~NullConvergenceMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NullConvergenceMethod::getMethodData () const
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void NullConvergenceMethod::setMethodImpl()
{
  ConvergenceMethod::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullConvergenceMethod::unsetMethodImpl()
{
  ConvergenceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullConvergenceMethod::takeStepImpl()
{
  CFLog(VERBOSE,"NullConvergenceMethod::takeStep() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConvergenceMethodData> NullConvergenceMethod::getConvergenceMethodData()
{
  return CFNULL;
}
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
