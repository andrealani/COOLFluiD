// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"
#include "Framework/NullCouplerMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullCouplerMethod,
               CouplerMethod,
               FrameworkLib,
               1>
nullCouplerMethodProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullCouplerMethod::NullCouplerMethod(const std::string& name)
  : CouplerMethod(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullCouplerMethod::~NullCouplerMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NullCouplerMethod::getMethodData () const
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::setMethodImpl()
{
  CouplerMethod::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::unsetMethodImpl()
{
  CouplerMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::preProcessReadImpl()
{
  CFLogDebugMed("NullCouplerMethod::preProcessReadImpl() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::preProcessWriteImpl()
{
  CFLogDebugMed("NullCouplerMethod::preProcessWriteImpl() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::meshMatchingReadImpl()
{
  CFLogDebugMed("NullCouplerMethod::meshMatchingReadImpl() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::meshMatchingWriteImpl()
{
  CFLogDebugMed("NullCouplerMethod::meshMatchingWriteImpl() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::dataTransferReadImpl()
{
  CFLogDebugMed("NullCouplerMethod::dataTransferReadImpl() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::dataTransferWriteImpl()
{
  CFLogDebugMed("NullCouplerMethod::dataTransferWriteImpl() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullCouplerMethod::finalizeImpl()
{
  CFLogDebugMed("NullCouplerMethod::finalizeImpl() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
