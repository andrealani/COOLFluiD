// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullDataProcessing.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullDataProcessing,
	       DataProcessingMethod,
               FrameworkLib,
	       1>
nullDataProcessingProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullDataProcessing::NullDataProcessing(const std::string& name)
  : DataProcessingMethod(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullDataProcessing::~NullDataProcessing()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NullDataProcessing::getMethodData () const
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void NullDataProcessing::setMethodImpl()
{
  DataProcessingMethod::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullDataProcessing::unsetMethodImpl()
{
  DataProcessingMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullDataProcessing::setCollaborator(MultiMethodHandle<SpaceMethod> spaceMtd)
{
  CFLog(VERBOSE,"NullDataProcessing::setSpaceMethod() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullDataProcessing::setCollaborator(MultiMethodHandle<ConvergenceMethod> convergenceMtd)
{
  CFLog(VERBOSE,"NullDataProcessing::setConvergenceMethod() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullDataProcessing::processDataImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
