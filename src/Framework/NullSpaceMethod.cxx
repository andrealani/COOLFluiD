// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullSpaceMethod.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullSpaceMethod,
               SpaceMethod,
               FrameworkLib,
               1>
nullSpaceMethodProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullSpaceMethod::NullSpaceMethod(const std::string& name)
  : SpaceMethod(name)
{
  m_builder = "Null";
}

//////////////////////////////////////////////////////////////////////////////

NullSpaceMethod::~NullSpaceMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::initializeSolution ()
{
  CFLog(VERBOSE,"NullSpaceMethod::initializeSolution() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NullSpaceMethod::getMethodData () const
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<SpaceMethodData> NullSpaceMethod::getSpaceMethodData()
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::setMethodImpl()
{
  SpaceMethod::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::unsetMethodImpl()
{
  SpaceMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::setCollaborator(MultiMethodHandle<LinearSystemSolver> lss)
{
  CFLog(VERBOSE,"NullSpaceMethod::setLinearSystemSolver() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  CFLog(VERBOSE,"NullSpaceMethod::setConvergenceMethod() called!" << "\n");
}


//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::initializeSolutionImpl(bool isRestart)
{
  CFLog(VERBOSE,"NullSpaceMethod::initializeSolution() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::extrapolateStatesToNodesImpl()
{
  CFLog(VERBOSE,"NullSpaceMethod::setNodalStates() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::prepareComputationImpl()
{
  CFLog(VERBOSE,"NullSpaceMethod::prepareComputation() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::computeSpaceResidualImpl(CFreal factor)
{
  CFLog(VERBOSE,"NullSpaceMethod::computeResidual() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::computeTimeResidualImpl(CFreal factor)
{
  CFLog(VERBOSE,"NullSpaceMethod::computeTimeResidual() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t NullSpaceMethod::beforeMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFLog(VERBOSE,"NullSpaceMethod::beforeMeshUpdateAction() called!" << "\n");
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t NullSpaceMethod::afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
{
  CFLog(VERBOSE,"NullSpaceMethod::afterMeshUpdateAction() called!" << "\n");
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

void NullSpaceMethod::applyBCImpl()
{
  CFLog(VERBOSE,"NullSpaceMethod::computeApplyBC() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
