// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullLinearSystemSolver.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< NullLinearSystemSolver,
                LinearSystemSolver,
                FrameworkLib,
                1 >

nullLinearSystemSolverProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullLinearSystemSolver::NullLinearSystemSolver(const std::string& name)
  : LinearSystemSolver(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullLinearSystemSolver::~NullLinearSystemSolver()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NullLinearSystemSolver::getMethodData () const
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void NullLinearSystemSolver::setMethodImpl()
{
  LinearSystemSolver::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullLinearSystemSolver::unsetMethodImpl()
{
  LinearSystemSolver::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullLinearSystemSolver::solveSysImpl()
{
  CFLog(VERBOSE,"NullLinearSystemSolver::solveSys() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullLinearSystemSolver::printToFile(const std::string prefix, const std::string suffix)
{
  CFLog(VERBOSE,"NullLinearSystemSolver::printToFile() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

BlockAccumulator* NullLinearSystemSolver::createBlockAccumulator(const CFuint nbRows,
								 const CFuint nbCols,
								 const CFuint subBlockSize, 
								 CFreal* ptr) const
{
  CFLog(VERBOSE,"NullLinearSystemSolver::createBlockAccumulator() called!" << "\n");
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<LSSMatrix> NullLinearSystemSolver::getMatrix() const
{
  CFLog(VERBOSE,"NullLinearSystemSolver::getMatrix() called!" << "\n");
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<LSSVector> NullLinearSystemSolver::getSolVector() const
{
  CFLog(VERBOSE,"NullLinearSystemSolver::getSolVector() called!" << "\n");
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<LSSVector> NullLinearSystemSolver::getRhsVector() const
{
  CFLog(VERBOSE,"NullLinearSystemSolver::getRhsVector() called!" << "\n");
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
