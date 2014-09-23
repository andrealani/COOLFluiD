// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "EmptySpaceMethod/StdSolve.hh"
#include "EmptySpaceMethod/Empty.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSolve,EmptySolverData,EmptyModule >
  stdSolveProvider("StdSolve");

//////////////////////////////////////////////////////////////////////////////

void StdSolve::execute()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

