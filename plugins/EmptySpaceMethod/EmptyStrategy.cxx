// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodStrategyProvider.hh"

#include "EmptySpaceMethod/Empty.hh"
#include "EmptySpaceMethod/EmptyStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    EmptyStrategy,EmptySolverData,EmptyStrategy,EmptyModule >
  emptyStrategyProvider("EmptyStrategy");

//////////////////////////////////////////////////////////////////////////////

EmptyStrategy::EmptyStrategy(const std::string& name) :
  EmptySolverStrategy(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

EmptyStrategy::~EmptyStrategy()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void EmptyStrategy::compute()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void EmptyStrategy::setup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod

}  // namespace COOLFluiD

