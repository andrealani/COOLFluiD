// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    FluxReconstructionStrategy,FluxReconstructionSolverData,FluxReconstructionStrategy,FluxReconstructionModule >
  fluxReconstructionStrategyProvider("FluxReconstructionStrategy");

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionStrategy::FluxReconstructionStrategy(const std::string& name) :
  FluxReconstructionSolverStrategy(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionStrategy::~FluxReconstructionStrategy()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionStrategy::compute()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionStrategy::setup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

