// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/NullFluxPntDistribution.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<NullFluxPntDistribution,
				  FluxReconstructionSolverData,
				  BaseFluxPntDistribution,
				  FluxReconstructionModule >
nullFluxPntDistributionStrategyProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullFluxPntDistribution::NullFluxPntDistribution(const std::string& name) :
  BaseFluxPntDistribution(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NullFluxPntDistribution::~NullFluxPntDistribution()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void NullFluxPntDistribution::compute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "NullFluxPntDistribution::compute()\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullFluxPntDistribution::setup()
{
  CFAUTOTRACE;
  
  BaseFluxPntDistribution::setup();

  CFLog(VERBOSE, "NullFluxPntDistribution::setup()\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullFluxPntDistribution::unsetup()
{
  CFAUTOTRACE;
  
  BaseFluxPntDistribution::unsetup();
  
  CFLog(VERBOSE, "NullFluxPntDistribution::unsetup()\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

