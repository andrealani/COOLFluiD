// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/NullPointDistribution.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<NullPointDistribution,
				  FluxReconstructionSolverData,
				  BasePointDistribution,
				  FluxReconstructionModule >
nullPointDistributionStrategyProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullPointDistribution::NullPointDistribution(const std::string& name) :
  BasePointDistribution(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NullPointDistribution::~NullPointDistribution()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<CFreal> NullPointDistribution::getLocalCoords1D(CFPolyOrder::Type solOrder)
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "NullPointDistribution::getLocalCoords1D()\n");
  std::vector<CFreal> coords;
  coords.resize(solOrder);
  return coords;
}

//////////////////////////////////////////////////////////////////////////////

void NullPointDistribution::setup()
{
  CFAUTOTRACE;
  
  BasePointDistribution::setup();

  CFLog(VERBOSE, "NullPointDistribution::setup()\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullPointDistribution::unsetup()
{
  CFAUTOTRACE;
  
  BasePointDistribution::unsetup();
  
  CFLog(VERBOSE, "NullPointDistribution::unsetup()\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

