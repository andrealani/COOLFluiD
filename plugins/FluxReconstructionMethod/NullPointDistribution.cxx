// Copyright (C) 2016 KU Leuven, von Karman Institute for Fluid Dynamics, Belgium
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

  std::vector<CFreal> coords;
  coords.resize(solOrder);
  return coords;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NullPointDistribution::getSubcellResolution(CFPolyOrder::Type solOrder)
{
  return 2.0;
}

//////////////////////////////////////////////////////////////////////////////

void NullPointDistribution::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  BasePointDistribution::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NullPointDistribution::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  BasePointDistribution::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

