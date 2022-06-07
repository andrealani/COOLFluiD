// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodStrategyProvider.hh"
#include "Common/StringOps.hh"
#include "Common/NotImplementedException.hh"
#include "Common/CFLog.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/Equidistant.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<Equidistant,
				  FluxReconstructionSolverData,
				  BasePointDistribution,
				  FluxReconstructionModule >
EquidistantStrategyProvider("Equidistant");

//////////////////////////////////////////////////////////////////////////////

Equidistant::Equidistant(const std::string& name) :
  BasePointDistribution(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

Equidistant::~Equidistant()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<CFreal> Equidistant::getLocalCoords1D(CFPolyOrder::Type solOrder)
{
  CFAUTOTRACE;
  
  // initialize return variable
  const CFuint orderPlus1 = solOrder+1;
  std::vector<CFreal> coords;
  coords.resize(orderPlus1);
  
  // for Equidistant just a formula is used for the coordinates
  for (CFuint i = 0; i < orderPlus1; ++i)
  {
    coords[i] = (2.*i - solOrder)/(solOrder + 2.);
  }

  return coords;
}

//////////////////////////////////////////////////////////////////////////////

CFreal Equidistant::getSubcellResolution(CFPolyOrder::Type solOrder)
{
  CFAUTOTRACE;

  return 2./(solOrder + 2.);
}

//////////////////////////////////////////////////////////////////////////////

void Equidistant::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  BasePointDistribution::setup();
}

//////////////////////////////////////////////////////////////////////////////

void Equidistant::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  BasePointDistribution::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

