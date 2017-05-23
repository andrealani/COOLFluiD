// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodStrategyProvider.hh"
#include "Common/StringOps.hh"
#include "Common/NotImplementedException.hh"
#include "Common/CFLog.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/Lobatto.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<Lobatto,
				  FluxReconstructionSolverData,
				  BasePointDistribution,
				  FluxReconstructionModule >
LobattoStrategyProvider("Lobatto");

//////////////////////////////////////////////////////////////////////////////

Lobatto::Lobatto(const std::string& name) :
  BasePointDistribution(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

Lobatto::~Lobatto()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<CFreal> Lobatto::getLocalCoords1D(CFPolyOrder::Type solOrder)
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "Lobatto::getLocalCoords1D()\n");
  
  std::vector<CFreal> coords;
  coords.resize(solOrder+1);
  switch(solOrder)
    {
      case CFPolyOrder::ORDER0:
      {
	coords[0] = 0.0;
      } break;
      case CFPolyOrder::ORDER1:
      {
	coords[0] = -1.;
	coords[1] = +1.;
      } break;
      case CFPolyOrder::ORDER2:
      {
	coords[0] = -1.;
	coords[1] = 0.0;
	coords[2] = +1.;
      } break;
      case CFPolyOrder::ORDER3:
      {
	coords[0] = -1.;
	coords[1] = -0.4472135954999579392818;
	coords[2] = +0.4472135954999579392818;
	coords[3] = +1.;
      } break;
      case CFPolyOrder::ORDER4:
      {
	coords[0] = -1.;
	coords[1] = -0.6546536707079771437983;
	coords[2] = 0.0;
	coords[3] = +0.6546536707079771437983;
	coords[4] = +1.;
      } break;
      case CFPolyOrder::ORDER5:
      {
        coords[0] = -1.0;
	coords[1] = -0.765055323929464692851;
	coords[2] = -0.2852315164806450963142;
	coords[3] = 0.2852315164806450963142;
	coords[4] = 0.765055323929464692851;
	coords[5] = 1.;
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Gauss Legendre not implemented for order "
                                      + StringOps::to_str(solOrder) + ".");
      }
    }

  return coords;
}

//////////////////////////////////////////////////////////////////////////////

void Lobatto::setup()
{
  CFAUTOTRACE;
  
  BasePointDistribution::setup();

  CFLog(VERBOSE, "Lobatto::setup()\n");
}

//////////////////////////////////////////////////////////////////////////////

void Lobatto::unsetup()
{
  CFAUTOTRACE;
  
  BasePointDistribution::unsetup();
  
  CFLog(VERBOSE, "Lobatto::unsetup()\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

