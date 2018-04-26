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
	coords[5] = 1.0;
      } break;
      case CFPolyOrder::ORDER6:
      {
        coords[0] = -1.0;
	coords[1] = -0.830223896278566929872;
	coords[2] = -0.4688487934707142138038;
	coords[3] = 0.0;
	coords[4] = 0.4688487934707142138038;
	coords[5] = 0.830223896278566929872;
	coords[6] = 1.0;
      } break;
      case CFPolyOrder::ORDER7:
      {
        coords[0] = -1.0;
	coords[1] = -0.8717401485096066153375;
	coords[2] = -0.5917001814331423021445;
	coords[3] = -0.2092992179024788687687;
	coords[4] = 0.2092992179024788687687;
	coords[5] = 0.5917001814331423021445;
	coords[6] = 0.8717401485096066153375;
	coords[7] = 1.0;
      } break;
      case CFPolyOrder::ORDER8:
      {
        coords[0] = -1.0;
	coords[1] = -0.8997579954114601573124;
	coords[2] = -0.6771862795107377534459;
	coords[3] = -0.3631174638261781587108;
	coords[4] = 0.0;
	coords[5] = 0.3631174638261781587108;
	coords[6] = 0.6771862795107377534459;
	coords[7] = 0.8997579954114601573124;
	coords[8] = 1.0;
      } break;
      case CFPolyOrder::ORDER9:
      {
        coords[0] = -1.0;
	coords[1] = -0.9195339081664588138289;
	coords[2] = -0.7387738651055050750031;
	coords[3] = -0.4779249498104444956612;
	coords[4] = -0.1652789576663870246262;
	coords[5] = 0.1652789576663870246262;
	coords[6] = 0.4779249498104444956612;
	coords[7] = 0.7387738651055050750031;
	coords[8] = 0.9195339081664588138289;
	coords[9] = 1.0;
      } break;
      case CFPolyOrder::ORDER10:
      {
        coords[0] = -1.0;
	coords[1] = -0.9340014304080591343323;
	coords[2] = -0.7844834736631444186224;
	coords[3] = -0.565235326996205006471;
	coords[4] = -0.2957581355869393914319;
	coords[5] = 0.0;
	coords[6] = 0.2957581355869393914319;
	coords[7] = 0.565235326996205006471;
	coords[8] = 0.7844834736631444186224;
	coords[9] = 0.9340014304080591343323;
	coords[10] = 1.0;
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

CFreal Lobatto::getSubcellResolution(CFPolyOrder::Type solOrder)
{
  CFAUTOTRACE;
  
  CFreal result;

  switch(solOrder)
    {
      case CFPolyOrder::ORDER0:
      {
	result = 1.;
      } break;
      case CFPolyOrder::ORDER1:
      {
	result = 2.;
      } break;
      case CFPolyOrder::ORDER2:
      {
	result = 1.;
      } break;
      case CFPolyOrder::ORDER3:
      {
	result = 2.*0.4472135954999579392818;
      } break;
      case CFPolyOrder::ORDER4:
      {
	result = 0.6546536707079771437983;
      } break;
      case CFPolyOrder::ORDER5:
      {
	result = 2.*0.2852315164806450963142;
      } break;
      case CFPolyOrder::ORDER6:
      {
	result = 0.4688487934707142138038;
      } break;
      case CFPolyOrder::ORDER7:
      {
	result = 2.*0.2092992179024788687687;
      } break;
      case CFPolyOrder::ORDER8:
      {
	result = 0.3631174638261781587108;
      } break;
      case CFPolyOrder::ORDER9:
      {
	result = 2.*0.1652789576663870246262;
      } break;
      case CFPolyOrder::ORDER10:
      {
	result = 0.2957581355869393914319;
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Gauss Legendre not implemented for order "
                                      + StringOps::to_str(solOrder) + ".");
      }
    }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Lobatto::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  BasePointDistribution::setup();
}

//////////////////////////////////////////////////////////////////////////////

void Lobatto::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  BasePointDistribution::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

