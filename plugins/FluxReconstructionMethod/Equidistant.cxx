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
  CFLog(VERBOSE, "Equidistant::getLocalCoords1D()\n");
  
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
	coords[0] = -1./3.;
	coords[1] = +1./3.;
      } break;
      case CFPolyOrder::ORDER2:
      {
	coords[0] = -1./2.;
	coords[1] = 0.0;
	coords[2] = +1./2.;
      } break;
      case CFPolyOrder::ORDER3:
      {
	coords[0] = -0.6;
	coords[1] = -0.2;
	coords[2] = +0.2;
	coords[3] = +0.6;
      } break;
      case CFPolyOrder::ORDER4:
      {
	coords[0] = -2./3.;
	coords[1] = -1./3.;
	coords[2] = 0.0;
	coords[3] = +1./3.;
	coords[4] = +2./3.;
      } break;
      case CFPolyOrder::ORDER5:
      {
        coords[0] = -0.7142857143;
	coords[1] = -0.4285714286;
	coords[2] = -0.1428571429;
	coords[3] = 0.1428571429;
	coords[4] = 0.4285714286;
	coords[5] = 0.7142857143;
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

void Equidistant::setup()
{
  CFAUTOTRACE;
  
  BasePointDistribution::setup();

  CFLog(VERBOSE, "Equidistant::setup()\n");
}

//////////////////////////////////////////////////////////////////////////////

void Equidistant::unsetup()
{
  CFAUTOTRACE;
  
  BasePointDistribution::unsetup();
  
  CFLog(VERBOSE, "Equidistant::unsetup()\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

