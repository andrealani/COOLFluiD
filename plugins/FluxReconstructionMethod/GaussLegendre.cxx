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
#include "FluxReconstructionMethod/GaussLegendre.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<GaussLegendre,
				  FluxReconstructionSolverData,
				  BasePointDistribution,
				  FluxReconstructionModule >
GaussLegendreStrategyProvider("GaussLegendre");

//////////////////////////////////////////////////////////////////////////////

GaussLegendre::GaussLegendre(const std::string& name) :
  BasePointDistribution(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

GaussLegendre::~GaussLegendre()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<CFreal> GaussLegendre::getLocalCoords1D(CFPolyOrder::Type solOrder)
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "GaussLegendre::getLocalCoords1D()\n");
  
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
	coords[0] = -1./sqrt(3.);
	coords[1] = +1./sqrt(3.);
      } break;
      case CFPolyOrder::ORDER2:
      {
	coords[0] = -sqrt(3./5.);
	coords[1] = 0.0;
	coords[2] = +sqrt(3./5.);
      } break;
      case CFPolyOrder::ORDER3:
      {
	coords[0] = -sqrt((3.+2.*sqrt(6./5.))/7.);
	coords[1] = -sqrt((3.-2.*sqrt(6./5.))/7.);
	coords[2] = +sqrt((3.-2.*sqrt(6./5.))/7.);
	coords[3] = +sqrt((3.+2.*sqrt(6./5.))/7.);
      } break;
      case CFPolyOrder::ORDER4:
      {
	coords[0] = -sqrt(5.+2.*sqrt(10./7.))/3.;
	coords[1] = -sqrt(5.-2.*sqrt(10./7.))/3.;
	coords[2] = 0.0;
	coords[3] = +sqrt(5.-2.*sqrt(10./7.))/3.;
	coords[4] = +sqrt(5.+2.*sqrt(10./7.))/3.;
      } break;
      case CFPolyOrder::ORDER5:
      {
        coords[0] = -0.9324695142031521;
	coords[1] = -0.6612093864662645;
	coords[2] = -0.2386191860831969;
	coords[3] = 0.2386191860831969;
	coords[4] = 0.6612093864662645;
	coords[5] = 0.9324695142031521;
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

void GaussLegendre::setup()
{
  CFAUTOTRACE;
  
  BasePointDistribution::setup();

  CFLog(VERBOSE, "GaussLegendre::setup()\n");
}

//////////////////////////////////////////////////////////////////////////////

void GaussLegendre::unsetup()
{
  CFAUTOTRACE;
  
  BasePointDistribution::unsetup();
  
  CFLog(VERBOSE, "GaussLegendre::unsetup()\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

