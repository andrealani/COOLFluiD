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
      case CFPolyOrder::ORDER6:
      {
        coords[0] = -0.949107912342759;
	coords[1] = -0.741531185599394;
	coords[2] = -0.405845151377397;
	coords[3] = 0.0;
	coords[4] = 0.405845151377397;
	coords[5] = 0.741531185599394;
	coords[6] = 0.949107912342759;
      } break;
      case CFPolyOrder::ORDER7:
      {
        coords[0] = -0.960289856497536;
	coords[1] = -0.796666477413627;
	coords[2] = -0.525532409916329;
	coords[3] = -0.183434642495650;
	coords[4] = 0.183434642495650;
	coords[5] = 0.525532409916329;
	coords[6] = 0.796666477413627;
	coords[7] = 0.960289856497536;
      } break;
      case CFPolyOrder::ORDER8:
      {
        coords[0] = -0.968160239507626;
	coords[1] = -0.836031107326636;
	coords[2] = -0.613371432700590;
	coords[3] = -0.324253423403809;
	coords[4] = 0.0;
	coords[5] = 0.324253423403809;
	coords[6] = 0.613371432700590;
	coords[7] = 0.836031107326636;
	coords[8] = 0.968160239507626;
      } break;
      case CFPolyOrder::ORDER9:
      {
        coords[0] = -0.973906528517172;
	coords[1] = -0.865063366688985;
	coords[2] = -0.679409568299024;
	coords[3] = -0.433395394129247;
	coords[4] = -0.148874338981631;
	coords[5] = 0.148874338981631;
	coords[6] = 0.433395394129247;
	coords[7] = 0.679409568299024;
	coords[8] = 0.865063366688985;
	coords[9] = 0.973906528517172;
      } break;
      case CFPolyOrder::ORDER10:
      {
        coords[0] = -0.978228658146057;
	coords[1] = -0.887062599768095;
	coords[2] = -0.730152005574049;
	coords[3] = -0.519096129110681;
	coords[4] = -0.269543155952345;
	coords[5] = 0.0;
	coords[6] = 0.269543155952345;
	coords[7] = 0.519096129110681;
	coords[8] = 0.730152005574049;
	coords[9] = 0.887062599768095;
	coords[10] = 0.978228658146057;
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
  
  // setup parent class
  BasePointDistribution::setup();
}

//////////////////////////////////////////////////////////////////////////////

void GaussLegendre::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  BasePointDistribution::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

