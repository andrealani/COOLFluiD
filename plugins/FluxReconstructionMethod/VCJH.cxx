#include <iterator>

#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/VCJH.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    VCJH,FluxReconstructionSolverData,BaseCorrectionFunction,FluxReconstructionModule >
VCJHProvider("VCJH");

//////////////////////////////////////////////////////////////////////////////

VCJH::VCJH(const std::string& name) :
  BaseCorrectionFunction(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

VCJH::~VCJH()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////
      
void VCJH::computeCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< RealVector > > corrcts)
{
    CFAUTOTRACE;
    const CFGeoShape::Type elemShape = frElemData->getShape();
    switch(elemShape)
    {
      case CFGeoShape::QUAD:
      {
      const CFuint nbrSolPnts1D = frElemData->getSolPntsLocalCoord1D()->size();
      for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
      {
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
          {
              
          }
      }
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"VCJH Correction Functions not implemented for elements of type "
                                               + StringOps::to_str(elemShape) + ".");
      }
    }
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFreal VCJH::computeCorrectionFunction1D(CFPolyOrder::Type solOrder, CFreal ksi, CFreal cfactor)
{
    CFAUTOTRACE;
    CFreal corrfct;
    switch(solOrder)
    {
        case CFPolyOrder::ORDER1:
        {
            corrfct = -0.5*(ksi-(1.5*cfactor+0.5*(3.*ksi*ksi-1.))/(1.+22.5*cfactor));
        } break;
        case CFPolyOrder::ORDER2:
        {
            corrfct = 0.5*(0.5*(3.*ksi*ksi-1.)-(22.5*cfactor*ksi+0.5*(5.*ksi*ksi*ksi-3.*ksi))/(1.+22.5*cfactor));
        } break;
        case CFPolyOrder::ORDER3:
        {
            corrfct = -0.5*(0.5*(5.*ksi*ksi*ksi-3.*ksi)-(787.5*cfactor*0.5*(3.*ksi*ksi-1.)+0.125*(35.*ksi*ksi*ksi*ksi-30.*ksi*ksi+3.))/(1.+787.5*cfactor));
        } break;
        case CFPolyOrder::ORDER4:
        {
            corrfct = 0.5*(0.125*(35.*ksi*ksi*ksi*ksi-30.*ksi*ksi+3.)-(49612.5*cfactor*0.5*(5.*ksi*ksi*ksi-3.*ksi)+0.125*(63.*ksi*ksi*ksi*ksi*ksi-70.*ksi*ksi*ksi+15.*ksi))/(1.+49612.5*cfactor));
        } break;
        case CFPolyOrder::ORDER5:
        {
            corrfct = -0.5*(0.125*(63.*ksi*ksi*ksi*ksi*ksi-70.*ksi*ksi*ksi+15.*ksi)-(4911637.5*cfactor*0.125*(35.*ksi*ksi*ksi*ksi-30.*ksi*ksi+3.)+0.0625*(231.*ksi*ksi*ksi*ksi*ksi*ksi-315.*ksi*ksi*ksi*ksi+105.*ksi*ksi-5.))/(1.+4911637.5*cfactor));
        } break;
        default:
        {
            throw Common::NotImplementedException (FromHere(),"VCJH1D Correction Function not implemented for order " + StringOps::to_str(solOrder) + ".");
        }
    }
    return corrfct;
}
      
//////////////////////////////////////////////////////////////////////////////

void VCJH::setup()
{
  CFAUTOTRACE;

  BaseCorrectionFunction::setup();

}

void VCJH::unsetup()
{
  CFAUTOTRACE;

  BaseCorrectionFunction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
