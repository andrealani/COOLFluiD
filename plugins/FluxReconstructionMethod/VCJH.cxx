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
void VCJH::defineConfigOptions(Config::OptionList& options)
{
    options.addConfigOption< CFreal >("CFactor","Value of the C factor for VCJH 1D correction function");
}

//////////////////////////////////////////////////////////////////////////////
      
VCJH::VCJH(const std::string& name) :
  BaseCorrectionFunction(name)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  
  m_cfactor = 0;
  setParameter("CFactor",&m_cfactor);
}

//////////////////////////////////////////////////////////////////////////////

VCJH::~VCJH()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////
      
void VCJH::computeCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< RealVector > >& corrfct)
{
    CFAUTOTRACE;
    const CFGeoShape::Type elemShape = frElemData->getShape();
    const CFPolyOrder::Type solOrder = frElemData->getPolyOrder();
    switch(elemShape)
    {
      case CFGeoShape::QUAD:
      {
      CFuint iSol = 0;
      const CFuint nbrSolPnts1D = frElemData->getSolPntsLocalCoord1D()->size();
      const CFuint nbrSolPnts = frElemData->getNbrOfSolPnts();
      const CFuint nbrFlxPnts = frElemData->getNbrOfFlxPnts();
      const CFuint dim = frElemData->getDimensionality();
      std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
      corrfct.resize(nbrSolPnts);
      for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
      {
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta, ++iSol)
          {
              corrfct[iSol].resize(nbrFlxPnts);
              for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
              {
                  corrfct[iSol][iFlx].resize(dim);
              }
              corrfct[iSol][4*iEta][0] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][0], m_cfactor);
              corrfct[iSol][1+4*iEta][0] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][0], m_cfactor);
              corrfct[iSol][2+4*iKsi][1] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][1], m_cfactor);
              corrfct[iSol][3+4*iKsi][1] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][1], m_cfactor);
          }
      }
          break;
      }
      case CFGeoShape::HEXA:
      {
      CFuint iSol = 0;
      const CFuint nbrSolPnts1D = frElemData->getSolPntsLocalCoord1D()->size();
      const CFuint nbrSolPnts = frElemData->getNbrOfSolPnts();
      const CFuint nbrFlxPnts = frElemData->getNbrOfFlxPnts();
      const CFuint dim = frElemData->getDimensionality();
      std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
      corrfct.resize(nbrSolPnts);
      for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
      {
          CFuint iSolPlane = 0;
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
          {
              for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++iSol, ++iSolPlane)
              {
                  corrfct[iSol].resize(nbrFlxPnts);
                  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
                  {
                      corrfct[iSol][iFlx].resize(dim);
                  }
                  corrfct[iSol][2*nbrSolPnts1D*nbrSolPnts1D+iSolPlane][0] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][KSI], m_cfactor);
                  corrfct[iSol][3*nbrSolPnts1D*nbrSolPnts1D+iSolPlane][0] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][KSI], m_cfactor);
                  corrfct[iSol][4*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iEta+iKsi][1] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][ETA], m_cfactor);
                  corrfct[iSol][5*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iEta+iKsi][1] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][ETA], m_cfactor);
                  corrfct[iSol][iZta+nbrSolPnts1D*iKsi][2] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][ZTA], m_cfactor);
                  corrfct[iSol][iZta+nbrSolPnts1D*iKsi+nbrSolPnts1D*nbrSolPnts1D][2] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][ZTA], m_cfactor);
                  
              }
          }
      }
          break;
      }
      default:
      {
        throw Common::NotImplementedException (FromHere(),"VCJH Correction Functions not implemented for elements of type "
                                               + StringOps::to_str(elemShape) + ".");
      }
    }
}
      
//////////////////////////////////////////////////////////////////////////////

void VCJH::computeDivCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< CFreal > >& corrfct)
{
    CFAUTOTRACE;
    const CFGeoShape::Type elemShape = frElemData->getShape();
    const CFPolyOrder::Type solOrder = frElemData->getPolyOrder();
    switch(elemShape)
    {
      case CFGeoShape::QUAD:
      {
        CFuint iSol = 0;
        const CFuint nbrSolPnts1D = frElemData->getSolPntsLocalCoord1D()->size();
        const CFuint nbrSolPnts = frElemData->getNbrOfSolPnts();
        const CFuint nbrFlxPnts = frElemData->getNbrOfFlxPnts();
        const CFuint dim = frElemData->getDimensionality();
        std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
        corrfct.resize(nbrSolPnts);
        for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
        {
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta, ++iSol)
          {
            corrfct[iSol].resize(nbrFlxPnts);
            corrfct[iSol][4*iEta] = computeDerivativeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][0], m_cfactor);
            corrfct[iSol][1+4*iEta] = -computeDerivativeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][0], m_cfactor);
            corrfct[iSol][2+4*iKsi] = computeDerivativeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][1], m_cfactor);
            corrfct[iSol][3+4*iKsi] = -computeDerivativeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][1], m_cfactor);
          }
        }
        break;
      }
      case CFGeoShape::HEXA:
      {
        CFuint iSol = 0;
        const CFuint nbrSolPnts1D = frElemData->getSolPntsLocalCoord1D()->size();
        const CFuint nbrSolPnts = frElemData->getNbrOfSolPnts();
        const CFuint nbrFlxPnts = frElemData->getNbrOfFlxPnts();
        const CFuint dim = frElemData->getDimensionality();
        std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
        corrfct.resize(nbrSolPnts);
        for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
        {
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
          {
            for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++iSol)
            {
              corrfct[iSol].resize(nbrFlxPnts);
              corrfct[iSol][2*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iEta+iZta] = computeDerivativeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][KSI], m_cfactor);
              corrfct[iSol][3*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iEta+iZta] = -computeDerivativeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][KSI], m_cfactor);
              corrfct[iSol][4*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iZta+iKsi] = computeDerivativeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][ETA], m_cfactor);
              corrfct[iSol][5*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iZta+iKsi] = -computeDerivativeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][ETA], m_cfactor);
              corrfct[iSol][0*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iKsi+iEta] = computeDerivativeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][ZTA], m_cfactor);
              corrfct[iSol][1*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iKsi+iEta] = -computeDerivativeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][ZTA], m_cfactor);
            }
	  }
        }
        break;
      }
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Divergence of VCJH Correction Functions not implemented for elements of type "
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
            corrfct = -0.5*(ksi-(1.5*cfactor+0.5*(3.*pow(ksi,2.)-1.))/(1.+1.5*cfactor));
        } break;
        case CFPolyOrder::ORDER2:
        {
            corrfct = 0.5*(0.5*(3.*pow(ksi,2.0)-1.)-(22.5*cfactor*ksi+0.5*(5.*pow(ksi,3.)-3.*ksi))/(1.+22.5*cfactor));
        } break;
        case CFPolyOrder::ORDER3:
        {
            corrfct = -0.5*(0.5*(5.*pow(ksi,3.)-3.*ksi)-(787.5*cfactor*0.5*(3.*pow(ksi,2.)-1.)+0.125*(35.*pow(ksi,4.)-30.*pow(ksi,2.)+3.))/(1.+787.5*cfactor));
        } break;
        case CFPolyOrder::ORDER4:
        {
            corrfct = 0.5*(0.125*(35.*pow(ksi,4.)-30.*pow(ksi,2.)+3.)-(49612.5*cfactor*0.5*(5.*pow(ksi,3.)-3.*ksi)+0.125*(63.*pow(ksi,5.)-70.*pow(ksi,3.)+15.*ksi))/(1.+49612.5*cfactor));
        } break;
        case CFPolyOrder::ORDER5:
        {
            corrfct = -0.5*(0.125*(63.*pow(ksi,5.)-70.*pow(ksi,3.)+15.*ksi)-(4911637.5*cfactor*0.125*(35.*pow(ksi,4.)-30.*pow(ksi,2.)+3.)+0.0625*(231.*pow(ksi,6.)-315.*pow(ksi,4.)+105.*pow(ksi,2.)-5.))/(1.+4911637.5*cfactor));
        } break;
        default:
        {
            throw Common::NotImplementedException (FromHere(),"VCJH1D Correction Function not implemented for order " + StringOps::to_str(solOrder) + ".");
        }
    }
    return corrfct;
}
      
//////////////////////////////////////////////////////////////////////////////

CFreal VCJH::computeDerivativeCorrectionFunction1D(CFPolyOrder::Type solOrder, CFreal ksi, CFreal cfactor)
{
    CFAUTOTRACE;
    CFreal corrfct;
    switch(solOrder)
    {
        case CFPolyOrder::ORDER1:
        {
            corrfct = -0.5+1.5*ksi/(1.+1.5*cfactor);
        } break;
        case CFPolyOrder::ORDER2:
        {
            corrfct = 0.5*(0.5*6.*ksi-(22.5*cfactor+0.5*(15.*pow(ksi,2.)-3.))/(1.+22.5*cfactor));
        } break;
        case CFPolyOrder::ORDER3:
        {
            corrfct = -0.5*(0.5*(15.*pow(ksi,2.)-3.)-(787.5*cfactor*0.5*6.*ksi+0.125*(140.*pow(ksi,3.)-60.*ksi))/(1.+787.5*cfactor));
        } break;
        case CFPolyOrder::ORDER4:
        {
            corrfct = 0.5*(0.125*(140.*pow(ksi,3.)-60.*ksi)-(49612.5*cfactor*0.5*(15.*pow(ksi,2.)-3.)+0.125*(315.*pow(ksi,4.)-210.*pow(ksi,2.)+15.))/(1.+49612.5*cfactor));
        } break;
        case CFPolyOrder::ORDER5:
        {
            corrfct = -0.5*(0.125*(315.*pow(ksi,4.)-210.*pow(ksi,2.)+15.)-(4911637.5*cfactor*0.125*(140.*pow(ksi,3.)-60.*ksi)+0.0625*(1386.*pow(ksi,5.)-1260.*pow(ksi,3.)+210.*ksi))/(1.+4911637.5*cfactor));
        } break;
        default:
        {
            throw Common::NotImplementedException (FromHere(),"Derivative of VCJH1D Correction Function not implemented for order " + StringOps::to_str(solOrder) + ".");
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
      
void VCJH::configure ( Config::ConfigArgs& args )
{
    FluxReconstructionSolverStrategy::configure(args);
}
      
//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
