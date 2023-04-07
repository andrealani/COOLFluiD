#include <iterator>
#include <sstream>
#include <fstream>
#include <string>


#include "Common/OSystem.hh"
#include "Common/StringOps.hh" 

#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/VCJH.hh"

#include "MathTools/RealMatrix.hh"
#include "MathTools/MathConsts.hh"

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
    //options.addConfigOption< CFreal >("kappa","Value of the k factor for VCJH 2D correction function");
}

//////////////////////////////////////////////////////////////////////////////
      
VCJH::VCJH(const std::string& name) :
  BaseCorrectionFunction(name)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  
  m_cfactor = 0.0;
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
    
    // get the element shape and order
    const CFGeoShape::Type elemShape = frElemData->getShape();
    const CFPolyOrder::Type solOrder = frElemData->getPolyOrder();
    
    // compute corr fct for correct element shape
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
	cf_assert(corrfct.size() == nbrSolPnts);
	
	// loop over sol pnt coordinates in the standard domain
        for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
        {
            for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta, ++iSol)
            {
	      cf_assert(corrfct[iSol].size() == nbrFlxPnts);
              cf_assert(corrfct[iSol][0].size() == dim);
              
              // add the contriubutions of the flx pnts related to the current sol pnt
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
        cf_assert(corrfct.size() == nbrSolPnts);
	
	// loop over sol pnt coordinates in the standard domain
        for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
        {
          CFuint iSolPlane = 0;
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
          {
            for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++iSol, ++iSolPlane)
            {
              cf_assert(corrfct[iSol].size() == nbrFlxPnts);
	      cf_assert(corrfct[iSol][0].size() == dim);
              
              // add the contriubutions of the flx pnts related to the current sol pnt
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
      
void VCJH::computeCorrectionFunction(const CFPolyOrder::Type solOrder, const CFreal factor, Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< CFreal > >& corrfct)
{
    CFAUTOTRACE;
    
    // get the element shape and order
    const CFGeoShape::Type elemShape = frElemData->getShape();
    
    // compute corr fct for correct element shape
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
	cf_assert(corrfct.size() == nbrSolPnts);
	
	// loop over sol pnt coordinates in the standard domain
        for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
        {
            for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta, ++iSol)
            {
	      cf_assert(corrfct[iSol].size() == nbrFlxPnts);
              
              // add the contriubutions of the flx pnts related to the current sol pnt
              corrfct[iSol][4*iEta] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][0], factor);
              corrfct[iSol][1+4*iEta] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][0], factor);
              corrfct[iSol][2+4*iKsi] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][1], factor);
              corrfct[iSol][3+4*iKsi] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][1], factor);
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
        cf_assert(corrfct.size() == nbrSolPnts);
	
	// loop over sol pnt coordinates in the standard domain
        for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
        {
          CFuint iSolPlane = 0;
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
          {
            for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++iSol, ++iSolPlane)
            {
              cf_assert(corrfct[iSol].size() == nbrFlxPnts);
              
              // add the contriubutions of the flx pnts related to the current sol pnt
              corrfct[iSol][2*nbrSolPnts1D*nbrSolPnts1D+iSolPlane] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][KSI], factor);
              corrfct[iSol][3*nbrSolPnts1D*nbrSolPnts1D+iSolPlane] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][KSI], factor);
              corrfct[iSol][4*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iEta+iKsi] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][ETA], factor);
              corrfct[iSol][5*nbrSolPnts1D*nbrSolPnts1D+nbrSolPnts1D*iEta+iKsi] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][ETA], factor);
              corrfct[iSol][iZta+nbrSolPnts1D*iKsi] = computeCorrectionFunction1D(solOrder, solPntsLocalCoord[iSol][ZTA], factor);
              corrfct[iSol][iZta+nbrSolPnts1D*iKsi+nbrSolPnts1D*nbrSolPnts1D] = computeCorrectionFunction1D(solOrder, -solPntsLocalCoord[iSol][ZTA], factor);        
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
    
    // get the element shape and order
    const CFGeoShape::Type elemShape = frElemData->getShape();
    const CFPolyOrder::Type solOrder = frElemData->getPolyOrder();

    std::stringstream ss;
    ss << solOrder;
    std::string p=ss.str();
    // compute corr fct for correct element shape
    switch(elemShape)
    {
    case CFGeoShape::TRIAG:
      {	
      	CFuint nbrSolPnts = frElemData->getNbrOfSolPnts();
        const CFuint nbrFlxPnts = frElemData->getNbrOfFlxPnts();	
        std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
        cf_assert(corrfct.size() == nbrSolPnts);

        CFuint iFlx;
        RealMatrix phi;
        phi.resize(3,solOrder+1);
        phi=0.;
        for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
        {
          iFlx = 0;
          phi=computeTriagDivCorrFct(solOrder, m_cfactor, solPntsLocalCoord[iSol][0], solPntsLocalCoord[iSol][1]);
          for (CFuint f = 0; f < 3; ++f)
          {
            for (CFuint j = 0; j < (solOrder+1); ++j, ++iFlx)
            {
                corrfct[iSol][iFlx] = phi(f,j);
            }
          }
        }
        break;	
      }
    
      case CFGeoShape::QUAD:
      {
        CFuint iSol = 0;
        const CFuint nbrSolPnts1D = frElemData->getSolPntsLocalCoord1D()->size();
        const CFuint nbrSolPnts = frElemData->getNbrOfSolPnts();
        const CFuint nbrFlxPnts = frElemData->getNbrOfFlxPnts();
        std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
        cf_assert(corrfct.size() == nbrSolPnts);
	
	// loop over sol pnt coordinates in the standard domain
        for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
        {
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta, ++iSol)
          {
            cf_assert(corrfct[iSol].size() == nbrFlxPnts);
	    
	    // add the contriubutions of the flx pnts related to the current sol pnt
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
        std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
        cf_assert(corrfct.size() == nbrSolPnts);
	
	// loop over sol pnt coordinates in the standard domain
        for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
        {
          for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
          {
            for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++iSol)
            {
              cf_assert(corrfct[iSol].size() == nbrFlxPnts);
	      
	      // add the contriubutions of the flx pnts related to the current sol pnt
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
      case CFGeoShape::TETRA:
      {	
      	CFuint nbrSolPnts = frElemData->getNbrOfSolPnts();
        const CFuint nbrFlxPnts = frElemData->getNbrOfFlxPnts();	
        std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
        cf_assert(corrfct.size() == nbrSolPnts);

        CFuint iFlx;
        RealMatrix phi;
        phi.resize(4,(solOrder+1)*(solOrder+2)/2);
        phi=0.;
        for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
        {
          iFlx = 0;
          phi=computeTetraDivCorrFct(solOrder, m_cfactor, solPntsLocalCoord[iSol][0], solPntsLocalCoord[iSol][1], solPntsLocalCoord[iSol][2]);
          for (CFuint f = 0; f < 4; ++f)
          {
            for (CFuint j = 0; j < (solOrder+1)*(solOrder+2)/2; ++j, ++iFlx)
            {
              corrfct[iSol][iFlx] = phi(f,j);
            }
          }
        }
        break;	
      }
      case CFGeoShape::PRISM:
      {	
      	CFuint nbrSolPnts = frElemData->getNbrOfSolPnts();
        const CFuint nbrFlxPnts = frElemData->getNbrOfFlxPnts();	
        std::vector< RealVector > solPntsLocalCoord = *(frElemData->getSolPntsLocalCoords());
        cf_assert(corrfct.size() == nbrSolPnts);

        CFuint iFlx;
        RealMatrix phi;
        phi.resize(5,(solOrder+1)*(solOrder+1));
        phi=0.;
        for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
        {
          iFlx = 0;
          //phi=computeTetraDivCorrFct(solOrder, m_cfactor, solPntsLocalCoord[iSol][0], solPntsLocalCoord[iSol][1], solPntsLocalCoord[iSol][2]);
          for (CFuint f = 0; f < 5; ++f)
          {
            for (CFuint j = 0; j < (solOrder+1)*(solOrder+1); ++j, ++iFlx)
            {
              corrfct[iSol][iFlx] = 1.; //phi(f,j);
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
        case CFPolyOrder::ORDER0:
        {
            corrfct = -0.5*ksi+0.5;
        } break;
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
        case CFPolyOrder::ORDER0:
        {
            corrfct = -0.5;
        } break;
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
	case CFPolyOrder::ORDER6:
        {
            corrfct = 0.5*(1./16.*(1386.*pow(ksi,5.)-1260.*pow(ksi,3.)+210.*ksi)-(702364162.5*cfactor*0.125*(315.*pow(ksi,4.)-210.*pow(ksi,2.)+15.)+1./16.*(3003.*pow(ksi,6.)-3465.*pow(ksi,4.)+945.*pow(ksi,2.)-35.))/(1.+702364162.5*cfactor));
        } break;
	case CFPolyOrder::ORDER7:
        {
            corrfct = -0.5*(1./16.*(3003.*pow(ksi,6.)-3465.*pow(ksi,4.)+945.*pow(ksi,2.)-35.)-(136961011687.5*cfactor*1./16.*(1386.*pow(ksi,5.)-1260.*pow(ksi,3.)+210.*ksi)+1./128.*(51480.*pow(ksi,7.)-72072.*pow(ksi,5.)+27720.*pow(ksi,3.)-2520.*ksi))/(1.+136961011687.5*cfactor));
        } break;
	case CFPolyOrder::ORDER8:
        {
            corrfct = 0.5*(1./128.*(51480.*pow(ksi,7.)-72072.*pow(ksi,5.)+27720.*pow(ksi,3.)-2520.*ksi)-(34925057980312.5*cfactor*1./16.*(3003.*pow(ksi,6.)-3465.*pow(ksi,4.)+945.*pow(ksi,2.)-35.)+1./128.*(109395.*pow(ksi,8.)-180180.*pow(ksi,6.)+90090.*pow(ksi,4.)-13860.*pow(ksi,2.0)+315.))/(1.+34925057980312.5*cfactor));
        } break;
	case CFPolyOrder::ORDER9:
        {
            corrfct = -0.5*(1./128.*(109395.*pow(ksi,8.)-180180.*pow(ksi,6.)+90090.*pow(ksi,4.)-13860.*pow(ksi,2.0)+315.)-(11280793727640937.5*cfactor*1./128.*(51480.*pow(ksi,7.)-72072.*pow(ksi,5.)+27720.*pow(ksi,3.)-2520.*ksi)+1./256.*(461890.*pow(ksi,9.)-875160.*pow(ksi,7.)+540540.*pow(ksi,5.)-120120.*pow(ksi,3.)+6930.*ksi))/(1.+11280793727640937.5*cfactor));
        } break;
	case CFPolyOrder::ORDER10:
        {
            corrfct = 0.5*(1./256.*(461890.*pow(ksi,9.)-875160.*pow(ksi,7.)+540540.*pow(ksi,5.)-120120.*pow(ksi,3.)+6930.*ksi)-(4501036697328734062.5*cfactor*1./128.*(109395.*pow(ksi,8.)-180180.*pow(ksi,6.)+90090.*pow(ksi,4.)-13860.*pow(ksi,2.0)+315.)+1./256.*(969969.*pow(ksi,10.)-2078586.*pow(ksi,8.)+1531530.*pow(ksi,6.)-450450.*pow(ksi,4.)+45045*pow(ksi,2.)-693.))/(1.+4501036697328734062.5*cfactor));
        } break;
        default:
        {
            throw Common::NotImplementedException (FromHere(),"Derivative of VCJH1D Correction Function not implemented for order " + StringOps::to_str(solOrder) + ".");
        }
    }
    return corrfct;
}

//////////////////////////////////////////////////////////////////////////////
RealVector VCJH::computeIntRHSTriag(const CFPolyOrder::Type solOrder, CFuint f, CFuint j)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)/2; 

  RealVector flux_points, rhs;
  RealVector weights;
  RealVector G, ksi, eta;
  RealVector dubasis;
  flux_points.resize(solOrder+1);
  weights.resize(solOrder+1);
  dubasis.resize(nbrSolPnts);
  rhs.resize(nbrSolPnts);
  G.resize(3);
  ksi.resize(3);
  eta.resize(3);

  switch(solOrder)
    {
      case 0:
      {
        flux_points[0]=0.;
        weights[0]=2.;
        break;
      }
      case 1:
      {
        flux_points[0]=-1./sqrt(3.);
        weights[0]=1.;
        flux_points[1]=1./sqrt(3);
        weights[1]=1.;
        break;
      }
      case 2:
      {
        flux_points[0]=-sqrt(3./5.);
        weights[0]=5./9.;
        flux_points[1]=0.;
        weights[1]=8./9.;
        flux_points[2]=sqrt(3./5.);
        weights[2]=5./9.;
        break;
      }
      case 3:
      {
        flux_points[0]=-sqrt((3./7.)+2./7.*sqrt(6./5.));
        weights[0]=(18.-sqrt(30.))/36.;
        flux_points[1]=-sqrt((3./7.)-2./7.*sqrt(6./5.));
        weights[1]=(18.+sqrt(30.))/36.;
        flux_points[2]=sqrt((3./7.)-2./7.*sqrt(6./5.));
        weights[2]=(18.+sqrt(30.))/36.;
        flux_points[3]=sqrt((3./7.)+2./7.*sqrt(6./5.));
        weights[3]=(18.-sqrt(30.))/36.;
        break;
      }
      case 4:
      {
        flux_points[0]=-0.9061798459386640;
        weights[0]=0.2369268850561891;
        flux_points[1]=-0.5384693101056831;
        weights[1]=0.4786286704993665;
        flux_points[2]=0.00000000000000;
        weights[2]=0.5688888888888889;
        flux_points[3]=0.5384693101056831;
        weights[3]=0.4786286704993665;
        flux_points[4]=0.9061798459386640;
        weights[4]=0.2369268850561891;
        break;
      }
      case 5:
      {
        flux_points[0]=-0.9324695142031521;
        weights[0]=0.1713244923791704;
        flux_points[1]=-0.6612093864662645;
        weights[1]=0.3607615730481386;
        flux_points[2]=-0.2386191860831969;
        weights[2]=0.4679139345726910;
        flux_points[3]=0.2386191860831969;
        weights[3]=0.4679139345726910;
        flux_points[4]=0.6612093864662645;
        weights[4]=0.3607615730481386;
        flux_points[5]=0.9324695142031521;
        weights[5]=0.1713244923791704;
        break;
      }
      case 6:
      {
        flux_points[0]=-0.9491079123427585;
        weights[0]=0.1294849661688697;
        flux_points[1]=-0.7415311855993945;
        weights[1]=0.2797053914892766;
        flux_points[2]=-0.4058451513773972;
        weights[2]=0.3818300505051189;
        flux_points[3]=0.0000000000000000;
        weights[3]=0.4179591836734694;
        flux_points[4]=0.4058451513773972;
        weights[4]=0.3818300505051189;
        flux_points[5]=0.7415311855993945;
        weights[5]=0.2797053914892766;
        flux_points[6]=0.9491079123427585;
        weights[6]=0.1294849661688697;
        break;
      }
      case 7:
      {
        flux_points[0] = -0.960289856497536;
        flux_points[1] = -0.796666477413627;
        flux_points[2] = -0.525532409916329;
        flux_points[3] = -0.183434642495650;
        flux_points[4] = 0.183434642495650;
        flux_points[5] = 0.525532409916329;
        flux_points[6] = 0.796666477413627;
        flux_points[7] = 0.960289856497536;
        weights[0]= 0.1012285362903763 ;
        weights[1]= 0.2223810344533745 ;
        weights[2]= 0.3137066458778873 ;
        weights[3]= 0.3626837833783620 ;
        weights[4]= 0.3626837833783620 ;
        weights[5]= 0.3137066458778873 ;
        weights[6]= 0.2223810344533745 ;
        weights[7]= 0.1012285362903763 ;
      } break;
      case 8:
      {
        flux_points[0] = -0.968160239507626;
        flux_points[1] = -0.836031107326636;
        flux_points[2] = -0.613371432700590;
        flux_points[3] = -0.324253423403809;
        flux_points[4] = 0.0;
        flux_points[5] = 0.324253423403809;
        flux_points[6] = 0.613371432700590;
        flux_points[7] = 0.836031107326636;
        flux_points[8] = 0.968160239507626;
        weights[0]= 0.0812743883615744 ;
        weights[1]= 0.1806481606948574 ;
        weights[2]= 0.2606106964029354 ;
        weights[3]= 0.3123470770400029 ;
        weights[4]= 0.3302393550012598 ;
        weights[5]= 0.3123470770400029 ;
        weights[6]= 0.2606106964029354 ;
        weights[7]= 0.1806481606948574 ;
        weights[8]= 0.0812743883615744 ;
      } break;
      case 9:
      {
        flux_points[0] = -0.973906528517172;
        flux_points[1] = -0.865063366688985;
        flux_points[2] = -0.679409568299024;
        flux_points[3] = -0.433395394129247;
        flux_points[4] = -0.148874338981631;
        flux_points[5] = 0.148874338981631;
        flux_points[6] = 0.433395394129247;
        flux_points[7] = 0.679409568299024;
        flux_points[8] = 0.865063366688985;
        flux_points[9] = 0.973906528517172;
        weights[0]= 0.0666713443086881 ;
        weights[1]= 0.1494513491505806 ;
        weights[2]= 0.2190863625159820 ;
        weights[3]= 0.2692667193099963 ;
        weights[4]= 0.2955242247147529 ;
        weights[5]= 0.2955242247147529 ;
        weights[6]= 0.2692667193099963 ;
        weights[7]= 0.2190863625159820 ;
        weights[8]= 0.1494513491505806 ;
        weights[9]= 0.0666713443086881 ;
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"VCJH Correction Functions not implemented for this order "
                                               + StringOps::to_str(solOrder) + ".");
      }
    }

  G[0]=0.5;
  G[1]=sqrt(2.)/2.;
  G[2]=0.5;

  ksi[0]=(flux_points[j]+1.)/2.;
  ksi[1]=(-flux_points[j]+1.)/2.;
  ksi[2]=0.;

  eta[0]= 0.;
  eta[1]=(flux_points[j]+1.)/2.;
  eta[2]=(-flux_points[j]+1.)/2.;

  dubasis=computeDubiner2D(solOrder,ksi[f],eta[f]);

  for (CFuint i = 0; i < nbrSolPnts; ++i)            
    {
      rhs[i]  = dubasis[i]*weights[j]*G[f]; //computeDubiner2D(P,ksi,eta) in which we define a=f(ksi,eta) b=f(ksi,eta)
    }

  return rhs;
}

//////////////////////////////////////////////////////////////////////////////
RealVector VCJH::computeIntRHSTetra(const CFPolyOrder::Type solOrder, CFuint f, CFuint j)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)*(solOrder+3)/6; 

  std::vector< RealVector > flux_points;
  RealVector rhs, weights;
  RealVector G, ksi, eta, zta;
  RealVector dubasis;

  flux_points.resize((solOrder+1)*(solOrder+2)/2);
  for (CFuint i=0; i<(solOrder+1)*(solOrder+2)/2  ; ++i)
  {
    flux_points[i].resize(2);
  }

  weights.resize((solOrder+1)*(solOrder+2)/2);
  dubasis.resize(nbrSolPnts);
  rhs.resize(nbrSolPnts);
  G.resize(4);
  ksi.resize(4);
  eta.resize(4);
  zta.resize(4);


  switch(solOrder)
    {
      case 0:
      {
        flux_points[0][0]=0.333333333333333;
        flux_points[0][1]=0.333333333333333;
        weights[0]=1.0;
      } break;
      case 1:
      {
        flux_points[0][0]=0.1666666666667;
        flux_points[0][1]=0.1666666666667;
        weights[0]=0.333333333333333;
        flux_points[1][0]=0.6666666666667;
        flux_points[1][1]=0.1666666666667;
        weights[1]=0.333333333333333;
        flux_points[2][0]=0.1666666666667;
        flux_points[2][1]=0.6666666666667;
        weights[2]=0.333333333333333;
      } break;
      case 2:
      {
        flux_points[0][0]=0.091576213509780;
        flux_points[0][1]=0.091576213509780;
        weights[0]=0.109951743655333;
        flux_points[1][0]=0.445948490915964;
        flux_points[1][1]=0.108103018168071;
        weights[1]=0.223381589678000;
        flux_points[2][0]=0.816847572980440;
        flux_points[2][1]=0.091576213509780;
        weights[2]=0.109951743655333;
        flux_points[3][0]=0.108103018168071;
        flux_points[3][1]=0.445948490915964;
        weights[3]=0.223381589678000;
        flux_points[4][0]=0.445948490915964;
        flux_points[4][1]=0.445948490915964;
        weights[4]=0.223381589678000;
        flux_points[5][0]=0.091576213509780;
        flux_points[5][1]=0.816847572980440;
        weights[5]=0.109951743655333;
      } break;
      case 3:
      {
        flux_points[0][0]=0.055564052669793;
        flux_points[0][1]=0.055564052669793;

        flux_points[1][0]=0.295533711735893;
        flux_points[1][1]=0.070255540518384;

        flux_points[2][0]=0.634210747745723;
        flux_points[2][1]=0.070255540518384;

        flux_points[3][0]=0.888871894660413;
        flux_points[3][1]=0.055564052669793;

        flux_points[4][0]=0.070255540518384;
        flux_points[4][1]=0.295533711735893;

        flux_points[5][0]=0.333333333333333;
        flux_points[5][1]=0.333333333333333;

        flux_points[6][0]=0.634210747745723;
        flux_points[6][1]=0.295533711735893;

        flux_points[7][0]=0.070255540518384;
        flux_points[7][1]=0.634210747745723;

        flux_points[8][0]=0.295533711735893;
        flux_points[8][1]=0.634210747745723;

        flux_points[9][0]=0.055564052669793;
        flux_points[9][1]=0.888871894660413; 

        weights[0] = 0.041955512996649;
        weights[1]=0.112098412070887;
        weights[2]=0.112098412070887;
        weights[3]=0.041955512996649;
        weights[4]=0.112098412070887;
        weights[5]=0.201542988584730;
        weights[6]=0.112098412070887;
        weights[7]=0.112098412070887;
        weights[8]=0.112098412070887;
        weights[9]=0.041955512996649;

      } break;
      case 4:
      {
        flux_points[0][0] = 0.035870877695734; 
        flux_points[0][1] = 0.035870877695734;

        flux_points[1][0] = 0.201503881881800; 
        flux_points[1][1] =  0.047312487011716;

        flux_points[2][0] = 0.474308787777079;
        flux_points[2][1] = 0.051382424445843; 

        flux_points[3][0] = 0.751183631106484; 
        flux_points[3][1] = 0.047312487011716;

        flux_points[4][0] = 0.928258244608533; 
        flux_points[4][1] = 0.035870877695734;

        flux_points[5][0] = 0.047312487011716;
        flux_points[5][1] = 0.201503881881800;

        flux_points[6][0] = 0.241729395767967;
        flux_points[6][1] = 0.241729395767967;

        flux_points[7][0] = 0.516541208464066; 
        flux_points[7][1] = 0.241729395767967; 

        flux_points[8][0] = 0.751183631106484; 
        flux_points[8][1] = 0.201503881881800;

        flux_points[9][0] = 0.051382424445843;
        flux_points[9][1] = 0.474308787777079;

        flux_points[10][0] = 0.241729395767967; 
        flux_points[10][1] = 0.516541208464066; 

        flux_points[11][0] = 0.474308787777079; 
        flux_points[11][1] = 0.474308787777079;

        flux_points[12][0] = 0.047312487011716; 
        flux_points[12][1] = 0.751183631106484;

        flux_points[13][0] =  0.201503881881800;
        flux_points[13][1] =  0.751183631106484;

        flux_points[14][0] = 0.035870877695734; 
        flux_points[14][1] = 0.928258244608533; 
      
        weights[0] = 0.017915455012303;
        weights[1] = 0.055749810027115;
        weights[2] = 0.076206062385535;
        weights[3] = 0.055749810027115;
        weights[4] = 0.017915455012303;
        weights[5] = 0.055749810027115;
        weights[6] = 0.127712195881265;
        weights[7] = 0.127712195881265;
        weights[8] = 0.055749810027115;
        weights[9] = 0.076206062385535;
        weights[10] = 0.127712195881265;
        weights[11] = 0.076206062385535;
        weights[12] = 0.055749810027115;
        weights[13] = 0.055749810027115;
        weights[14] = 0.017915455012303;

      } break;
      case 5:
      {
        flux_points[0][0] = 0.028112952182664; 
        flux_points[0][1] = 0.028112952182664; 
        flux_points[1][0] = 0.148565812270887; 
        flux_points[1][1] = 0.033533207700614; 
        flux_points[2][0] = 0.357196298615681; 
        flux_points[2][1] = 0.037824789609186;
        flux_points[3][0] = 0.604978911775132; 
        flux_points[3][1] = 0.037824789609186;
        flux_points[4][0] = 0.817900980028499; 
        flux_points[4][1] = 0.033533207700614;
        flux_points[5][0] = 0.943774095634672;
        flux_points[5][1] = 0.028112952182664;
        flux_points[6][0] = 0.033533207700614;
        flux_points[6][1] = 0.148565812270887;
        flux_points[7][0] = 0.177139098469317;
        flux_points[7][1] = 0.177139098469317;
        flux_points[8][0] = 0.405508595867433;
        flux_points[8][1] = 0.188982808265134;
        flux_points[9][0] = 0.645721803061365; 
        flux_points[9][1] = 0.177139098469317;
        flux_points[10][0] = 0.817900980028499;
        flux_points[10][1] = 0.148565812270887;
        flux_points[11][0] = 0.037824789609186; 
        flux_points[11][1] = 0.357196298615681; 
        flux_points[12][0] = 0.188982808265134;
        flux_points[12][1] = 0.405508595867433;
        flux_points[13][0] = 0.405508595867433;
        flux_points[13][1] = 0.405508595867433;
        flux_points[14][0] = 0.604978911775132; 
        flux_points[14][1] = 0.357196298615681; 
        flux_points[15][0] = 0.037824789609186; 
        flux_points[15][1] = 0.604978911775132;
        flux_points[16][0] = 0.177139098469317;
        flux_points[16][1] = 0.645721803061365;
        flux_points[17][0] = 0.357196298615681;
        flux_points[17][1] = 0.604978911775132;
        flux_points[18][0] = 0.033533207700614; 
        flux_points[18][1] = 0.817900980028499;
        flux_points[19][0] = 0.148565812270887;
        flux_points[19][1] = 0.817900980028499;
        flux_points[20][0] = 0.028112952182664;
        flux_points[20][1] = 0.943774095634672;

        weights[0] = 0.010359374696538;
        weights[1] = 0.028969269372473;
        weights[2] = 0.046046366595935;
        weights[3] = 0.046046366595935;
        weights[4] = 0.028969269372473;
        weights[5] = 0.010359374696538;
        weights[6] = 0.028969269372473;
        weights[7] = 0.075394884326738;
        weights[8] = 0.097547802373242;
        weights[9] = 0.075394884326738;
        weights[10] = 0.028969269372473;
        weights[11] = 0.046046366595935;
        weights[12] = 0.097547802373242;
        weights[13] = 0.097547802373242;
        weights[14] = 0.046046366595935;
        weights[15] = 0.046046366595935;
        weights[16] = 0.075394884326738;
        weights[17] = 0.046046366595935;
        weights[18] = 0.028969269372473;
        weights[19] = 0.028969269372473;
        weights[20] = 0.010359374696538;
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"VCJH Correction Functions not implemented for this order "
                                               + StringOps::to_str(solOrder) + ".");
      }
    }

  //G[0]=0.5;
  //G[1]=0.5;
  //G[2]=0.5;
  //G[3]=sqrt(3.)/2.;

  G[0]=1.;
  G[1]=1.;
  G[2]=sqrt(3.);
  G[3]=1.;


/*Param =   [ x         y      0;                   
            x            0       y;
           0             x      y;      swapped
            1-x-y        x      y];     swapped   */

  ksi[0]= flux_points[j][KSI];
  ksi[1]= flux_points[j][ETA];
  ksi[2]= 1.-flux_points[j][KSI]-flux_points[j][ETA];
  ksi[3]= 0.;

  eta[0]= flux_points[j][ETA];
  eta[1]= 0.;
  eta[2]= flux_points[j][ETA];
  eta[3]= flux_points[j][KSI];

  zta[0]= 0.;
  zta[1]= flux_points[j][KSI];
  zta[2]= flux_points[j][KSI];
  zta[3]= flux_points[j][ETA];

  dubasis=computeDubiner3D(solOrder,ksi[f],eta[f],zta[f]);

  for (CFuint i = 0; i < nbrSolPnts; ++i)            
    {
      rhs[i]  = dubasis[i]*weights[j]*G[f]; //computeDubiner3D(P,ksi,eta,zta) in which we define a=f(ksi,eta,zta) b=f(ksi,eta,zta) c=f(ksi,eta,zta)
    }

  return rhs;
}

//////////////////////////////////////////////////////////////////////////////
RealVector VCJH::computeDubiner2D(const CFPolyOrder::Type solOrder,CFreal ksi, CFreal eta)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)/2; 

  RealVector mat_v, mat_w, Qav, Qbw, dubasis;
  CFuint k;
  CFreal A, B;
  Qav.resize(nbrSolPnts+1);
  Qbw.resize(nbrSolPnts+1);
  mat_v.resize(nbrSolPnts+1);
  mat_w.resize(nbrSolPnts+1);
  dubasis.resize(nbrSolPnts);
  //setup
  for (CFuint w = 0; w < solOrder+1; ++w)
  {
    for (CFuint v = 0; v < solOrder+1; ++v)
    {
      if ((w+v) <= solOrder)
      {
        k=int(w+(solOrder+1)*v+1-(v*(v-1))/2.);
        mat_v[k]=v;
        mat_w[k]=w;
    
      }
    }
  }

  A = ((2.*ksi)/(1.-eta))-1.;
  B = 2.*eta-1.; 

  //normalized nth order Jacobi polynomials
  for (CFuint i = 0; i < nbrSolPnts; ++i)
    {
      Qav[i] = ComputeJacobi(mat_v[i+1],  0.,  0. , A);
      Qbw[i] = ComputeJacobi(mat_w[i+1],  2.*mat_v[i+1] +1.,  0. , B);
    }
  //the 2D orthonormal dubinner formula

  for (CFuint i = 0; i < nbrSolPnts; ++i)
    {
      dubasis[i] = Qav[i]*pow((1-B), mat_v[i+1])*Qbw[i]*sqrt(2.)*0.25;
    }
  return dubasis;
}

//////////////////////////////////////////////////////////////////////////////
RealVector VCJH::computeDubiner3D(const CFPolyOrder::Type solOrder,CFreal ksi, CFreal eta, CFreal zta)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)*(solOrder+3)/6; 

  RealVector mat_u, mat_v, mat_w, Qau, Qbv, Qcw, dubasis;
  CFuint k;
  CFreal A, B, C;
  Qau.resize(nbrSolPnts+1);
  Qbv.resize(nbrSolPnts+1);
  Qcw.resize(nbrSolPnts+1);
  mat_u.resize(nbrSolPnts+1);
  mat_v.resize(nbrSolPnts+1);
  mat_w.resize(nbrSolPnts+1);
  dubasis.resize(nbrSolPnts);

  //setup
  for (CFuint w = 0; w < solOrder+1; ++w)
  {
    for (CFuint v = 0; v < solOrder+1; ++v)
    {
      for (CFuint u = 0; u < solOrder+1; ++u)
      {
        if ((w+v+u) <= solOrder)
        {
          k= round(1.+(11.+12.*solOrder+3.*pow(solOrder,2.))*u/6.+(2.*solOrder+3.)*v/2.+w-((2.+solOrder)*(pow(u,2.))/2.)-u*v-(pow(v,2.))/2.+(pow(u,3.))/6.);
          mat_u[k]=u;
          mat_v[k]=v;
          mat_w[k]=w;
        }
      }
    }
  }

  A = ((2.*ksi)/(1.-eta-zta))-1.;
  B = ((2.*eta)/(1.-zta))-1.;
  C = (2.*zta)-1.;

  //normalized nth order Jacobi polynomials
  for (CFuint i = 0; i < nbrSolPnts; ++i)
    {
      Qau[i] = ComputeJacobi(mat_u[i+1],  0.,  0. , A);
      Qbv[i] = ComputeJacobi(mat_v[i+1],  2.*mat_u[i+1] +1.,  0. , B);
      Qcw[i] = ComputeJacobi(mat_w[i+1],  2.*mat_u[i+1] +2.*mat_v[i+1] +2.,  0. , C);
    }
  //the 2D orthonormal dubinner formula

  for (CFuint i = 0; i < nbrSolPnts; ++i)
    {
      dubasis[i] = Qau[i]*Qbv[i]*pow((1-B), mat_u[i+1])*Qcw[i]*pow((1-C), mat_u[i+1]+mat_v[i+1])*sqrt(8.)*0.125;
    }
  return dubasis;
}

//////////////////////////////////////////////////////////////////////////////
CFreal VCJH::ComputeJacobi(CFuint N,CFuint alpha, CFuint beta,CFreal x)
{
  // Evaluate Jacobi Polynomial (normalized to be orthonormal) of type (alpha,beta) > -1
  // (alpha+beta <> -1) at points x for order N   

  RealVector PL;
  CFreal gamma0, gamma1, aold, H, anew, bnew;

  PL.resize(N+1);
  // Initial values P_0(x) and P_1(x)
  gamma0 = pow(2.,(alpha+beta+1.))/(alpha+beta+1.)*factorial(alpha)*factorial(beta)/factorial(alpha+beta);
  PL[0] = 1.0/sqrt(gamma0);
//  PL[N]=PL[0];
  if (N>=1) {
    gamma1 = (alpha+1.)*(beta+1.)/(alpha+beta+3.)*gamma0;
    PL[1] = ((alpha+beta+2.)*x/2. + (alpha-beta)/2.)/sqrt(gamma1);
  }
  if (N>1) {
    // Repeat value in recurrence.
    aold = 2./(2.+alpha+beta)*sqrt((alpha+1.)*(beta+1.)/(alpha+beta+3.));
    // Forward recurrence using the symmetry of the recurrence.
    for (CFuint i = 0; i < N-1; ++i){
        H = 2.*(i+1.)+alpha+beta;
        anew = 2./(H+2.)*sqrt( ((i+1.)+1.)*((i+1.)+1.+alpha+beta)*((i+1.)+1.+alpha)*((i+1.)+1.+beta)/(H+1.)/(H+3.));
        bnew = - (pow(alpha,2.)-pow(beta,2.))/H/(H+2);
        PL[i+2] = 1./anew*( -aold*PL[i] + (x-bnew)*PL[i+1]);
        aold =anew;
    }
  }
  return PL[N];
}
//////////////////////////////////////////////////////////////////////////////
RealVector VCJH::computeSigmasTriag(const CFPolyOrder::Type solOrder,CFuint f, CFuint j, CFreal ksi, CFreal eta, CFreal cfactor)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)/2; 
  RealMatrix L, LInv;
  RealVector R, sigmas;
  L.resize(nbrSolPnts,nbrSolPnts);
  LInv.resize(nbrSolPnts,nbrSolPnts);
  R.resize(nbrSolPnts);
  sigmas.resize(nbrSolPnts);
  LInv=0.;
  sigmas=0.;
  L=computeLhsTriag(solOrder, ksi,eta);

  //L = cfactor*lhs + I;

  for (CFuint i = 0; i < nbrSolPnts; ++i)
    {
      for (CFuint j = 0; j < nbrSolPnts; ++j)
      {
        L(i,j)= cfactor*L(i,j);
      }
      L(i,i) += 1.;
    }
  
  R=computeIntRHSTriag(solOrder, f, j);
  InvertMatrix(L,LInv);
  //sigmas=LInv*R;
  for (CFuint i = 0; i < nbrSolPnts; ++i)
  {
    sigmas[i]=0.;
    for (CFuint j = 0; j < nbrSolPnts; ++j)
    {
      sigmas[i]= sigmas[i] + LInv(i,j)*R[j];
    }
  }

  return sigmas;
}

//////////////////////////////////////////////////////////////////////////////
RealVector VCJH::computeSigmasTetra(const CFPolyOrder::Type solOrder,CFuint f, CFuint j, CFreal ksi, CFreal eta, CFreal zta, CFreal cfactor)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)*(solOrder+3)/6; 
  RealMatrix L, LInv;
  RealVector R, sigmas;
  L.resize(nbrSolPnts,nbrSolPnts);
  LInv.resize(nbrSolPnts,nbrSolPnts);
  R.resize(nbrSolPnts);
  sigmas.resize(nbrSolPnts);
  LInv=0.;
  sigmas=0.;
  L=computeLhsTetra(solOrder);

  //L = cfactor*lhs + I;

  for (CFuint i = 0; i < nbrSolPnts; ++i)
    {
      for (CFuint j = 0; j < nbrSolPnts; ++j)
      {
        L(i,j)= cfactor*L(i,j);
      }
      L(i,i) += 1.;
    }
  
  R=computeIntRHSTetra(solOrder, f, j);
  InvertMatrix(L,LInv);

  for (CFuint i = 0; i < nbrSolPnts; ++i)
  {
    sigmas[i]=0.;
    for (CFuint j = 0; j < nbrSolPnts; ++j)
    {
      sigmas[i]= sigmas[i] + LInv(i,j)*R[j];
    }
  }

  return sigmas;
}

//////////////////////////////////////////////////////////////////////////////
CFreal VCJH::factorial(CFreal n)
{
  return (n==1. || n==0.) ? 1. : factorial(n-1.)*n;
}
//////////////////////////////////////////////////////////////////////////////
RealMatrix VCJH::computeTriagDivCorrFct(const CFPolyOrder::Type solOrder, CFreal cfactor, CFreal ksi, CFreal eta)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)/2; 
  RealVector dubasis, sigmas;
  RealMatrix  phi;
  dubasis.resize(nbrSolPnts);
  sigmas.resize(nbrSolPnts);
  phi.resize(3,solOrder+1);
  dubasis=computeDubiner2D(solOrder,ksi,eta);
  for (CFuint f = 0; f < 3; ++f)
    {
			for (CFuint j = 0; j < solOrder+1; ++j)
      {
        sigmas=computeSigmasTriag(solOrder,f, j, ksi,eta, cfactor);
        dubasis=computeDubiner2D(solOrder,ksi,eta);

        for (CFuint k = 0; k < nbrSolPnts; ++k)
        {
          phi(f,j)=phi(f,j) + sigmas[k]*dubasis[k];
        }        
      }
    }
  return (phi*16.*4.);
}

//////////////////////////////////////////////////////////////////////////////
RealMatrix VCJH::computeTetraDivCorrFct(const CFPolyOrder::Type solOrder, CFreal cfactor, CFreal ksi, CFreal eta, CFreal zta)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)*(solOrder+3)/6; 
  RealVector dubasis, sigmas;
  RealMatrix  phi;
  dubasis.resize(nbrSolPnts);
  sigmas.resize(nbrSolPnts);
  phi.resize(4,(solOrder+1)*(solOrder+2)/2);
  dubasis=computeDubiner3D(solOrder,ksi,eta,zta);
  for (CFuint f = 0; f < 4; ++f)
    {
			for (CFuint j = 0; j < (solOrder+1)*(solOrder+2)/2; ++j)
      {
        sigmas=computeSigmasTetra(solOrder,f, j, ksi, eta, zta, cfactor);
        dubasis=computeDubiner3D(solOrder,ksi,eta,zta);

        for (CFuint k = 0; k < nbrSolPnts; ++k)
        {
          phi(f,j)=phi(f,j) + sigmas[k]*dubasis[k];
        }

      }
    }
  return (phi*16.*8.*2.);
}

//////////////////////////////////////////////////////////////////////////////
void VCJH::InvertMatrix(RealMatrix A, RealMatrix& AI)
{
  cf_assert(A.nbRows() == AI.nbRows());
  cf_assert(A.nbCols() == AI.nbCols());
  
  const CFuint n = A.nbRows();

  for (CFuint i = 0; i < n; ++i)
  {
    for (CFuint j = 0; j < n; ++j)
    {
      AI(i,j) = 0.;
    }
    AI(i,i) = 1.;
  }

  for (CFuint i = 0; i < n; ++i)
  {
    CFreal fac = fabs(A(i,i));
    CFuint k = i;
    for (CFuint j = i+1; j < n; ++j)
    {
      if (fabs(A(j,i)) > fac)
      {
        fac = fabs(A(j,i));
        k = j;
      }
    }

    if (fac < MathTools::MathConsts::CFrealEps()) throw MathTools::ZeroDeterminantException (FromHere(),"Matrix is singular to working precision!!!");

    SwapRows(A ,i,k);
    SwapRows(AI,i,k);

    fac = 1./A(i,i);
    A(i,i) = 1.;
    for (CFuint j = i+1; j < n; ++j)
    {
      A(i,j) = A(i,j)*fac;
    }
    for (CFuint j = 0; j < n; ++j)
    {
      AI(i,j) = AI(i,j)*fac;
    }
    for (CFuint k = i+1; k < n; ++k)
    {
      fac = A(k,i);
      for (CFuint j = i+1; j < n; ++j)
      {
        A(k,j) = A(k,j) - A(i,j)*fac;
      }
      for (CFuint j = 0; j < n; ++j)
      {
        AI(k,j) = AI(k,j) - AI(i,j)*fac;
      }
    }
  }

  for (CFuint i = 0; i < n-1; ++i)
  {
    const CFuint ii = n-2-i;
    for (CFuint j = 0; j < n; ++j)
    {
      for (CFuint k = ii+1; k < n; ++k)
      {
        AI(ii,j) = AI(ii,j) - AI(k,j)*A(ii,k);
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

void VCJH::SwapRows(RealMatrix& A, CFuint row1, CFuint row2)
{
  RealVector swapRow = A.getRow<RealVector>(row1);
  A.setRow(A.getRow<RealVector>(row2),row1);
  A.setRow(swapRow,row2);
}

//////////////////////////////////////////////////////////////////////////////
CFuint VCJH::round(double x)
{
    return (x > 0.0) ? (int)(x + 0.5) : (int)(x - 0.5);
}

//////////////////////////////////////////////////////////////////////////////
RealMatrix VCJH::computeLhsTriag(const CFPolyOrder::Type solOrder, CFreal ksi, CFreal eta)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)/2; 
  RealMatrix lhs;
  lhs.resize(nbrSolPnts,nbrSolPnts);
  lhs=0.;
  switch(solOrder)
    {
        case 0:
        {
            lhs(0,0) = 1./32;
        } break;
        case 1:
        {
          lhs(1,1) = 9./16.;
          lhs(1,2) = (3.*sqrt(3.))/16.;
          lhs(2,1) = (3.*sqrt(3.))/16.;
          lhs(2,2) = 15./16.;
        } break;
        case 2:
        {

          lhs(2,2) = 75./2.;
          lhs(2,4) = (75.*sqrt(3.))/4.;
          lhs(2,5) = (15.*sqrt(5.))/4.;
          lhs(4,2) = (75.*sqrt(3.))/4.;
          lhs(4,4) = 675/8.;
          lhs(4,5) = (105.*sqrt(15.))/8.;
          lhs(5,2) = (15.*sqrt(5.))/4.;
          lhs(5,4) = (105.*sqrt(15.))/8.;
          lhs(5,5) = 825./8.;


        } break;
        case 3:
        {
          lhs(3,3) = 11025./2.;
          lhs(3,6) = (6615.*sqrt(3.))/2.;
          lhs(3,8) = (2205.*sqrt(5.))/2.;
          lhs(3,9) = (315.*sqrt(7.))/2.;
          lhs(6,3) = (6615.*sqrt(3.))/2.;
          lhs(6,6) = 27783./2.;
          lhs(6,8) = (6615.*sqrt(15.))/2.;
          lhs(6,9) = (1701.*sqrt(21.))/2.;
          lhs(8,3) = (2205.*sqrt(5.))/2.;
          lhs(8,6) = (6615.*sqrt(15.))/2.;
          lhs(8,8) = 55125./2.;
          lhs(8,9) = (5355.*sqrt(35.))/2.;
          lhs(9,3) = (315.*sqrt(7.))/2.;
          lhs(9,6) = (1701.*sqrt(21.))/2.;
          lhs(9,8) = (5355.*sqrt(35.))/2.;
          lhs(9,9) = 47187./2.;
           
        } break;
        case 4:
        {
          lhs(4,4) = 1428840.;
          lhs(4,8) = 952560.*sqrt(3.);
          lhs(4,11) = 408240.*sqrt(5.);
          lhs(4,13) = 102060.*sqrt(7);
          lhs(4,14) = 34020.;

          lhs(8,4) = 952560.*sqrt(3.);
          lhs(8,8) = 3810240.;
          lhs(8,11) = 1088640.*sqrt(15.);
          lhs(8,13) = 476280.*sqrt(21.);
          lhs(8,14) = 249480.*sqrt(3.);

          lhs(11,4) = 408240.*sqrt(5.);
          lhs(11,8) = 1088640.*sqrt(15.);
          lhs(11,11) = 9331200.;
          lhs(11,13) = 1428840.*sqrt(35.);
          lhs(11,14) = 1176120.*sqrt(5.);

          lhs(13,4) = 102060.*sqrt(7.);
          lhs(13,8) = 476280.*sqrt(21.);
          lhs(13,11) = 1428840.*sqrt(35.);
          lhs(13,13) = 14645610.;
          lhs(13,14) = 2942730.*sqrt(7.);

          lhs(14,4) = 34020.;
          lhs(14,8) = 249480.*sqrt(3.);
          lhs(14,11) = 1176120.*sqrt(5.);
          lhs(14,13) = 2942730.*sqrt(7.);
          lhs(14,14) = 9113310.;

        } break;
        case 5:
        {
          lhs(5,5)=576298800.;
          lhs(5,10)=712984858.5293;
          lhs(5,14)=460229747.197;
          lhs(5,17)=181517060.1982;
          lhs(5,19)=41164200.;
          lhs(5,20)=4137157.7635;
          lhs(10,5)=712984858.5293;
          lhs(10,10)=1587762000;
          lhs(10,14)=1935914598.5851;
          lhs(10,17)=1302500907.2016;
          lhs(10,19)=458347409.0545;
          lhs(10,20)=66539269.1348;
          lhs(14,5)=460229747.197;
          lhs(14,10)=1935914598.5851;
          lhs(14,14)=4336942500.;
          lhs(14,17)=4841622079.1723;
          lhs(14,19)=2597010716.3258;
          lhs(14,20)=538538377.4382;
          lhs(17,5)=181517060.1982;
          lhs(17,10)=1302500907.2016;
          lhs(17,14)=4841622079.1723;
          lhs(17,17)=9136165500.;
          lhs(17,19)=7766337075.6233;
          lhs(17,20)=2375518871.6151;
          lhs(19,5)=41164200.;
          lhs(19,10)=458347409.0545;
          lhs(19,14)=2597010716.3258;
          lhs(19,17)=7766337075.6233;
          lhs(19,19)=11264289300.;
          lhs(19,20)=5517490900.1507;
          lhs(20,5)=4137157.7635;
          lhs(20,10)=66539269.1348;
          lhs(20,14)=538538377.4382;
          lhs(20,17)=2375518871.6151;
          lhs(20,19)=5517490900.1507;
          lhs(20,20)=5311399500.;

        } break;
        case 6:
        {
          lhs(6,6)=333923990400;
          lhs(6,12)=433779987929.2065;
          lhs(6,17)=311115309105.1613;
          lhs(6,21)=147246639232.7866;
          lhs(6,24)=45535089600;
          lhs(6,26)=8390155944.3564;
          lhs(6,27)=701620087.1181;
          lhs(12,6)=433779987929.2065;
          lhs(12,12)=939161223000;
          lhs(12,17)=1212451925360.93;
          lhs(12,21)=956394976481.067;
          lhs(12,24)=453497260107.8068;
          lhs(12,26)=119890455127.1727;
          lhs(12,27)=13671468433.1198;
          lhs(17,6)=311115309105.1613;
          lhs(17,12)=1212451925360.93;
          lhs(17,17)=2724727005000;
          lhs(17,21)=3429723903988.525;
          lhs(17,24)=2418214448044.663;
          lhs(17,26)=898962620350.6442;
          lhs(17,27)=137929839223.3467;
          lhs(21,6)=147246639232.7866;
          lhs(21,12)=956394976481.067;
          lhs(21,17)=3429723903988.525;
          lhs(21,21)=6817614804000;
          lhs(21,24)=7248550467686.724;
          lhs(21,26)=3851399166583.282;
          lhs(21,27)=805948908544.7104;
          lhs(24,6)=45535089600;
          lhs(24,12)=453497260107.8068;
          lhs(24,17)=2418214448044.663;
          lhs(24,21)=7248550467686.724;
          lhs(24,24)=11886728162400;
          lhs(24,26)=9405746184348.215;
          lhs(24,27)=2786165257768.314;
          lhs(26,6)=8390155944.3564;
          lhs(26,12)=119890455127.1727;
          lhs(26,17)=898962620350.6442;
          lhs(26,21)=3851399166583.282;
          lhs(26,24)=9405746184348.215;
          lhs(26,26)=11766393639000;
          lhs(26,27)=5321361946360.8;
          lhs(27,6)=701620087.1181;
          lhs(27,12)=13671468433.1198;
          lhs(27,17)=137929839223.3467;
          lhs(27,21)=805948908544.7104;
          lhs(27,24)=2786165257768.314;
          lhs(27,26)=5321361946360.8;
          lhs(27,27)=4347585333000;
        } break;
	case 7:
        {
          lhs(7,7)=262965142440000;
          lhs(7,14)=354253656808852;
          lhs(7,20)=274403702630752.3;
          lhs(7,25)=147581290685588.4;
          lhs(7,29)=55780484760000;
          lhs(7,32)=14230995274850.58;
          lhs(7,34)=2210103274422.171;
          lhs(7,35)=158268782797.6971;
          lhs(14,7)=354253656808852;
          lhs(14,14)=749937628440000;
          lhs(14,20)=1003371326676467;
          lhs(14,25)=880463011681661.1;
          lhs(14,29)=504543086970183.1;
          lhs(14,32)=183496673278558.5;
          lhs(14,34)=38705446186210.27;
          lhs(14,35)=3624602132651.74;
          lhs(20,7)=274403702630752.3;
          lhs(20,14)=1003371326676467;
          lhs(20,20)=2249812885320000;
          lhs(20,25)=3058016663135892;
          lhs(20,29)=2552785960837605;
          lhs(20,32)=1287708933532056;
          lhs(20,34)=362079555575615.3;
          lhs(20,35)=43765603288739.29;
          lhs(25,7)=147581290685588.4;
          lhs(25,14)=880463011681661.1;
          lhs(25,20)=3058016663135892;
          lhs(25,25)=6282910965240000;
          lhs(25,29)=7633977672736346;
          lhs(25,32)=5370497922119272;
          lhs(25,34)=2025498252783803;
          lhs(25,35)=317722025558606.7;
          lhs(29,7)=55780484760000;
          lhs(29,14)=504543086970183.1;
          lhs(29,20)=2552785960837605;
          lhs(29,25)=7633977672736346;
          lhs(29,29)=1.353437398404e+16;
          lhs(29,32)=1.354402632112825e+16;
          lhs(29,34)=7004487004908838;
          lhs(29,35)=1453708362044291;
          lhs(32,7)=14230995274850.58;
          lhs(32,14)=183496673278558.5;
          lhs(32,20)=1287708933532056;
          lhs(32,25)=5370497922119272;
          lhs(32,29)=1.354402632112825e+16;
          lhs(32,32)=1.9679949135e+16;
          lhs(32,34)=1.44991269093228e+16;
          lhs(32,35)=4132531463693933;
          lhs(34,7)=2210103274422.171;
          lhs(34,14)=38705446186210.27;
          lhs(34,20)=362079555575615.3;
          lhs(34,25)=2025498252783803;
          lhs(34,29)=7004487004908838;
          lhs(34,32)=1.44991269093228e+16;
          lhs(34,34)=1.59867835218e+16;
          lhs(34,35)=6700856371453843;
          lhs(35,7)=158268782797.6971;
          lhs(35,14)=3624602132651.74;
          lhs(35,20)=43765603288739.29;
          lhs(35,25)=317722025558606.7;
          lhs(35,29)=1453708362044291;
          lhs(35,32)=4132531463693933;
          lhs(35,34)=6700856371453843;
          lhs(35,35)=4753934047800000;
        } break;
        case 8:
        {
        } break;
        case 9:
        {
        } break;
        default:
        {
            throw Common::NotImplementedException (FromHere(),"computeLhsTriag Correction Function not implemented for order " + StringOps::to_str(solOrder) + ".");
        }
    }
  return lhs;
}

//////////////////////////////////////////////////////////////////////////////
RealMatrix VCJH::computeLhsTetra(const CFPolyOrder::Type solOrder)
{
  const CFuint nbrSolPnts = (solOrder+1)*(solOrder+2)*(solOrder+3)/6; 
  RealMatrix lhs;
  lhs.resize(nbrSolPnts,nbrSolPnts);
  lhs=0.;
  switch(solOrder)
    {
        case 0:
        {
          lhs(0,0)=0.011719;
        } break;
        case 1:
        {
          lhs(1,1)=0.3125;
          lhs(1,2)=0.11049;
          lhs(1,3)=0.19137;
          lhs(2,1)=0.11049;
          lhs(2,2)=0.39062;
          lhs(2,3)=0.27063;
          lhs(3,1)=0.19137;
          lhs(3,2)=0.27063;
          lhs(3,3)=0.70312;
        } break;
        case 2:
        {
          lhs(2,2)=24.6094;
          lhs(2,4)=13.9212;
          lhs(2,5)=2.8416;
          lhs(2,7)=24.1122;
          lhs(2,8)=4.9219;
          lhs(2,9)=6.3541;
          lhs(4,2)=13.9212;
          lhs(4,4)=43.3125;
          lhs(4,5)=20.8972;
          lhs(4,7)=34.0998;
          lhs(4,8)=27.8423;
          lhs(4,9)=14.3777;
          lhs(5,2)=2.8416;
          lhs(5,4)=20.8972;
          lhs(5,5)=43.6406;
          lhs(5,7)=13.9212;
          lhs(5,8)=42.6247;
          lhs(5,9)=13.9405;
          lhs(7,2)=24.1122;
          lhs(7,4)=34.0998;
          lhs(7,5)=13.9212;
          lhs(7,7)=82.6875;
          lhs(7,8)=28.9346;
          lhs(7,9)=49.8059;
          lhs(8,2)=4.9219;
          lhs(8,4)=27.8423;
          lhs(8,5)=42.6247;
          lhs(8,7)=28.9346;
          lhs(8,8)=94.5;
          lhs(8,9)=60.9995;
          lhs(9,2)=6.3541;
          lhs(9,4)=14.3777;
          lhs(9,5)=13.9405;
          lhs(9,7)=49.8059;
          lhs(9,8)=60.9995;
          lhs(9,9)=124.6875;
        } break;
        case 3:
        {
          lhs(3,3)=3969;
          lhs(3,6)=2806.5068;
          lhs(3,8)=982.0728;
          lhs(3,9)=141.75;
          lhs(3,12)=4861.0124;
          lhs(3,14)=1701;
          lhs(3,15)=245.5182;
          lhs(3,17)=2195.9816;
          lhs(3,18)=316.9626;
          lhs(3,19)=375.0352;
          lhs(6,3)=2806.5068;
          lhs(6,6)=7938;
          lhs(6,8)=6249.8731;
          lhs(6,9)=1603.7182;
          lhs(6,12)=6874.5097;
          lhs(6,14)=8419.5204;
          lhs(6,15)=2430.5062;
          lhs(6,17)=4658.3804;
          lhs(6,18)=2241.2643;
          lhs(6,19)=1060.7599;
          lhs(8,3)=982.0728;
          lhs(8,6)=6249.8731;
          lhs(8,8)=13527;
          lhs(8,9)=6699.1395;
          lhs(8,12)=4410.225;
          lhs(8,14)=14169.9077;
          lhs(8,15)=8849.25;
          lhs(8,17)=5252.5237;
          lhs(8,18)=5882.0935;
          lhs(8,19)=1763.146;
          lhs(9,3)=141.75;
          lhs(9,6)=1603.7182;
          lhs(9,8)=6699.1395;
          lhs(9,9)=10003.5;
          lhs(9,12)=1041.6455;
          lhs(9,14)=6439.5;
          lhs(9,15)=11363.9853;
          lhs(9,17)=2039.1257;
          lhs(9,18)=5841.1686;
          lhs(9,19)=1285.8351;
          lhs(12,3)=4861.0124;
          lhs(12,6)=6874.5097;
          lhs(12,8)=4410.225;
          lhs(12,9)=1041.6455;
          lhs(12,12)=15876;
          lhs(12,14)=9027.5944;
          lhs(12,15)=2004.6477;
          lhs(12,17)=15240.5972;
          lhs(12,18)=3105.5869;
          lhs(12,19)=4593.225;
          lhs(14,3)=1701;
          lhs(14,6)=8419.5204;
          lhs(14,8)=14169.9077;
          lhs(14,9)=6439.5;
          lhs(14,12)=9027.5944;
          lhs(14,14)=28674;
          lhs(14,15)=12977.3907;
          lhs(14,17)=19450.1224;
          lhs(14,18)=17387.6646;
          lhs(14,19)=9643.7635;
          lhs(15,3)=245.5182;
          lhs(15,6)=2430.5062;
          lhs(15,8)=8849.25;
          lhs(15,9)=11363.9853;
          lhs(15,12)=2004.6477;
          lhs(15,14)=12977.3907;
          lhs(15,15)=23490;
          lhs(15,17)=6610.935;
          lhs(15,18)=22273.5272;
          lhs(15,19)=7887.7584;
          lhs(17,3)=2195.9816;
          lhs(17,6)=4658.3804;
          lhs(17,8)=5252.5237;
          lhs(17,9)=2039.1257;
          lhs(17,12)=15240.5972;
          lhs(17,14)=19450.1224;
          lhs(17,15)=6610.935;
          lhs(17,17)=39690;
          lhs(17,18)=14731.0921;
          lhs(17,19)=23240.0843;
          lhs(18,3)=316.9626;
          lhs(18,6)=2241.2643;
          lhs(18,8)=5882.0935;
          lhs(18,9)=5841.1686;
          lhs(18,12)=3105.5869;
          lhs(18,14)=17387.6646;
          lhs(18,15)=22273.5272;
          lhs(18,17)=14731.0921;
          lhs(18,18)=43942.5;
          lhs(18,19)=26835.3379;
          lhs(19,3)=375.0352;
          lhs(19,6)=1060.7599;
          lhs(19,8)=1763.146;
          lhs(19,9)=1285.8351;
          lhs(19,12)=4593.225;
          lhs(19,14)=9643.7635;
          lhs(19,15)=7887.7584;
          lhs(19,17)=23240.0843;
          lhs(19,18)=26835.3379;
          lhs(19,19)=42525;
        } break;
        case 4:
        {
          lhs(4,4)=1091475;
          lhs(4,8)=882044.9989;
          lhs(4,11)=405105.0333;
          lhs(4,13)=103950;
          lhs(4,14)=11621.9633;
          lhs(4,18)=1527746.7526;
          lhs(4,21)=701662.5;
          lhs(4,23)=180046.6814;
          lhs(4,24)=20129.8309;
          lhs(4,27)=905842.3924;
          lhs(4,29)=232439.2663;
          lhs(4,30)=25987.5;
          lhs(4,32)=275025.8488;
          lhs(4,33)=30748.8247;
          lhs(4,34)=34865.8899;
          lhs(8,4)=882044.9989;
          lhs(8,8)=2316600;
          lhs(8,11)=2291620.1289;
          lhs(8,13)=1029052.4987;
          lhs(8,14)=178447.3284;
          lhs(8,18)=2160560.1774;
          lhs(8,21)=3118659.1031;
          lhs(8,23)=1564121.6753;
          lhs(8,24)=284678.7993;
          lhs(8,27)=1830077.9953;
          lhs(8,29)=1455754.5209;
          lhs(8,30)=304515.5353;
          lhs(8,32)=722326.9585;
          lhs(8,33)=248488.0279;
          lhs(8,34)=112703.5758;
          lhs(11,4)=405105.0333;
          lhs(11,8)=2291620.1289;
          lhs(11,11)=5061993.75;
          lhs(11,13)=4089631.7643;
          lhs(11,14)=1117205.6173;
          lhs(11,18)=1701086.7835;
          lhs(11,21)=5555726.1704;
          lhs(11,23)=5479650;
          lhs(11,24)=1643677.6686;
          lhs(11,27)=2241378.6389;
          lhs(11,29)=3795910.9776;
          lhs(11,30)=1466094.4061;
          lhs(11,32)=1224922.4833;
          lhs(11,33)=855940.6062;
          lhs(11,34)=245871.5065;
          lhs(13,4)=103950;
          lhs(13,8)=1029052.4987;
          lhs(13,11)=4089631.7643;
          lhs(13,13)=6939900;
          lhs(13,14)=3316133.532;
          lhs(13,18)=691123.5309;
          lhs(13,21)=4076325;
          lhs(13,23)=8247852.7406;
          lhs(13,24)=4487993.7367;
          lhs(13,27)=1380331.2646;
          lhs(13,29)=4615579.7158;
          lhs(13,30)=3356100;
          lhs(13,32)=1139392.8021;
          lhs(13,33)=1511085.0982;
          lhs(13,34)=318773.8509;
          lhs(14,4)=11621.9633;
          lhs(14,8)=178447.3284;
          lhs(14,11)=1117205.6173;
          lhs(14,13)=3316133.532;
          lhs(14,14)=3862361.25;
          lhs(14,18)=113871.5197;
          lhs(14,21)=1060919.2224;
          lhs(14,23)=3734562.9214;
          lhs(14,24)=4786240.9485;
          lhs(14,27)=327942.1698;
          lhs(14,29)=1885950;
          lhs(14,30)=3049935.2294;
          lhs(14,32)=404127.41;
          lhs(14,33)=1114836.9228;
          lhs(14,34)=185996.25;
          lhs(18,4)=1527746.7526;
          lhs(18,8)=2160560.1774;
          lhs(18,11)=1701086.7835;
          lhs(18,13)=691123.5309;
          lhs(18,14)=113871.5197;
          lhs(18,18)=4811400;
          lhs(18,21)=3437430.1933;
          lhs(18,23)=1323067.4983;
          lhs(18,24)=211319.2046;
          lhs(18,27)=5705618.5253;
          lhs(18,29)=2033419.9947;
          lhs(18,30)=309186.8428;
          lhs(18,32)=2983408.9805;
          lhs(18,33)=430393.8894;
          lhs(18,34)=585624.9585;
          lhs(21,4)=701662.5;
          lhs(21,8)=3118659.1031;
          lhs(21,11)=5555726.1704;
          lhs(21,13)=4076325;
          lhs(21,14)=1060919.2224;
          lhs(21,18)=3437430.1933;
          lhs(21,21)=10600115.625;
          lhs(21,23)=7986356.3699;
          lhs(21,24)=2038145.3829;
          lhs(21,27)=7473199.7372;
          lhs(21,29)=10459766.9817;
          lhs(21,30)=2898534.375;
          lhs(21,32)=6011279.2663;
          lhs(21,33)=3162736.252;
          lhs(21,34)=1613792.62;
          lhs(23,4)=180046.6814;
          lhs(23,8)=1564121.6753;
          lhs(23,11)=5479650;
          lhs(23,13)=8247852.7406;
          lhs(23,14)=3734562.9214;
          lhs(23,18)=1323067.4983;
          lhs(23,21)=7986356.3699;
          lhs(23,23)=15770700;
          lhs(23,24)=7321836.8872;
          lhs(23,27)=4183906.7927;
          lhs(23,29)=15164472.6429;
          lhs(23,30)=9431016.6472;
          lhs(23,32)=5829723.6703;
          lhs(23,33)=7202581.6934;
          lhs(23,34)=2346563.1498;
          lhs(24,4)=20129.8309;
          lhs(24,8)=284678.7993;
          lhs(24,11)=1643677.6686;
          lhs(24,13)=4487993.7367;
          lhs(24,14)=4786240.9485;
          lhs(24,18)=211319.2046;
          lhs(24,21)=2038145.3829;
          lhs(24,23)=7321836.8872;
          lhs(24,24)=9514209.375;
          lhs(24,27)=893784.375;
          lhs(24,29)=5855803.9728;
          lhs(24,30)=10461521.0687;
          lhs(24,32)=1876729.0328;
          lhs(24,33)=6079244.9169;
          lhs(24,34)=1466094.4061;
          lhs(27,4)=905842.3924;
          lhs(27,8)=1830077.9953;
          lhs(27,11)=2241378.6389;
          lhs(27,13)=1380331.2646;
          lhs(27,14)=327942.1698;
          lhs(27,18)=5705618.5253;
          lhs(27,21)=7473199.7372;
          lhs(27,23)=4183906.7927;
          lhs(27,24)=893784.375;
          lhs(27,27)=15160921.875;
          lhs(27,29)=9066636.4586;
          lhs(27,30)=1887171.6508;
          lhs(27,32)=14379802.1836;
          lhs(27,33)=3215421.5186;
          lhs(27,34)=4456155.3658;
          lhs(29,4)=232439.2663;
          lhs(29,8)=1455754.5209;
          lhs(29,11)=3795910.9776;
          lhs(29,13)=4615579.7158;
          lhs(29,14)=1885950;
          lhs(29,18)=2033419.9947;
          lhs(29,21)=10459766.9817;
          lhs(29,23)=15164472.6429;
          lhs(29,24)=5855803.9728;
          lhs(29,27)=9066636.4586;
          lhs(29,29)=28314000;
          lhs(29,30)=12280541.2341;
          lhs(29,32)=18039310.4746;
          lhs(29,33)=16868252.0588;
          lhs(29,34)=9147600;
          lhs(30,4)=25987.5;
          lhs(30,8)=304515.5353;
          lhs(30,11)=1466094.4061;
          lhs(30,13)=3356100;
          lhs(30,14)=3049935.2294;
          lhs(30,18)=309186.8428;
          lhs(30,21)=2898534.375;
          lhs(30,23)=9431016.6472;
          lhs(30,24)=10461521.0687;
          lhs(30,27)=1887171.6508;
          lhs(30,29)=12280541.2341;
          lhs(30,30)=20610253.125;
          lhs(30,32)=5729705.183;
          lhs(30,33)=19320511.5027;
          lhs(30,34)=6775604.6115;
          lhs(32,4)=275025.8488;
          lhs(32,8)=722326.9585;
          lhs(32,11)=1224922.4833;
          lhs(32,13)=1139392.8021;
          lhs(32,14)=404127.41;
          lhs(32,18)=2983408.9805;
          lhs(32,21)=6011279.2663;
          lhs(32,23)=5829723.6703;
          lhs(32,24)=1876729.0328;
          lhs(32,27)=14379802.1836;
          lhs(32,29)=18039310.4746;
          lhs(32,30)=5729705.183;
          lhs(32,32)=29521800;
          lhs(32,33)=12086841.8456;
          lhs(32,34)=16235379.4272;
          lhs(33,4)=30748.8247;
          lhs(33,8)=248488.0279;
          lhs(33,11)=855940.6062;
          lhs(33,13)=1511085.0982;
          lhs(33,14)=1114836.9228;
          lhs(33,18)=430393.8894;
          lhs(33,21)=3162736.252;
          lhs(33,23)=7202581.6934;
          lhs(33,24)=6079244.9169;
          lhs(33,27)=3215421.5186;
          lhs(33,29)=16868252.0588;
          lhs(33,30)=19320511.5027;
          lhs(33,32)=12086841.8456;
          lhs(33,33)=32224500;
          lhs(33,34)=18151706.0198;
          lhs(34,4)=34865.8899;
          lhs(34,8)=112703.5758;
          lhs(34,11)=245871.5065;
          lhs(34,13)=318773.8509;
          lhs(34,14)=185996.25;
          lhs(34,18)=585624.9585;
          lhs(34,21)=1613792.62;
          lhs(34,23)=2346563.1498;
          lhs(34,24)=1466094.4061;
          lhs(34,27)=4456155.3658;
          lhs(34,29)=9147600;
          lhs(34,30)=6775604.6115;
          lhs(34,32)=16235379.4272;
          lhs(34,33)=18151706.0198;
          lhs(34,34)=23295195;
        } break;
        case 5:
        {
          lhs(5,5)=458686800;
          lhs(5,10)=405425683.4009;
          lhs(5,14)=220685789.5448;
          lhs(5,17)=76447800;
          lhs(5,19)=15540225.23;
          lhs(5,20)=1418621.9845;
          lhs(5,25)=702217882.3438;
          lhs(5,29)=382239000;
          lhs(5,32)=132411473.7269;
          lhs(5,34)=26916459.6595;
          lhs(5,35)=2457125.3539;
          lhs(5,39)=493468427.0903;
          lhs(5,42)=170942477.5303;
          lhs(5,44)=34749000;
          lhs(5,45)=3172135.1918;
          lhs(5,48)=202261867.078;
          lhs(5,50)=41115571.2766;
          lhs(5,51)=3753320.9755;
          lhs(5,53)=46620675.6901;
          lhs(5,54)=4255865.9536;
          lhs(5,55)=4705036.842;
          lhs(10,5)=405425683.4009;
          lhs(10,10)=1003377375;
          lhs(10,14)=1131351032.665;
          lhs(10,17)=675709472.3349;
          lhs(10,19)=211530523.8454;
          lhs(10,20)=27585723.6931;
          lhs(10,25)=993086052.9515;
          lhs(10,29)=1554131786.3703;
          lhs(10,32)=1029919560.7709;
          lhs(10,34)=337832397.9258;
          lhs(10,35)=45173700;
          lhs(10,39)=959570895.5168;
          lhs(10,42)=966996680.4361;
          lhs(10,44)=362425989.7069;
          lhs(10,45)=51589881.014;
          lhs(10,48)=500572582.2521;
          lhs(10,50)=297999267.4358;
          lhs(10,51)=49098976.2041;
          lhs(10,53)=140104632.6768;
          lhs(10,54)=37616895.9451;
          lhs(10,55)=16634817.2837;
          lhs(14,5)=220685789.5448;
          lhs(14,10)=1131351032.665;
          lhs(14,14)=2527024500;
          lhs(14,17)=2574667544.689;
          lhs(14,19)=1236661785.4658;
          lhs(14,20)=230014234.5221;
          lhs(14,25)=878422314.0354;
          lhs(14,29)=2868915264.0821;
          lhs(14,32)=3491116200;
          lhs(14,34)=1828566502.0666;
          lhs(14,35)=356783429.11;
          lhs(14,39)=1234584559.9411;
          lhs(14,42)=2500240030.591;
          lhs(14,44)=1651799697.5018;
          lhs(14,45)=362318516.3124;
          lhs(14,48)=856356918.2981;
          lhs(14,50)=993043270.375;
          lhs(14,51)=281346350.3423;
          lhs(14,53)=300567132.8641;
          lhs(14,54)=153570334.6216;
          lhs(14,55)=43010526.2156;
          lhs(17,5)=76447800;
          lhs(17,10)=675709472.3349;
          lhs(17,14)=2574667544.689;
          lhs(17,17)=4727022300;
          lhs(17,19)=3688213454.5933;
          lhs(17,20)=1018097710.8998;
          lhs(17,25)=468145254.8958;
          lhs(17,29)=2643819750;
          lhs(17,32)=5781967686.0731;
          lhs(17,34)=5058051377.6758;
          lhs(17,35)=1495570298.7679;
          lhs(17,39)=945814485.2565;
          lhs(17,42)=3418849550.6062;
          lhs(17,44)=3900575250;
          lhs(17,45)=1353444348.4976;
          lhs(17,48)=910178401.851;
          lhs(17,50)=1850200707.4464;
          lhs(17,51)=870770466.312;
          lhs(17,53)=419586081.2108;
          lhs(17,54)=366004472.0095;
          lhs(17,55)=75280589.4727;
          lhs(19,5)=15540225.23;
          lhs(19,10)=211530523.8454;
          lhs(19,14)=1236661785.4658;
          lhs(19,17)=3688213454.5933;
          lhs(19,19)=5210559900;
          lhs(19,20)=2353772374.9635;
          lhs(19,25)=137987880.8429;
          lhs(19,29)=1201777417.7888;
          lhs(19,32)=4259978348.7725;
          lhs(19,34)=6665805153.9524;
          lhs(19,35)=3270489995.3185;
          lhs(19,39)=387871993.7454;
          lhs(19,42)=2258685000;
          lhs(19,44)=4515141803.1972;
          lhs(19,45)=2643902987.9919;
          lhs(19,48)=518056198.085;
          lhs(19,50)=1797790397.0446;
          lhs(19,51)=1439802547.0897;
          lhs(19,53)=332642700;
          lhs(19,54)=490960014.9113;
          lhs(19,55)=79862301.3628;
          lhs(20,5)=1418621.9845;
          lhs(20,10)=27585723.6931;
          lhs(20,14)=230014234.5221;
          lhs(20,17)=1018097710.8998;
          lhs(20,19)=2353772374.9635;
          lhs(20,20)=2248988625;
          lhs(20,25)=17374500;
          lhs(20,29)=216339852.6413;
          lhs(20,32)=1138468080.6612;
          lhs(20,34)=2900906794.3488;
          lhs(20,35)=2950076566.8485;
          lhs(20,39)=65626352.9968;
          lhs(20,42)=566754820.9334;
          lhs(20,44)=1828195346.8855;
          lhs(20,45)=2140822997.4084;
          lhs(20,48)=118855164.2236;
          lhs(20,50)=657045180.5073;
          lhs(20,51)=1006549959.5003;
          lhs(20,53)=107131656.7046;
          lhs(20,54)=288495675;
          lhs(20,55)=38823912.3022;
          lhs(25,5)=702217882.3438;
          lhs(25,10)=993086052.9515;
          lhs(25,14)=878422314.0354;
          lhs(25,17)=468145254.8958;
          lhs(25,19)=137987880.8429;
          lhs(25,20)=17374500;
          lhs(25,25)=2150094375;
          lhs(25,29)=1755544705.8594;
          lhs(25,32)=891936503.4821;
          lhs(25,34)=255484918.4107;
          lhs(25,35)=31598192.5939;
          lhs(25,39)=2870771395.0446;
          lhs(25,42)=1360845997.2784;
          lhs(25,44)=372388270.9399;
          lhs(25,45)=44678147.5363;
          lhs(25,48)=1981753475.055;
          lhs(25,50)=516149871.8172;
          lhs(25,51)=59759188.0003;
          lhs(25,53)=699455809.79;
          lhs(25,54)=78185250;
          lhs(25,55)=100843220.4851;
          lhs(29,5)=382239000;
          lhs(29,10)=1554131786.3703;
          lhs(29,14)=2868915264.0821;
          lhs(29,17)=2643819750;
          lhs(29,19)=1201777417.7888;
          lhs(29,20)=216339852.6413;
          lhs(29,25)=1755544705.8594;
          lhs(29,29)=5255786250;
          lhs(29,32)=5064738870.0525;
          lhs(29,34)=2281169956.1403;
          lhs(29,35)=404197120.724;
          lhs(29,39)=3824380309.95;
          lhs(29,42)=6510059352.6126;
          lhs(29,44)=3199803750;
          lhs(29,45)=574156469.7142;
          lhs(29,48)=3859830630.0719;
          lhs(29,50)=3467413177.6587;
          lhs(29,51)=701871022.4153;
          lhs(29,53)=1818206351.9133;
          lhs(29,54)=680938552.5758;
          lhs(29,55)=329352578.9431;
          lhs(32,5)=132411473.7269;
          lhs(32,10)=1029919560.7709;
          lhs(32,14)=3491116200;
          lhs(32,17)=5781967686.0731;
          lhs(32,19)=4259978348.7725;
          lhs(32,20)=1138468080.6612;
          lhs(32,25)=891936503.4821;
          lhs(32,29)=5064738870.0525;
          lhs(32,32)=10557441180;
          lhs(32,34)=8111220558.8131;
          lhs(32,35)=2162831077.6188;
          lhs(32,39)=2720834434.0241;
          lhs(32,42)=10234535177.8534;
          lhs(32,44)=10250854924.3547;
          lhs(32,45)=2965456486.6706;
          lhs(32,48)=4157223585.1924;
          lhs(32,50)=7980744561.7547;
          lhs(32,51)=3097915859.6865;
          lhs(32,53)=2740095593.3343;
          lhs(32,54)=2093470801.5613;
          lhs(32,55)=646516120.3778;
          lhs(34,5)=26916459.6595;
          lhs(34,10)=337832397.9258;
          lhs(34,14)=1828566502.0666;
          lhs(34,17)=5058051377.6758;
          lhs(34,19)=6665805153.9524;
          lhs(34,20)=2900906794.3488;
          lhs(34,25)=255484918.4107;
          lhs(34,29)=2281169956.1403;
          lhs(34,32)=8111220558.8131;
          lhs(34,34)=12438404550;
          lhs(34,35)=5575662026.6113;
          lhs(34,39)=1022199750;
          lhs(34,42)=6561389770.0562;
          lhs(34,44)=13593016040.6067;
          lhs(34,45)=7188134077.1063;
          lhs(34,48)=2194586064.0518;
          lhs(34,50)=8339346503.7719;
          lhs(34,51)=6246399426.5981;
          lhs(34,53)=2198468187.746;
          lhs(34,54)=3212666657.3166;
          lhs(34,55)=734423151.8869;
          lhs(35,5)=2457125.3539;
          lhs(35,10)=45173700;
          lhs(35,14)=356783429.11;
          lhs(35,17)=1495570298.7679;
          lhs(35,19)=3270489995.3185;
          lhs(35,20)=2950076566.8485;
          lhs(35,25)=31598192.5939;
          lhs(35,29)=404197120.724;
          lhs(35,32)=2162831077.6188;
          lhs(35,34)=5575662026.6113;
          lhs(35,35)=5719859145;
          lhs(35,39)=160192827.1855;
          lhs(35,42)=1538770096.4061;
          lhs(35,44)=5394023047.8341;
          lhs(35,45)=6770855818.5976;
          lhs(35,48)=454199191.5153;
          lhs(35,50)=2890879220.1052;
          lhs(35,51)=4953471453.8346;
          lhs(35,53)=659465763.0061;
          lhs(35,54)=2053745729.8917;
          lhs(35,55)=373678486.4369;
          lhs(39,5)=493468427.0903;
          lhs(39,10)=959570895.5168;
          lhs(39,14)=1234584559.9411;
          lhs(39,17)=945814485.2565;
          lhs(39,19)=387871993.7454;
          lhs(39,20)=65626352.9968;
          lhs(39,25)=2870771395.0446;
          lhs(39,29)=3824380309.95;
          lhs(39,32)=2720834434.0241;
          lhs(39,34)=1022199750;
          lhs(39,35)=160192827.1855;
          lhs(39,39)=7697868750;
          lhs(39,42)=5756221010.6262;
          lhs(39,44)=2112194403.8336;
          lhs(39,45)=316696156.7307;
          lhs(39,48)=9160923077.9015;
          lhs(39,50)=3574051601.1613;
          lhs(39,51)=539467508.5288;
          lhs(39,53)=4975461437.0094;
          lhs(39,54)=796673401.4956;
          lhs(39,55)=1022486368.1032;
          lhs(42,5)=170942477.5303;
          lhs(42,10)=966996680.4361;
          lhs(42,14)=2500240030.591;
          lhs(42,17)=3418849550.6062;
          lhs(42,19)=2258685000;
          lhs(42,20)=566754820.9334;
          lhs(42,25)=1360845997.2784;
          lhs(42,29)=6510059352.6126;
          lhs(42,32)=10234535177.8534;
          lhs(42,34)=6561389770.0562;
          lhs(42,35)=1538770096.4061;
          lhs(42,39)=5756221010.6262;
          lhs(42,42)=18360213300;
          lhs(42,44)=13277827440.2898;
          lhs(42,45)=3148395057.6733;
          lhs(42,48)=12090718993.4013;
          lhs(42,50)=18013564464.9167;
          lhs(42,51)=4949443648.8831;
          lhs(42,53)=10042461000;
          lhs(42,54)=5627367830.2375;
          lhs(42,55)=2833597343.3621;
          lhs(44,5)=34749000;
          lhs(44,10)=362425989.7069;
          lhs(44,14)=1651799697.5018;
          lhs(44,17)=3900575250;
          lhs(44,19)=4515141803.1972;
          lhs(44,20)=1828195346.8855;
          lhs(44,25)=372388270.9399;
          lhs(44,29)=3199803750;
          lhs(44,32)=10250854924.3547;
          lhs(44,34)=13593016040.6067;
          lhs(44,35)=5394023047.8341;
          lhs(44,39)=2112194403.8336;
          lhs(44,42)=13277827440.2898;
          lhs(44,44)=25131161250;
          lhs(44,45)=10957516205.5351;
          lhs(44,48)=6388103968.5469;
          lhs(44,50)=24013985478.3307;
          lhs(44,51)=14968471524.2244;
          lhs(44,53)=8996377660.438;
          lhs(44,54)=11712014138.6685;
          lhs(44,55)=3816070032.6464;
          lhs(45,5)=3172135.1918;
          lhs(45,10)=51589881.014;
          lhs(45,14)=362318516.3124;
          lhs(45,17)=1353444348.4976;
          lhs(45,19)=2643902987.9919;
          lhs(45,20)=2140822997.4084;
          lhs(45,25)=44678147.5363;
          lhs(45,29)=574156469.7142;
          lhs(45,32)=2965456486.6706;
          lhs(45,34)=7188134077.1063;
          lhs(45,35)=6770855818.5976;
          lhs(45,39)=316696156.7307;
          lhs(45,42)=3148395057.6733;
          lhs(45,44)=10957516205.5351;
          lhs(45,45)=13322275200;
          lhs(45,48)=1267854319.2821;
          lhs(45,50)=8510257232.4101;
          lhs(45,51)=14983677037.83;
          lhs(45,53)=2625740327.7354;
          lhs(45,54)=8859341128.8642;
          lhs(45,55)=2155875551.0414;
          lhs(48,5)=202261867.078;
          lhs(48,10)=500572582.2521;
          lhs(48,14)=856356918.2981;
          lhs(48,17)=910178401.851;
          lhs(48,19)=518056198.085;
          lhs(48,20)=118855164.2236;
          lhs(48,25)=1981753475.055;
          lhs(48,29)=3859830630.0719;
          lhs(48,32)=4157223585.1924;
          lhs(48,34)=2194586064.0518;
          lhs(48,35)=454199191.5153;
          lhs(48,39)=9160923077.9015;
          lhs(48,42)=12090718993.4013;
          lhs(48,44)=6388103968.5469;
          lhs(48,45)=1267854319.2821;
          lhs(48,48)=20174574420;
          lhs(48,50)=13206083400.4778;
          lhs(48,51)=2712972683.2207;
          lhs(48,53)=18082628247.4424;
          lhs(48,54)=4549025022.2853;
          lhs(48,55)=5551967397.4776;
          lhs(50,5)=41115571.2766;
          lhs(50,10)=297999267.4358;
          lhs(50,14)=993043270.375;
          lhs(50,17)=1850200707.4464;
          lhs(50,19)=1797790397.0446;
          lhs(50,20)=657045180.5073;
          lhs(50,25)=516149871.8172;
          lhs(50,29)=3467413177.6587;
          lhs(50,32)=7980744561.7547;
          lhs(50,34)=8339346503.7719;
          lhs(50,35)=2890879220.1052;
          lhs(50,39)=3574051601.1613;
          lhs(50,42)=18013564464.9167;
          lhs(50,44)=24013985478.3307;
          lhs(50,45)=8510257232.4101;
          lhs(50,48)=13206083400.4778;
          lhs(50,50)=39104629200;
          lhs(50,51)=17025945402.864;
          lhs(50,53)=22857262565.4927;
          lhs(50,54)=21936941780.469;
          lhs(50,55)=11285978645.1114;
          lhs(51,5)=3753320.9755;
          lhs(51,10)=49098976.2041;
          lhs(51,14)=281346350.3423;
          lhs(51,17)=870770466.312;
          lhs(51,19)=1439802547.0897;
          lhs(51,20)=1006549959.5003;
          lhs(51,25)=59759188.0003;
          lhs(51,29)=701871022.4153;
          lhs(51,32)=3097915859.6865;
          lhs(51,34)=6246399426.5981;
          lhs(51,35)=4953471453.8346;
          lhs(51,39)=539467508.5288;
          lhs(51,42)=4949443648.8831;
          lhs(51,44)=14968471524.2244;
          lhs(51,45)=14983677037.83;
          lhs(51,48)=2712972683.2207;
          lhs(51,50)=17025945402.864;
          lhs(51,51)=26301791880;
          lhs(51,53)=7172537642.2452;
          lhs(51,54)=23148118473.5223;
          lhs(51,55)=7830007777.4355;
          lhs(53,5)=46620675.6901;
          lhs(53,10)=140104632.6768;
          lhs(53,14)=300567132.8641;
          lhs(53,17)=419586081.2108;
          lhs(53,19)=332642700;
          lhs(53,20)=107131656.7046;
          lhs(53,25)=699455809.79;
          lhs(53,29)=1818206351.9133;
          lhs(53,32)=2740095593.3343;
          lhs(53,34)=2198468187.746;
          lhs(53,35)=659465763.0061;
          lhs(53,39)=4975461437.0094;
          lhs(53,42)=10042461000;
          lhs(53,44)=8996377660.438;
          lhs(53,45)=2625740327.7354;
          lhs(53,48)=18082628247.4424;
          lhs(53,50)=22857262565.4927;
          lhs(53,51)=7172537642.2452;
          lhs(53,53)=30750969600;
          lhs(53,54)=13781485528.7;
          lhs(53,55)=15798389148.6235;
          lhs(54,5)=4255865.9536;
          lhs(54,10)=37616895.9451;
          lhs(54,14)=153570334.6216;
          lhs(54,17)=366004472.0095;
          lhs(54,19)=490960014.9113;
          lhs(54,20)=288495675;
          lhs(54,25)=78185250;
          lhs(54,29)=680938552.5758;
          lhs(54,32)=2093470801.5613;
          lhs(54,34)=3212666657.3166;
          lhs(54,35)=2053745729.8917;
          lhs(54,39)=796673401.4956;
          lhs(54,42)=5627367830.2375;
          lhs(54,44)=11712014138.6685;
          lhs(54,45)=8859341128.8642;
          lhs(54,48)=4549025022.2853;
          lhs(54,50)=21936941780.469;
          lhs(54,51)=23148118473.5223;
          lhs(54,53)=13781485528.7;
          lhs(54,54)=33267113100;
          lhs(54,55)=17306268217.8919;
          lhs(55,5)=4705036.842;
          lhs(55,10)=16634817.2837;
          lhs(55,14)=43010526.2156;
          lhs(55,17)=75280589.4727;
          lhs(55,19)=79862301.3628;
          lhs(55,20)=38823912.3022;
          lhs(55,25)=100843220.4851;
          lhs(55,29)=329352578.9431;
          lhs(55,32)=646516120.3778;
          lhs(55,34)=734423151.8869;
          lhs(55,35)=373678486.4369;
          lhs(55,39)=1022486368.1032;
          lhs(55,42)=2833597343.3621;
          lhs(55,44)=3816070032.6464;
          lhs(55,45)=2155875551.0414;
          lhs(55,48)=5551967397.4776;
          lhs(55,50)=11285978645.1114;
          lhs(55,51)=7830007777.4355;
          lhs(55,53)=15798389148.6235;
          lhs(55,54)=17306268217.8919;
          lhs(55,55)=18624305700;
        } break;
        default:
        {
            throw Common::NotImplementedException (FromHere(),"computeLhsTetra Correction Function not implemented for order " + StringOps::to_str(solOrder) + ".");
        }
    }
  return lhs;
}
//////////////////////////////////////////////////////////////////////////////
void VCJH::setup()
{
  CFAUTOTRACE;

  // setup parent class
  BaseCorrectionFunction::setup();

}

void VCJH::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  BaseCorrectionFunction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////
      
void VCJH::configure ( Config::ConfigArgs& args )
{
  // configure parent class
  BaseCorrectionFunction::configure(args);
}
      
//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
