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
RealVector VCJH::computeIntRHS(const CFPolyOrder::Type solOrder, CFuint f, CFuint j)
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

  dubasis=computeDubiner(solOrder,ksi[f],eta[f]);

  for (CFuint i = 0; i < nbrSolPnts; ++i)            
    {
      rhs[i]  = dubasis[i]*weights[j]*G[f]; //computeDubiner(P,ksi,eta) in which we define a=f(ksi,eta) b=f(ksi,eta)
    }

  return rhs;
}
//////////////////////////////////////////////////////////////////////////////
RealVector VCJH::computeDubiner(const CFPolyOrder::Type solOrder,CFreal ksi, CFreal eta)
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
RealVector VCJH::computeSigmas(const CFPolyOrder::Type solOrder,CFuint f, CFuint j, CFreal ksi, CFreal eta, CFreal cfactor)
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
  L=computeLhs(solOrder, ksi,eta);

  //L = cfactor*lhs + I;

  for (CFuint i = 0; i < nbrSolPnts; ++i)
    {
      for (CFuint j = 0; j < nbrSolPnts; ++j)
      {
        L(i,j)= cfactor*L(i,j);
      }
      L(i,i) += 1.;
    }
  
  R=computeIntRHS(solOrder, f, j);
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
RealMatrix VCJH::computeLhs(const CFPolyOrder::Type solOrder, CFreal ksi, CFreal eta)
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
            throw Common::NotImplementedException (FromHere(),"computeLhs Correction Function not implemented for order " + StringOps::to_str(solOrder) + ".");
        }
    }
  return lhs;
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
  dubasis=computeDubiner(solOrder,ksi,eta);
  for (CFuint f = 0; f < 3; ++f)
    {
			for (CFuint j = 0; j < solOrder+1; ++j)
      {
       /* if (f==0) {
          sigmas=computeSigmas(solOrder,2, 2-j, eta,ksi, cfactor);
          dubasis=computeDubiner(solOrder,eta,ksi);
          for (CFuint k = 0; k < nbrSolPnts; ++k)
          {
            phi(f,j)=phi(f,j) + sigmas[k]*dubasis[k];
          }
        }*/
        //else
        {
          sigmas=computeSigmas(solOrder,f, j, ksi,eta, cfactor);
            dubasis=computeDubiner(solOrder,ksi,eta);

          for (CFuint k = 0; k < nbrSolPnts; ++k)
          {
            phi(f,j)=phi(f,j) + sigmas[k]*dubasis[k];
          }
        }
      }
    }
  return (phi*16*4);
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
