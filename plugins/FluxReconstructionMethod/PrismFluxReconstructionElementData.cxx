#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/PrismFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"

#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"


//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////

PrismFluxReconstructionElementData::PrismFluxReconstructionElementData() :
  FluxReconstructionElementData()
{
  m_shape = CFGeoShape::PRISM;
  m_dimensionality = DIM_3D;
}

//////////////////////////////////////////////////////////////////////

PrismFluxReconstructionElementData::PrismFluxReconstructionElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::PRISM;
  m_dimensionality = DIM_3D;
  m_polyOrder = polyOrder;
  
  m_solPntsLocalCoord1D.resize(polyOrder+1);
  m_flxPntsLocalCoord1D.resize(polyOrder+1);
  // Use a default solution and flux point distribution: Gauss Legendre.  
  std::vector<CFreal> coords;
  coords.resize(polyOrder+1);
  switch(polyOrder)
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
                                      + StringOps::to_str(polyOrder) + ".");
      }
    }
  m_solPntsLocalCoord1D = coords;
  m_flxPntsLocalCoord1D = coords;

CFuint rows = (0.5*(polyOrder+1)*(polyOrder+2)) ;
solPntsLocalCoord2D.resize(rows, std::vector< CFreal > (2));
switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
                        solPntsLocalCoord2D[0][0] = 0.333333333333333; 
                        solPntsLocalCoord2D[0][1] = 0.333333333333333;  
      } break;
      case CFPolyOrder::ORDER1:
      {
                        solPntsLocalCoord2D[0][0] = 0.1666666666667; 
                        solPntsLocalCoord2D[0][1] = 0.1666666666667;
                        solPntsLocalCoord2D[1][0] = 0.6666666666667; 
                        solPntsLocalCoord2D[1][1] = 0.1666666666667;  
                        solPntsLocalCoord2D[2][0] = 0.1666666666667; 
                        solPntsLocalCoord2D[2][1] = 0.6666666666667;

    } break;

      case CFPolyOrder::ORDER2:              
      {
                        solPntsLocalCoord2D[0][0] =  0.091576213509780; 
                      solPntsLocalCoord2D[0][1] =  0.091576213509780; 
                      
                      solPntsLocalCoord2D[1][0] =  0.445948490915964; 
                      solPntsLocalCoord2D[1][1] =  0.108103018168071; 

                      solPntsLocalCoord2D[2][0] =  0.816847572980440; 
                      solPntsLocalCoord2D[2][1] =  0.091576213509780;  

                      solPntsLocalCoord2D[3][0] =  0.108103018168071; 
                      solPntsLocalCoord2D[3][1] =  0.445948490915964;
                      
                      solPntsLocalCoord2D[4][0] =  0.445948490915964; 
                      solPntsLocalCoord2D[4][1] =  0.445948490915964;
                      
                      solPntsLocalCoord2D[5][0] =  0.091576213509780; 
                      solPntsLocalCoord2D[5][1] =  0.816847572980440;


      } break;

      case CFPolyOrder::ORDER3:          
      {
                        solPntsLocalCoord2D[0][0] = 0.055564052669793;
                        solPntsLocalCoord2D[0][1] = 0.055564052669793;

                        solPntsLocalCoord2D[1][0] = 0.295533711735893; 
                        solPntsLocalCoord2D[1][1] = 0.070255540518384;  
                        
                        solPntsLocalCoord2D[2][0] = 0.634210747745723; 
                        solPntsLocalCoord2D[2][1] = 0.070255540518384;  	

                        solPntsLocalCoord2D[3][0] = 0.888871894660413; 
                        solPntsLocalCoord2D[3][1] = 0.055564052669793;	

                        solPntsLocalCoord2D[4][0] = 0.070255540518384; 
                        solPntsLocalCoord2D[4][1] = 0.295533711735893; 	

                        solPntsLocalCoord2D[5][0] = 0.333333333333333; 
                        solPntsLocalCoord2D[5][1] = 0.333333333333333;

                        solPntsLocalCoord2D[6][0] = 0.634210747745723; 
                        solPntsLocalCoord2D[6][1] = 0.295533711735893; 

                        solPntsLocalCoord2D[7][0] = 0.070255540518384; 
                        solPntsLocalCoord2D[7][1] = 0.634210747745723;	
                        
                        solPntsLocalCoord2D[8][0] = 0.295533711735893; 
                        solPntsLocalCoord2D[8][1] = 0.634210747745723; 
                        
                        solPntsLocalCoord2D[9][0] = 0.055564052669793; 
                        solPntsLocalCoord2D[9][1] = 0.888871894660413; 


      } break;
        case CFPolyOrder::ORDER4:   
      {
                        solPntsLocalCoord2D[0][0] = 0.035870877695734; 
                        solPntsLocalCoord2D[0][1] = 0.035870877695734;

                        solPntsLocalCoord2D[1][0] = 0.201503881881800; 
                        solPntsLocalCoord2D[1][1] =  0.047312487011716;

                        solPntsLocalCoord2D[2][0] = 0.474308787777079;
                        solPntsLocalCoord2D[2][1] = 0.051382424445843; 

                        solPntsLocalCoord2D[3][0] = 0.751183631106484; 
                        solPntsLocalCoord2D[3][1] = 0.047312487011716;

                        solPntsLocalCoord2D[4][0] = 0.928258244608533; 
                        solPntsLocalCoord2D[4][1] = 0.035870877695734;


                        solPntsLocalCoord2D[5][0] = 0.047312487011716;
                        solPntsLocalCoord2D[5][1] = 0.201503881881800;

                        solPntsLocalCoord2D[6][0] = 0.241729395767967;
                        solPntsLocalCoord2D[6][1] = 0.241729395767967;

                        solPntsLocalCoord2D[7][0] = 0.516541208464066; 
                        solPntsLocalCoord2D[7][1] = 0.241729395767967; 

                        solPntsLocalCoord2D[8][0] = 0.751183631106484; 
                        solPntsLocalCoord2D[8][1] = 0.201503881881800;

                        solPntsLocalCoord2D[9][0] = 0.051382424445843;
                        solPntsLocalCoord2D[9][1] = 0.474308787777079;

                        solPntsLocalCoord2D[10][0] = 0.241729395767967; 
                        solPntsLocalCoord2D[10][1] = 0.516541208464066; 

                        solPntsLocalCoord2D[11][0] = 0.474308787777079; 
                        solPntsLocalCoord2D[11][1] = 0.474308787777079;

                        solPntsLocalCoord2D[12][0] = 0.047312487011716; 
                        solPntsLocalCoord2D[12][1] = 0.751183631106484;

                        solPntsLocalCoord2D[13][0] =  0.201503881881800;
                        solPntsLocalCoord2D[13][1] =  0.751183631106484;

                        solPntsLocalCoord2D[14][0] = 0.035870877695734; 
                        solPntsLocalCoord2D[14][1] = 0.928258244608533; 

      } break;

      case CFPolyOrder::ORDER5:
      {
                        solPntsLocalCoord2D[0][0] = 0.028112952182664; 
                        solPntsLocalCoord2D[0][1] = 0.028112952182664; 
 
                        solPntsLocalCoord2D[1][0] = 0.148565812270887; 
                        solPntsLocalCoord2D[1][1] = 0.033533207700614; 

                        solPntsLocalCoord2D[2][0] = 0.357196298615681; 
                        solPntsLocalCoord2D[2][1] = 0.037824789609186;

                        solPntsLocalCoord2D[3][0] = 0.604978911775132; 
                        solPntsLocalCoord2D[3][1] = 0.037824789609186;

                        solPntsLocalCoord2D[4][0] = 0.817900980028499; 
                        solPntsLocalCoord2D[4][1] = 0.033533207700614;

                        solPntsLocalCoord2D[5][0] = 0.943774095634672;
                        solPntsLocalCoord2D[5][1] = 0.028112952182664;

                        solPntsLocalCoord2D[6][0] = 0.033533207700614;
                        solPntsLocalCoord2D[6][1] = 0.148565812270887;

                        solPntsLocalCoord2D[7][0] = 0.177139098469317;
                        solPntsLocalCoord2D[7][1] = 0.177139098469317;

                        solPntsLocalCoord2D[8][0] = 0.405508595867433;
                        solPntsLocalCoord2D[8][1] = 0.188982808265134;

                        solPntsLocalCoord2D[9][0] = 0.645721803061365; 
                        solPntsLocalCoord2D[9][1] = 0.177139098469317;

                        solPntsLocalCoord2D[10][0] = 0.817900980028499;
                        solPntsLocalCoord2D[10][1] = 0.148565812270887;


                        solPntsLocalCoord2D[11][0] = 0.037824789609186; 
                        solPntsLocalCoord2D[11][1] = 0.357196298615681; 

                        solPntsLocalCoord2D[12][0] = 0.188982808265134;
                        solPntsLocalCoord2D[12][1] = 0.405508595867433;

                        solPntsLocalCoord2D[13][0] = 0.405508595867433;
                        solPntsLocalCoord2D[13][1] = 0.405508595867433;

                        solPntsLocalCoord2D[14][0] = 0.604978911775132; 
                        solPntsLocalCoord2D[14][1] = 0.357196298615681; 


                        solPntsLocalCoord2D[15][0] = 0.037824789609186; 
                        solPntsLocalCoord2D[15][1] = 0.604978911775132;

                        solPntsLocalCoord2D[16][0] = 0.177139098469317;
                        solPntsLocalCoord2D[16][1] = 0.645721803061365;

                        solPntsLocalCoord2D[17][0] = 0.357196298615681;
                        solPntsLocalCoord2D[17][1] = 0.604978911775132;

                        solPntsLocalCoord2D[18][0] = 0.033533207700614; 
                        solPntsLocalCoord2D[18][1] = 0.817900980028499;

                        solPntsLocalCoord2D[19][0] = 0.148565812270887;
                        solPntsLocalCoord2D[19][1] = 0.817900980028499;

                        solPntsLocalCoord2D[20][0] = 0.028112952182664;
                        solPntsLocalCoord2D[20][1] = 0.943774095634672;

      } break;

      default:
      {
        throw Common::NotImplementedException (FromHere(),"Quadrature Points not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
    }


  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////

PrismFluxReconstructionElementData::PrismFluxReconstructionElementData(CFPolyOrder::Type polyOrder, 
								     Common::SafePtr< BasePointDistribution > solPntDist, 
								     Common::SafePtr< BasePointDistribution > flxPntDist)
{
  m_shape = CFGeoShape::PRISM;
  m_dimensionality = DIM_3D;
  m_polyOrder = polyOrder;
  m_solPntsLocalCoord1D = solPntDist->getLocalCoords1D(polyOrder);
  m_flxPntsLocalCoord1D = flxPntDist->getLocalCoords1D(polyOrder);
  m_solPntDistribution = solPntDist;
  m_flxPntDistribution = flxPntDist;

CFuint rows = (0.5*(polyOrder+1)*(polyOrder+2)) ;
solPntsLocalCoord2D.resize(rows, std::vector< CFreal > (2));
switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
                        solPntsLocalCoord2D[0][0] = 0.333333333333333; 
                        solPntsLocalCoord2D[0][1] = 0.333333333333333;  
      } break;
      case CFPolyOrder::ORDER1:
      {
                        solPntsLocalCoord2D[0][0] = 0.1666666666667; 
                        solPntsLocalCoord2D[0][1] = 0.1666666666667;
                        solPntsLocalCoord2D[1][0] = 0.6666666666667; 
                        solPntsLocalCoord2D[1][1] = 0.1666666666667;  
                        solPntsLocalCoord2D[2][0] = 0.1666666666667; 
                        solPntsLocalCoord2D[2][1] = 0.6666666666667;

    } break;

      case CFPolyOrder::ORDER2:              
      {
                        solPntsLocalCoord2D[0][0] =  0.091576213509780; 
                      solPntsLocalCoord2D[0][1] =  0.091576213509780; 
                      
                      solPntsLocalCoord2D[1][0] =  0.445948490915964; 
                      solPntsLocalCoord2D[1][1] =  0.108103018168071; 

                      solPntsLocalCoord2D[2][0] =  0.816847572980440; 
                      solPntsLocalCoord2D[2][1] =  0.091576213509780;  

                      solPntsLocalCoord2D[3][0] =  0.108103018168071; 
                      solPntsLocalCoord2D[3][1] =  0.445948490915964;
                      
                      solPntsLocalCoord2D[4][0] =  0.445948490915964; 
                      solPntsLocalCoord2D[4][1] =  0.445948490915964;
                      
                      solPntsLocalCoord2D[5][0] =  0.091576213509780; 
                      solPntsLocalCoord2D[5][1] =  0.816847572980440;


      } break;

      case CFPolyOrder::ORDER3:          
      {
                        solPntsLocalCoord2D[0][0] = 0.055564052669793;
                        solPntsLocalCoord2D[0][1] = 0.055564052669793;

                        solPntsLocalCoord2D[1][0] = 0.295533711735893; 
                        solPntsLocalCoord2D[1][1] = 0.070255540518384;  
                        
                        solPntsLocalCoord2D[2][0] = 0.634210747745723; 
                        solPntsLocalCoord2D[2][1] = 0.070255540518384;  	

                        solPntsLocalCoord2D[3][0] = 0.888871894660413; 
                        solPntsLocalCoord2D[3][1] = 0.055564052669793;	

                        solPntsLocalCoord2D[4][0] = 0.070255540518384; 
                        solPntsLocalCoord2D[4][1] = 0.295533711735893; 	

                        solPntsLocalCoord2D[5][0] = 0.333333333333333; 
                        solPntsLocalCoord2D[5][1] = 0.333333333333333;

                        solPntsLocalCoord2D[6][0] = 0.634210747745723; 
                        solPntsLocalCoord2D[6][1] = 0.295533711735893; 

                        solPntsLocalCoord2D[7][0] = 0.070255540518384; 
                        solPntsLocalCoord2D[7][1] = 0.634210747745723;	
                        
                        solPntsLocalCoord2D[8][0] = 0.295533711735893; 
                        solPntsLocalCoord2D[8][1] = 0.634210747745723; 
                        
                        solPntsLocalCoord2D[9][0] = 0.055564052669793; 
                        solPntsLocalCoord2D[9][1] = 0.888871894660413; 


      } break;
        case CFPolyOrder::ORDER4:   
      {
                        solPntsLocalCoord2D[0][0] = 0.035870877695734; 
                        solPntsLocalCoord2D[0][1] = 0.035870877695734;

                        solPntsLocalCoord2D[1][0] = 0.201503881881800; 
                        solPntsLocalCoord2D[1][1] =  0.047312487011716;

                        solPntsLocalCoord2D[2][0] = 0.474308787777079;
                        solPntsLocalCoord2D[2][1] = 0.051382424445843; 

                        solPntsLocalCoord2D[3][0] = 0.751183631106484; 
                        solPntsLocalCoord2D[3][1] = 0.047312487011716;

                        solPntsLocalCoord2D[4][0] = 0.928258244608533; 
                        solPntsLocalCoord2D[4][1] = 0.035870877695734;


                        solPntsLocalCoord2D[5][0] = 0.047312487011716;
                        solPntsLocalCoord2D[5][1] = 0.201503881881800;

                        solPntsLocalCoord2D[6][0] = 0.241729395767967;
                        solPntsLocalCoord2D[6][1] = 0.241729395767967;

                        solPntsLocalCoord2D[7][0] = 0.516541208464066; 
                        solPntsLocalCoord2D[7][1] = 0.241729395767967; 

                        solPntsLocalCoord2D[8][0] = 0.751183631106484; 
                        solPntsLocalCoord2D[8][1] = 0.201503881881800;

                        solPntsLocalCoord2D[9][0] = 0.051382424445843;
                        solPntsLocalCoord2D[9][1] = 0.474308787777079;

                        solPntsLocalCoord2D[10][0] = 0.241729395767967; 
                        solPntsLocalCoord2D[10][1] = 0.516541208464066; 

                        solPntsLocalCoord2D[11][0] = 0.474308787777079; 
                        solPntsLocalCoord2D[11][1] = 0.474308787777079;

                        solPntsLocalCoord2D[12][0] = 0.047312487011716; 
                        solPntsLocalCoord2D[12][1] = 0.751183631106484;

                        solPntsLocalCoord2D[13][0] =  0.201503881881800;
                        solPntsLocalCoord2D[13][1] =  0.751183631106484;

                        solPntsLocalCoord2D[14][0] = 0.035870877695734; 
                        solPntsLocalCoord2D[14][1] = 0.928258244608533; 

      } break;

      case CFPolyOrder::ORDER5:
      {
                        solPntsLocalCoord2D[0][0] = 0.028112952182664; 
                        solPntsLocalCoord2D[0][1] = 0.028112952182664; 
 
                        solPntsLocalCoord2D[1][0] = 0.148565812270887; 
                        solPntsLocalCoord2D[1][1] = 0.033533207700614; 

                        solPntsLocalCoord2D[2][0] = 0.357196298615681; 
                        solPntsLocalCoord2D[2][1] = 0.037824789609186;

                        solPntsLocalCoord2D[3][0] = 0.604978911775132; 
                        solPntsLocalCoord2D[3][1] = 0.037824789609186;

                        solPntsLocalCoord2D[4][0] = 0.817900980028499; 
                        solPntsLocalCoord2D[4][1] = 0.033533207700614;

                        solPntsLocalCoord2D[5][0] = 0.943774095634672;
                        solPntsLocalCoord2D[5][1] = 0.028112952182664;

                        solPntsLocalCoord2D[6][0] = 0.033533207700614;
                        solPntsLocalCoord2D[6][1] = 0.148565812270887;

                        solPntsLocalCoord2D[7][0] = 0.177139098469317;
                        solPntsLocalCoord2D[7][1] = 0.177139098469317;

                        solPntsLocalCoord2D[8][0] = 0.405508595867433;
                        solPntsLocalCoord2D[8][1] = 0.188982808265134;

                        solPntsLocalCoord2D[9][0] = 0.645721803061365; 
                        solPntsLocalCoord2D[9][1] = 0.177139098469317;

                        solPntsLocalCoord2D[10][0] = 0.817900980028499;
                        solPntsLocalCoord2D[10][1] = 0.148565812270887;


                        solPntsLocalCoord2D[11][0] = 0.037824789609186; 
                        solPntsLocalCoord2D[11][1] = 0.357196298615681; 

                        solPntsLocalCoord2D[12][0] = 0.188982808265134;
                        solPntsLocalCoord2D[12][1] = 0.405508595867433;

                        solPntsLocalCoord2D[13][0] = 0.405508595867433;
                        solPntsLocalCoord2D[13][1] = 0.405508595867433;

                        solPntsLocalCoord2D[14][0] = 0.604978911775132; 
                        solPntsLocalCoord2D[14][1] = 0.357196298615681; 


                        solPntsLocalCoord2D[15][0] = 0.037824789609186; 
                        solPntsLocalCoord2D[15][1] = 0.604978911775132;

                        solPntsLocalCoord2D[16][0] = 0.177139098469317;
                        solPntsLocalCoord2D[16][1] = 0.645721803061365;

                        solPntsLocalCoord2D[17][0] = 0.357196298615681;
                        solPntsLocalCoord2D[17][1] = 0.604978911775132;

                        solPntsLocalCoord2D[18][0] = 0.033533207700614; 
                        solPntsLocalCoord2D[18][1] = 0.817900980028499;

                        solPntsLocalCoord2D[19][0] = 0.148565812270887;
                        solPntsLocalCoord2D[19][1] = 0.817900980028499;

                        solPntsLocalCoord2D[20][0] = 0.028112952182664;
                        solPntsLocalCoord2D[20][1] = 0.943774095634672;

      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Quadrature Points not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
    }

  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////

PrismFluxReconstructionElementData::~PrismFluxReconstructionElementData()
{
}

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
std::vector<CFreal> PrismFluxReconstructionElementData::getPercentage(CFPolyOrder::Type solOrder){
    CFAUTOTRACE;
    std::vector<CFreal> alpha;
    alpha.resize(solOrder+1);

    //std::vector< CFreal > m_flxPntsLocalCoord1D = getLocalCoords1D(solOrder);
    for (CFuint iFlx=0; iFlx< (solOrder+1) ; ++iFlx){
        alpha[iFlx] = (m_flxPntsLocalCoord1D[iFlx]+1)/2;
    }
		return alpha;
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFlxPntsLocalCoords()
{
  CFAUTOTRACE;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
  const CFuint nbrFlxPnts2DTriag = (m_polyOrder+1)*(m_polyOrder+2)/2;
  // set flux point local coordinates
  m_flxPntsLocalCoords.resize(0);
  std::vector<CFreal> alpha = getPercentage(m_polyOrder);

    RealVector flxCoords(3);
    //faces ZTA=-1 and 1 --> triag faces
    flxCoords[ZTA] = -1;
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts2DTriag; ++iFlx)
    {
      flxCoords[KSI] = solPntsLocalCoord2D[iFlx][KSI];
      flxCoords[ETA] = solPntsLocalCoord2D[iFlx][ETA];
      m_flxPntsLocalCoords.push_back(flxCoords);
    }

    flxCoords[ZTA] = 1;
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts2DTriag; ++iFlx)
    {
      flxCoords[KSI] = solPntsLocalCoord2D[iFlx][KSI];
      flxCoords[ETA] = solPntsLocalCoord2D[iFlx][ETA];
      m_flxPntsLocalCoords.push_back(flxCoords);
    }
    // 3 rectangular faces 

    //face at eta=0
    flxCoords[ETA] = 0;
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
        for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
        {
            flxCoords[KSI] = (m_flxPntsLocalCoord1D[iKsi]+1.)/2.;
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
    // (Oblique) face at ksi+eta=1
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
        for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
        {
            flxCoords[KSI] = 1-alpha[iEta];
            flxCoords[ETA] = 0+alpha[iEta];
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }

    //face at ksi=0
    flxCoords[KSI] = 0;
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
        for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
        {
            flxCoords[ETA] = (m_flxPntsLocalCoord1D[nbrFlxPnts1D-1-iEta]+1.)/2.;
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
  cf_assert(m_flxPntsLocalCoords.size() == (3*nbrFlxPnts1D*nbrFlxPnts1D)+((m_polyOrder+1)*(m_polyOrder+2)));
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createSolPntsLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  const CFuint nbrSolPnts2D = (m_polyOrder+1)*(m_polyOrder+2)/2;
  // set solution point local coordinates
  m_solPntsLocalCoords.resize(0);
  for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
  {
    for (CFuint iFlx = 0; iFlx < nbrSolPnts2D; ++iFlx)
    {
      RealVector solCoords(3);
      solCoords[KSI] = solPntsLocalCoord2D[iFlx][KSI];
      solCoords[ETA] = solPntsLocalCoord2D[iFlx][ETA];
      solCoords[ZTA] = m_solPntsLocalCoord1D[iZta];
      m_solPntsLocalCoords.push_back(solCoords);
    }
  }
    cf_assert(m_solPntsLocalCoords.size() == ((m_polyOrder+1)*(m_polyOrder+1)*(m_polyOrder+2)/2));
}

//////////////////////////////////////////////////////////////////////

// !!!!!!!!! To be chceked

void PrismFluxReconstructionElementData::createFaceFlxPntsFaceLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // set face flux point face local coordinates
  m_faceFlxPntsFaceLocalCoords.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      RealVector flxCoord(2);
      flxCoord[KSI] = m_flxPntsLocalCoord1D[iKsi];
      flxCoord[ETA] = m_flxPntsLocalCoord1D[iEta];
      m_faceFlxPntsFaceLocalCoords.push_back(flxCoord);
    }
  }

}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFaceFlxPntsLocalCoordsPerType()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
  // number of solution points in triag
  const CFuint nbrFlxPntsTriag = (m_polyOrder+1)*(m_polyOrder+2)/2;

  // set face flux point face local coordinates on Triag Face
  m_faceFlxPntsLocalCoordsPerType.resize(2);
  m_faceFlxPntsLocalCoordsPerType[0].resize(0);
  for (CFuint iKsi = 0; iKsi < nbrFlxPntsTriag; ++iKsi)
  {
      RealVector flxCoord(2);
      flxCoord[KSI] = solPntsLocalCoord2D[iKsi][KSI];
      flxCoord[ETA] = solPntsLocalCoord2D[iKsi][ETA];
      (m_faceFlxPntsLocalCoordsPerType[0]).push_back(flxCoord);    
  }


  // set face flux point face local coordinates on Quad Face (reference -1 1 quad)
  m_faceFlxPntsLocalCoordsPerType[1].resize(0);
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      RealVector flxCoord(2);
      flxCoord[KSI] = m_flxPntsLocalCoord1D[iKsi];
      flxCoord[ETA] = m_flxPntsLocalCoord1D[iEta];
      (m_faceFlxPntsLocalCoordsPerType[1]).push_back(flxCoord);
    }
  }

}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFaceIntegrationCoefsPerType()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
  // number of solution points in triag
  const CFuint nbrFlxPntsTriag = (m_polyOrder+1)*(m_polyOrder+2)/2;

  // set face flux point face local coordinates on Triag Face
  m_faceIntegrationCoefsPerType.resize(2);

  m_faceIntegrationCoefsPerType[0].resize(nbrFlxPntsTriag);

  switch (m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
        m_faceIntegrationCoefsPerType[0][0] = 0.5;
    } break;
    
    case CFPolyOrder::ORDER1:
    {
        m_faceIntegrationCoefsPerType[0][0] = 0.166666666666667;
        m_faceIntegrationCoefsPerType[0][1] = 0.166666666666667;
        m_faceIntegrationCoefsPerType[0][2] = 0.166666666666667;
    } break;
    
    case CFPolyOrder::ORDER2:
    {
        m_faceIntegrationCoefsPerType[0][0] = 0.054975871827667;
        m_faceIntegrationCoefsPerType[0][1] = 0.111690794839000;
        m_faceIntegrationCoefsPerType[0][2] = 0.054975871827667;
        m_faceIntegrationCoefsPerType[0][3] = 0.111690794839000;
        m_faceIntegrationCoefsPerType[0][4] = 0.111690794839000;
        m_faceIntegrationCoefsPerType[0][5] = 0.054975871827667;
    } break;
    
    case CFPolyOrder::ORDER3:
    {
        m_faceIntegrationCoefsPerType[0][0] = 0.020977756498325;
        m_faceIntegrationCoefsPerType[0][1] = 0.056049206035444;
        m_faceIntegrationCoefsPerType[0][2] = 0.056049206035444;
        m_faceIntegrationCoefsPerType[0][3] = 0.020977756498325;
        m_faceIntegrationCoefsPerType[0][4] = 0.056049206035444;
        m_faceIntegrationCoefsPerType[0][5] = 0.100771494292365;
        m_faceIntegrationCoefsPerType[0][6] = 0.056049206035444;
        m_faceIntegrationCoefsPerType[0][7] = 0.056049206035444;
        m_faceIntegrationCoefsPerType[0][8] = 0.056049206035444;
        m_faceIntegrationCoefsPerType[0][9] = 0.020977756498325;
    } break;

    case CFPolyOrder::ORDER4:
    {
        m_faceIntegrationCoefsPerType[0][0] = 0.008957727506152;
        m_faceIntegrationCoefsPerType[0][1] = 0.027874905013558;
        m_faceIntegrationCoefsPerType[0][2] = 0.038103031192768;
        m_faceIntegrationCoefsPerType[0][3] = 0.055749810027115;
        m_faceIntegrationCoefsPerType[0][4] = 0.063856097940632;
        m_faceIntegrationCoefsPerType[0][5] = 0.027874905013558;
        m_faceIntegrationCoefsPerType[0][6] = 0.063856097940632;
        m_faceIntegrationCoefsPerType[0][7] = 0.127712195881265;
        m_faceIntegrationCoefsPerType[0][8] = 0.063856097940632;
        m_faceIntegrationCoefsPerType[0][9] = 0.038103031192768;
        m_faceIntegrationCoefsPerType[0][10] = 0.038103031192768;
        m_faceIntegrationCoefsPerType[0][11] = 0.055749810027115;
        m_faceIntegrationCoefsPerType[0][12] = 0.027874905013558;
        m_faceIntegrationCoefsPerType[0][13] = 0.055749810027115;
        m_faceIntegrationCoefsPerType[0][14] = 0.008957727506152;
    } break;

    case CFPolyOrder::ORDER5:
    {
        m_faceIntegrationCoefsPerType[0][0] = 0.005179687348269;
        m_faceIntegrationCoefsPerType[0][1] = 0.014484634686237;
        m_faceIntegrationCoefsPerType[0][2] = 0.023023183297967;
        m_faceIntegrationCoefsPerType[0][3] = 0.023023183297967;
        m_faceIntegrationCoefsPerType[0][4] = 0.028969269372473;
        m_faceIntegrationCoefsPerType[0][5] = 0.014484634686237;
        m_faceIntegrationCoefsPerType[0][6] = 0.037697442163369;
        m_faceIntegrationCoefsPerType[0][7] = 0.023023183297967;
        m_faceIntegrationCoefsPerType[0][8] = 0.037697442163369;
        m_faceIntegrationCoefsPerType[0][9] = 0.075394884326738;
        m_faceIntegrationCoefsPerType[0][10] = 0.048773901186621;
        m_faceIntegrationCoefsPerType[0][11] = 0.037697442163369;
        m_faceIntegrationCoefsPerType[0][12] = 0.014484634686237;
        m_faceIntegrationCoefsPerType[0][13] = 0.028969269372473;
        m_faceIntegrationCoefsPerType[0][14] = 0.028969269372473;
        m_faceIntegrationCoefsPerType[0][15] = 0.046046366595935;
        m_faceIntegrationCoefsPerType[0][16] = 0.023023183297967;
        m_faceIntegrationCoefsPerType[0][17] = 0.037697442163369;
        m_faceIntegrationCoefsPerType[0][18] = 0.014484634686237;
        m_faceIntegrationCoefsPerType[0][19] = 0.075394884326738;
        m_faceIntegrationCoefsPerType[0][20] = 0.005179687348269;
    } break;

    default:
    {
        throw Common::NotImplementedException(FromHere(), "Face Integration Coeff Per Face not implemented for order "
                                              + StringOps::to_str(m_polyOrder) + ".");
    }
  }

  // set face flux point face local coordinates on Quad Face (reference -1 1 quad)

  // number of flux points on a face
  const CFuint nbrFlxPnts = nbrFlxPnts1D*nbrFlxPnts1D;

  // resize m_faceIntegrationCoefsPerType
  m_faceIntegrationCoefsPerType[1].resize(nbrFlxPnts);

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_2D,m_polyOrder);

  // create face node local coordinates
  vector< RealVector > nodeCoord(4);
  nodeCoord[0].resize(2);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[0][ETA] = -1.0;
  nodeCoord[1].resize(2);
  nodeCoord[1][KSI] = +1.0;
  nodeCoord[1][ETA] = -1.0;
  nodeCoord[2].resize(2);
  nodeCoord[2][KSI] = +1.0;
  nodeCoord[2][ETA] = +1.0;
  nodeCoord[3].resize(2);
  nodeCoord[3][KSI] = -1.0;
  nodeCoord[3][ETA] = +1.0;

  // get quadrature point coordinates and wheights
  vector< RealVector > quadPntCoords   = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  vector< CFreal     > quadPntWheights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWheights.size() == nbrQPnts);

  // compute the coefficients for integration over a face
  // loop over flux points
  CFuint iFlx = 0;
  for (CFuint iFlxKsi = 0; iFlxKsi < nbrFlxPnts1D; ++iFlxKsi)
  {
    const CFreal ksiFlx = m_flxPntsLocalCoord1D[iFlxKsi];
    for (CFuint iFlxEta = 0; iFlxEta < nbrFlxPnts1D; ++iFlxEta, ++iFlx)
    {
      const CFreal etaFlx = m_flxPntsLocalCoord1D[iFlxEta];

      m_faceIntegrationCoefsPerType[1][iFlx] = 0.0;
      for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
      {
        // quadrature point local coordinate on the face
        const CFreal ksiQPnt = quadPntCoords[iQPnt][KSI];
        const CFreal etaQPnt = quadPntCoords[iQPnt][ETA];

        // evaluate polynomial value in quadrature point
        CFreal quadPntPolyVal = 1.;
        for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
        {
          if (iFac != iFlxKsi)
          {
            const CFreal ksiFac = m_flxPntsLocalCoord1D[iFac];
            quadPntPolyVal *= (ksiQPnt-ksiFac)/(ksiFlx-ksiFac);
          }
        }
        for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
        {
          if (iFac != iFlxEta)
          {
            const CFreal etaFac = m_flxPntsLocalCoord1D[iFac];
            quadPntPolyVal *= (etaQPnt-etaFac)/(etaFlx-etaFac);
          }
        }

        // add contribution of quadrature point to integration coefficient
        m_faceIntegrationCoefsPerType[1][iFlx] += quadPntWheights[iQPnt]*quadPntPolyVal;
      }
    }
  }

}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createSolPolyExponents()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // define exponents
  m_solPolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D-iKsi; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
      {
        vector< CFint > solPolyExps(3);
        solPolyExps[KSI] = iKsi;
        solPolyExps[ETA] = iEta;
        solPolyExps[ZTA] = iZta;
        m_solPolyExponents.push_back(solPolyExps);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createNodePolyExponents()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrNodes1D = 2;

  // define exponents
  m_nodePolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrNodes1D-iKsi; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrNodes1D; ++iZta)
      {
        vector< CFint > nodePolyExps(3);
        nodePolyExps[KSI] = iKsi;
        nodePolyExps[ETA] = iEta;
        nodePolyExps[ZTA] = iZta;
        m_nodePolyExponents.push_back(nodePolyExps);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFaceFluxPntsConn()
{
  CFAUTOTRACE;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
  const CFuint nbrFlxPnts2D = (m_polyOrder+1)*(m_polyOrder+2)/2; //on triag face

  // resize m_faceFlxPntConn
  m_faceFlxPntConn.resize(5);

  // variable holding the face index
  CFuint faceIdx = 0;

  // zeroth face (triag face zta=-1)
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D-iKsi; ++iEta)
    {
      m_faceFlxPntConn[faceIdx].push_back((m_polyOrder +1)*iEta - (iEta)*(iEta-1)/2 + iKsi);
    }
  }
  ++faceIdx;

  // first face (triag face zta=1)
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx)
  {
    m_faceFlxPntConn[faceIdx].push_back(nbrFlxPnts2D + iFlx); 
  }
  ++faceIdx;

  // second face (eta=0)
  for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
  {
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
      m_faceFlxPntConn[faceIdx].push_back(2*nbrFlxPnts2D + nbrFlxPnts1D*iKsi + iEta);
    }
  }
  ++faceIdx;

  // third face (ksi+eta=1)
  for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
  {
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
      m_faceFlxPntConn[faceIdx].push_back(2*nbrFlxPnts2D + nbrFlxPnts1D*nbrFlxPnts1D + nbrFlxPnts1D*iZta + iEta);
    }
  }
  ++faceIdx;

  // fourth face (ksi=0)
  for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
  {
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
      m_faceFlxPntConn[faceIdx].push_back(2*nbrFlxPnts2D + 2*nbrFlxPnts1D*nbrFlxPnts1D + nbrFlxPnts1D*iZta + (nbrFlxPnts1D-1-iEta));
    }
  }
}

//////////////////////////////////////////////////////////////////////
// !!!! to be checked

void PrismFluxReconstructionElementData::createFaceFluxPntsConnPerOrient()
{
  CFAUTOTRACE;

     // number of faces
  const CFuint nbrFaces = 5;

  const CFuint nbrQuadFaces = 3;

  const CFuint nbrTriagFaces = 2;

  const CFuint nbrTriagFaceFlxPnts = (m_polyOrder+1)*(m_polyOrder+2)/2;
  // number of possible orientations
  const CFuint nbrOrient = 33;
  const CFuint nbrQuadOrient = 24;
  const CFuint nbrTriagOrient = 9;

 // flux point indexes for inverted face
  vector< CFuint > invFlxIdxs;
  for (CFuint iKsi = 0; iKsi < nbrTriagFaceFlxPnts; ++iKsi)
  {
      invFlxIdxs.push_back(iKsi);
  }
  //Build matrix reprensenting flx indices
  const CFuint p = m_polyOrder;
  RealMatrix flxIdx(p+1, p+1);
  CFuint id=0;
  for (CFuint i = 0; i<p+1 ; ++i)
  {
    for(CFuint j=0; j<p+1-i; ++j, ++id)
    {
      flxIdx((p-i),j)=id;
    }
  }
          // flip order flx pnts indices ref matrix
        id=0;
        for (CFuint i = 0; i<p+1 ; ++i)
        {
          for(CFuint j=i; j<p+1; ++j, ++id)
          {
            flxIdx((p-i),(p-j))=id;
          }
        }

        id=0;
        for (CFuint i = 0; i<p+1 ; ++i)
        {
          for(CFuint j=0; j<p+1-i; ++j, ++id)
          {
            invFlxIdxs[id]=flxIdx(p-i,j);
          }
        }

  // resize the variables
  m_faceFlxPntConnPerOrient.resize(nbrOrient);


  // fill the variable
  CFuint iOrient = 0;

/////////////////////Triag Faces

  for (CFuint iFaceL = 0; iFaceL < nbrTriagFaces; ++iFaceL)
  {
    vector< CFuint > faceFlxConnL = m_faceFlxPntConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrTriagFaces; ++iFaceR)
    {
      vector< CFuint > faceFlxConnR = m_faceFlxPntConn[iFaceR];
      for (CFuint iRot = 0; iRot < 3; ++iRot, ++iOrient)
      {
        m_faceFlxPntConnPerOrient[iOrient].resize(2);

        for (CFuint iFlx = 0; iFlx < nbrTriagFaceFlxPnts; ++iFlx)
        {
          m_faceFlxPntConnPerOrient[iOrient][LEFT ].push_back(faceFlxConnL[iFlx            ]);
          m_faceFlxPntConnPerOrient[iOrient][RIGHT].push_back(faceFlxConnR[invFlxIdxs[iFlx]]);
        }

        // rotate flx pnts indices
        id=0;
        for (CFuint j = 0; j<p+1 ; ++j)
        {
          for(CFuint i=j; i<p+1; ++i, ++id)
          {
            invFlxIdxs[id]=flxIdx(i,j);
          }
        }
        
        // rotate flx pnts indices ref matrix
        id=0;
        for (CFuint i = 0; i<p+1 ; ++i)
        {
          for(CFuint j=0; j<p+1-i; ++j, ++id)
          {
            flxIdx((p-i),j)=int(invFlxIdxs[id]);
          }
        }
      }

      // reset flx pnts indices ref matrix
      for (CFuint i = 0; i<p+1 ; ++i)
      {
        for(CFuint j=0; j<p+1-i; ++j, ++id)
        {
          flxIdx((p-i),j)=id;
        }
      }
      // flip order flx pnts indices ref matrix
      id=0;
      for (CFuint i = 0; i<p+1 ; ++i)
      {
        for(CFuint j=i; j<p+1; ++j, ++id)
        {
          flxIdx((p-i),(p-j))=id;
        }
      }
      id=0;
      for (CFuint i = 0; i<p+1 ; ++i)
      {
        for(CFuint j=0; j<p+1-i; ++j, ++id)
        {
          invFlxIdxs[id]=flxIdx(p-i,j);
        }
      }
    }
  }


//////////////////////QUAD FACES

  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  // number of face flux points
  const CFuint nbrQuadFaceFlxPnts = nbrSolPnts1D*nbrSolPnts1D;
  // flux point indexes for inverted face
  vector< CFuint > invFlxIdxsq;
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    const CFuint idxKsi = nbrSolPnts1D - iKsi - 1;
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
    {
      invFlxIdxsq.push_back(nbrSolPnts1D*idxKsi+iEta);
    }
  }
  // number of rotatable flux point groups
  const CFuint nbrRotFlxGroups = nbrQuadFaceFlxPnts/4;
  // maximum number of flux points in a line of flux points
  const CFuint maxNbrFlxPntsInLine = (nbrSolPnts1D+1)/2;
  // storage of flux point rotatable groups
  vector< vector< CFuint > > rotFlxIdxs(nbrRotFlxGroups);
  CFuint iRotGroup = 0;
  for (CFuint iFlxLine = 0; iRotGroup < nbrRotFlxGroups; ++iFlxLine)
  {
    for (CFuint iFlx = 0;
         iFlx < maxNbrFlxPntsInLine && iRotGroup < nbrRotFlxGroups;
         ++iFlx, ++iRotGroup)
    {
      const CFuint idxFlxLine = nbrSolPnts1D - 1 - iFlxLine;
      const CFuint idxFlx     = nbrSolPnts1D - 1 - iFlx    ;

      rotFlxIdxs[iRotGroup].push_back(nbrSolPnts1D*iFlx      +iFlxLine  );
      rotFlxIdxs[iRotGroup].push_back(nbrSolPnts1D*idxFlxLine+iFlx      );
      rotFlxIdxs[iRotGroup].push_back(nbrSolPnts1D*idxFlx    +idxFlxLine);
      rotFlxIdxs[iRotGroup].push_back(nbrSolPnts1D*iFlxLine  +idxFlx    );
    }
  }

  
  for (CFuint iFaceL = 2; iFaceL < nbrQuadFaces+2; ++iFaceL)
  {
    vector< CFuint > faceFlxConnL = m_faceFlxPntConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrQuadFaces+2; ++iFaceR)
    {
      vector< CFuint > faceFlxConnR = m_faceFlxPntConn[iFaceR];
      for (CFuint iRot = 0; iRot < 4; ++iRot, ++iOrient)
      {
        m_faceFlxPntConnPerOrient[iOrient].resize(2);
        for (CFuint iFlx = 0; iFlx < nbrQuadFaceFlxPnts; ++iFlx)
        {
          m_faceFlxPntConnPerOrient[iOrient][LEFT ]
              .push_back(faceFlxConnL[iFlx            ]);
          m_faceFlxPntConnPerOrient[iOrient][RIGHT]
              .push_back(faceFlxConnR[invFlxIdxsq[iFlx]]);
        }
        // rotate the right face
        for (CFuint iRotGroup = 0; iRotGroup < nbrRotFlxGroups; ++iRotGroup)
        {
          // indexes of flux points to be rotated
          const CFuint flx0Idx = rotFlxIdxs[iRotGroup][0];
          const CFuint flx1Idx = rotFlxIdxs[iRotGroup][1];
          const CFuint flx2Idx = rotFlxIdxs[iRotGroup][2];
          const CFuint flx3Idx = rotFlxIdxs[iRotGroup][3];
          //CFLog(VERBOSE,"flxrotIdx" << flx0Idx << flx1Idx << flx2Idx << flx3Idx << "\n");

          // rotate flux points
          const CFuint swap = invFlxIdxsq[flx3Idx];
          invFlxIdxsq[flx3Idx] = invFlxIdxsq[flx2Idx];
          invFlxIdxsq[flx2Idx] = invFlxIdxsq[flx1Idx];
          invFlxIdxsq[flx1Idx] = invFlxIdxsq[flx0Idx];
          invFlxIdxsq[flx0Idx] = swap;
        }
      }
    }
  }
  
  CFLog(VERBOSE,"Prismelemdata --- createFaceFluxPntsConnPerOrient end \n");

}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createCellNodeCoords()
{
  CFAUTOTRACE;

  m_cellNodeCoords.resize(6);

  // first node
  m_cellNodeCoords[0].resize(3);
  m_cellNodeCoords[0][KSI] = +0.;
  m_cellNodeCoords[0][ETA] = +0.;
  m_cellNodeCoords[0][ZTA] = -1.0;

  // second node
  m_cellNodeCoords[1].resize(3);
  m_cellNodeCoords[1][KSI] = +1.0;
  m_cellNodeCoords[1][ETA] = +0.0;
  m_cellNodeCoords[1][ZTA] = -1.0;

  // third node
  m_cellNodeCoords[2].resize(3);
  m_cellNodeCoords[2][KSI] = +0.0;
  m_cellNodeCoords[2][ETA] = +1.0;
  m_cellNodeCoords[2][ZTA] = -1.0;

  // fourth node
  m_cellNodeCoords[3].resize(3);
  m_cellNodeCoords[3][KSI] = +0.0;
  m_cellNodeCoords[3][ETA] = +0.0;
  m_cellNodeCoords[3][ZTA] = +1.0;

  // fifth node
  m_cellNodeCoords[4].resize(3);
  m_cellNodeCoords[4][KSI] = +1.0;
  m_cellNodeCoords[4][ETA] = +0.0;
  m_cellNodeCoords[4][ZTA] = +1.0;

  // sixth node
  m_cellNodeCoords[5].resize(3);
  m_cellNodeCoords[5][KSI] = +0.0;
  m_cellNodeCoords[5][ETA] = +1.0;
  m_cellNodeCoords[5][ZTA] = +1.0;

}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFaceNodeConnectivity()
{
  CFAUTOTRACE;

  m_faceNodeConn.resize(5);

  m_faceNodeConn[0].resize(3);
  m_faceNodeConn[0][0] = 0;
  m_faceNodeConn[0][1] = 2;
  m_faceNodeConn[0][2] = 1;


  m_faceNodeConn[1].resize(3);
  m_faceNodeConn[1][0] = 3;
  m_faceNodeConn[1][1] = 4;
  m_faceNodeConn[1][2] = 5;

  m_faceNodeConn[2].resize(4);
  m_faceNodeConn[2][0] = 0;
  m_faceNodeConn[2][1] = 1;
  m_faceNodeConn[2][2] = 4;
  m_faceNodeConn[2][3] = 3;

  m_faceNodeConn[3].resize(4);
  m_faceNodeConn[3][0] = 1;
  m_faceNodeConn[3][1] = 2;
  m_faceNodeConn[3][2] = 5;
  m_faceNodeConn[3][3] = 4;

  m_faceNodeConn[4].resize(4);
  m_faceNodeConn[4][0] = 0;
  m_faceNodeConn[4][1] = 3;
  m_faceNodeConn[4][2] = 5;
  m_faceNodeConn[4][3] = 2;

}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFaceMappedCoordDir()
{
  CFAUTOTRACE;

  m_faceMappedCoordDir.resize(5);

  m_faceMappedCoordDir[0] = -1;
  m_faceMappedCoordDir[1] = +1;
  m_faceMappedCoordDir[2] = -2;
  m_faceMappedCoordDir[3] = -2;
  m_faceMappedCoordDir[4] = -2;

}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFluxPntFluxDim()
{
  CFAUTOTRACE;

  m_flxPntFlxDim.resize(3*(m_flxPntsLocalCoord1D.size()*m_flxPntsLocalCoord1D.size())+2 * (solPntsLocalCoord2D.size()));
  
  for (CFuint iFace = 0; iFace < m_faceFlxPntConn.size(); ++iFace)
  {
    for (CFuint iFlx = 0; iFlx < m_faceFlxPntConn[iFace].size(); ++iFlx)
    {
     if (iFace == 0)
      {
        m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 2; // z direction
      }
      else if (iFace == 1)
      {
	      m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 2; // z direction
      }
      else if (iFace == 2)
      {
	      m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 1; // y direction
      }
      else if (iFace == 3)
      {
	      m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 3;
      }
      else if (iFace == 4)
      {
	      m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 0; // x direction
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFaceNormals()
{
  CFAUTOTRACE;

  m_faceNormals.resize(5);
  for (CFuint iFace = 0; iFace < 5; ++iFace)
  {
    m_faceNormals[iFace].resize(3);
  }

  m_faceNormals[0][0] = 0.;
  m_faceNormals[0][1] = 0.;
  m_faceNormals[0][2] = -1.;

  m_faceNormals[1][0] = 0.;
  m_faceNormals[1][1] = 0.;
  m_faceNormals[1][2] = 1.;

  m_faceNormals[2][0] = 0.;
  m_faceNormals[2][1] = -1.;
  m_faceNormals[2][2] = 0.;

  m_faceNormals[3][0] = sqrt(2.)/2.;
  m_faceNormals[3][1] = sqrt(2.)/2.;
  m_faceNormals[3][2] = 0.;

  m_faceNormals[4][0] = -1.;
  m_faceNormals[4][1] = 0.;
  m_faceNormals[4][2] = 0.;

}

//////////////////////////////////////////////////////////////////////
// !!!! to be checked

void PrismFluxReconstructionElementData::createFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;


   // number of faces
  const CFuint nbrFaces = 5;

  const CFuint nbrQuadFaces = 3;

  const CFuint nbrTriagFaces = 2;

  // number of possible orientations
  const CFuint nbrOrient = 33;
  const CFuint nbrQuadOrient = 24;
  const CFuint nbrTriagOrient = 9;



  // resize the variables
  m_faceNodeConnPerOrient.resize(nbrOrient);
  m_faceConnPerOrient.resize(nbrOrient);
  m_faceMappedCoordDirPerOrient.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_faceNodeConnPerOrient[iOrient].resize(2);
    m_faceConnPerOrient[iOrient].resize(2);
    m_faceMappedCoordDirPerOrient[iOrient].resize(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_faceNodeConnPerOrient[iOrient][iSide].resize(3);
    }
  }

  // fill the variable
  CFuint iOrient = 0;

  for (CFuint iFaceL = 0; iFaceL < nbrTriagFaces; ++iFaceL)
  {
    vector< CFuint > faceNodesL = m_faceNodeConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrTriagFaces; ++iFaceR)
    {
      vector< CFuint > faceNodesR(3);
      faceNodesR[0] = m_faceNodeConn[iFaceR][1];
      faceNodesR[1] = m_faceNodeConn[iFaceR][0];
      faceNodesR[2] = m_faceNodeConn[iFaceR][2];

      for (CFuint iRot = 0; iRot < 3; ++iRot, ++iOrient)
      {
        m_faceConnPerOrient[iOrient][LEFT ] = iFaceL;
        m_faceConnPerOrient[iOrient][RIGHT] = iFaceR;

        m_faceMappedCoordDirPerOrient[iOrient][LEFT ] = m_faceMappedCoordDir[iFaceL];
        m_faceMappedCoordDirPerOrient[iOrient][RIGHT] = -m_faceMappedCoordDir[iFaceR];

        for (CFuint iNode = 0; iNode < 3; ++iNode)
        {
          m_faceNodeConnPerOrient[iOrient][LEFT ][iNode] = faceNodesL[iNode];
          m_faceNodeConnPerOrient[iOrient][RIGHT][iNode] = faceNodesR[iNode];
        }

        // rotate nodes of right face to new orientation
        CFuint swap = faceNodesR[2];
        faceNodesR[2] = faceNodesR[1];
        faceNodesR[1] = faceNodesR[0];
        faceNodesR[0] = swap;
      }
    }
  }

  for (CFuint iFaceL = 0+2; iFaceL < nbrQuadFaces+2; ++iFaceL)
  {
    vector< CFuint > faceNodesL = m_faceNodeConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrQuadFaces+2; ++iFaceR)
    {
      vector< CFuint > faceNodesR(4);
      faceNodesR[0] = m_faceNodeConn[iFaceR][1];
      faceNodesR[1] = m_faceNodeConn[iFaceR][0];
      faceNodesR[2] = m_faceNodeConn[iFaceR][3];
      faceNodesR[3] = m_faceNodeConn[iFaceR][2];

      for (CFuint iRot = 0; iRot < 4; ++iRot, ++iOrient)
      {
            for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_faceNodeConnPerOrient[iOrient][iSide].resize(4);
    }

        m_faceConnPerOrient[iOrient][LEFT ] = iFaceL;
        m_faceConnPerOrient[iOrient][RIGHT] = iFaceR;

        m_faceMappedCoordDirPerOrient[iOrient][LEFT ] = m_faceMappedCoordDir[iFaceL];
        m_faceMappedCoordDirPerOrient[iOrient][RIGHT] = -m_faceMappedCoordDir[iFaceR];

        for (CFuint iNode = 0; iNode < 4; ++iNode)
        {
          m_faceNodeConnPerOrient[iOrient][LEFT ][iNode] = faceNodesL[iNode];
          m_faceNodeConnPerOrient[iOrient][RIGHT][iNode] = faceNodesR[iNode];
        }

        // rotate nodes of right face to new orientation
        CFuint swap = faceNodesR[3];
        faceNodesR[3] = faceNodesR[2];
        faceNodesR[2] = faceNodesR[1];
        faceNodesR[1] = faceNodesR[0];
        faceNodesR[0] = swap;
      }

      cf_assert(faceNodesR[0] == m_faceNodeConn[iFaceR][1]);
      cf_assert(faceNodesR[1] == m_faceNodeConn[iFaceR][0]);
      cf_assert(faceNodesR[2] == m_faceNodeConn[iFaceR][3]);
      cf_assert(faceNodesR[3] == m_faceNodeConn[iFaceR][2]);
    }
  }
  
//cout << "createFaceNodeConnectivityPerOrient -- Nb Orientation=  "<<iOrient<<endl;

  cf_assert(iOrient == nbrOrient);

}


//////////////////////////////////////////////////////////////////////
// !!!! to be checked

void PrismFluxReconstructionElementData::createFaceIntegrationCoefs()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  // number of flux points on a face
  const CFuint nbrFlxPnts = nbrFlxPnts1D*nbrFlxPnts1D;

  // resize m_faceIntegrationCoefs
  m_faceIntegrationCoefs.resize(nbrFlxPnts);

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_2D,m_polyOrder);

  // create face node local coordinates
  vector< RealVector > nodeCoord(4);
  nodeCoord[0].resize(2);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[0][ETA] = -1.0;
  nodeCoord[1].resize(2);
  nodeCoord[1][KSI] = +1.0;
  nodeCoord[1][ETA] = -1.0;
  nodeCoord[2].resize(2);
  nodeCoord[2][KSI] = +1.0;
  nodeCoord[2][ETA] = +1.0;
  nodeCoord[3].resize(2);
  nodeCoord[3][KSI] = -1.0;
  nodeCoord[3][ETA] = +1.0;

  // get quadrature point coordinates and wheights
  vector< RealVector > quadPntCoords   = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  vector< CFreal     > quadPntWheights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWheights.size() == nbrQPnts);

  // compute the coefficients for integration over a face
  // loop over flux points
  CFuint iFlx = 0;
  for (CFuint iFlxKsi = 0; iFlxKsi < nbrFlxPnts1D; ++iFlxKsi)
  {
    const CFreal ksiFlx = m_flxPntsLocalCoord1D[iFlxKsi];
    for (CFuint iFlxEta = 0; iFlxEta < nbrFlxPnts1D; ++iFlxEta, ++iFlx)
    {
      const CFreal etaFlx = m_flxPntsLocalCoord1D[iFlxEta];

      m_faceIntegrationCoefs[iFlx] = 0.0;
      for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
      {
        // quadrature point local coordinate on the face
        const CFreal ksiQPnt = quadPntCoords[iQPnt][KSI];
        const CFreal etaQPnt = quadPntCoords[iQPnt][ETA];

        // evaluate polynomial value in quadrature point
        CFreal quadPntPolyVal = 1.;
        for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
        {
          if (iFac != iFlxKsi)
          {
            const CFreal ksiFac = m_flxPntsLocalCoord1D[iFac];
            quadPntPolyVal *= (ksiQPnt-ksiFac)/(ksiFlx-ksiFac);
          }
        }
        for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
        {
          if (iFac != iFlxEta)
          {
            const CFreal etaFac = m_flxPntsLocalCoord1D[iFac];
            quadPntPolyVal *= (etaQPnt-etaFac)/(etaFlx-etaFac);
          }
        }

        // add contribution of quadrature point to integration coefficient
        m_faceIntegrationCoefs[iFlx] += quadPntWheights[iQPnt]*quadPntPolyVal;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// !!!! to be checked

void PrismFluxReconstructionElementData::createCellAvgSolCoefs()
{
  CFAUTOTRACE;

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // resize m_cellAvgSolCoefs
  m_cellAvgSolCoefs.resize(nbrSolPnts);

  const CFuint nbrSolPnts2D = (m_polyOrder+1)*(m_polyOrder+2)/2;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();

  RealVector oneDcoeffs;
  RealVector twoDcoeffs;

  // resize oneDcoeffs
  oneDcoeffs.resize(nbrFlxPnts1D);

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_1D,m_polyOrder);

  // create face node local coordinates
  vector< RealVector > nodeCoord(2);
  nodeCoord[0].resize(1);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[1].resize(1);
  nodeCoord[1][KSI] = +1.0;

  // get quadrature point coordinates and wheights
  vector< RealVector > quadPntCoords   = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  vector< CFreal     > quadPntWheights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWheights.size() == nbrQPnts);

  // compute the coefficients for integration over a face
  // loop over flux points
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts1D; ++iFlx)
  {
    oneDcoeffs[iFlx] = 0.0;

    const CFreal ksiFlx = m_flxPntsLocalCoord1D[iFlx];
    for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
    {
      // quadrature point local coordinate on the face
      const CFreal ksiQPnt = quadPntCoords[iQPnt][KSI];

      // evaluate polynomial value in quadrature point
      CFreal quadPntPolyVal = 1.;
      for (CFuint iFac = 0; iFac < nbrFlxPnts1D; ++iFac)
      {
        if (iFac != iFlx)
        {
          const CFreal ksiFac = m_flxPntsLocalCoord1D[iFac];
          quadPntPolyVal *= (ksiQPnt-ksiFac)/(ksiFlx-ksiFac);
        }
      }

      // add contribution of quadrature point to integration coefficient
      oneDcoeffs[iFlx] += quadPntWheights[iQPnt]*quadPntPolyVal;
    }
    
  }


  m_cellAvgSolCoefs.resize(nbrSolPnts);
  twoDcoeffs.resize(nbrSolPnts2D);


  switch(m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      m_cellAvgSolCoefs[0] = 1.0;
    } break;
    case CFPolyOrder::ORDER1:
    {
      twoDcoeffs[0] = 0.166666666666667;
      twoDcoeffs[1] = 0.166666666666667;
      twoDcoeffs[2] = 0.166666666666667;

      for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
      {
        for (CFuint iFlx = 0; iFlx < nbrSolPnts2D; ++iFlx)
        {
          m_cellAvgSolCoefs[iZta*nbrSolPnts2D+iFlx]=(twoDcoeffs[iFlx]*2.)*oneDcoeffs[iZta];  // Sum should be = 1 so twoDcoeffs * 2
        }  
      }

      } break;

    case CFPolyOrder::ORDER2:
    {
      twoDcoeffs[0] = 0.054975871827667;
      twoDcoeffs[1] = 0.111690794839000;
      twoDcoeffs[2] = 0.054975871827667;
      twoDcoeffs[3] = 0.111690794839000;
      twoDcoeffs[4] = 0.111690794839000;
      twoDcoeffs[5] = 0.054975871827667;

      for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
      {
        for (CFuint iFlx = 0; iFlx < nbrSolPnts2D; ++iFlx)
        {
          m_cellAvgSolCoefs[iZta*nbrSolPnts2D+iFlx]=(twoDcoeffs[iFlx]*2.)*oneDcoeffs[iZta];
        }  
      }
    } break;

    case CFPolyOrder::ORDER3:
    {
      twoDcoeffs[0] = 0.020977756498325;
      twoDcoeffs[1] = 0.056049206035444;
      twoDcoeffs[2] = 0.056049206035444;
      twoDcoeffs[3] = 0.020977756498325;
      twoDcoeffs[4] = 0.056049206035444;
      twoDcoeffs[5] = 0.100771494292365;
      twoDcoeffs[6] = 0.056049206035444;
      twoDcoeffs[7] = 0.056049206035444;
      twoDcoeffs[8] = 0.056049206035444;
      twoDcoeffs[9] = 0.020977756498325;

      for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
      {
        for (CFuint iFlx = 0; iFlx < nbrSolPnts2D; ++iFlx)
        {
          m_cellAvgSolCoefs[iZta*nbrSolPnts2D+iFlx]=(twoDcoeffs[iFlx]*2.)*oneDcoeffs[iZta];
        }  
      }
    } break;

    case CFPolyOrder::ORDER4:
    {
      twoDcoeffs[0] = 0.008957727506152;
      twoDcoeffs[1] = 0.027874905013558;
      twoDcoeffs[2] = 0.038103031192768;
      twoDcoeffs[3] = 0.055749810027115;
      twoDcoeffs[4] = 0.063856097940632;
      twoDcoeffs[5] = 0.027874905013558;
      twoDcoeffs[6] = 0.063856097940632;
      twoDcoeffs[7] = 0.127712195881265;
      twoDcoeffs[8] = 0.063856097940632;
      twoDcoeffs[9] = 0.038103031192768;
      twoDcoeffs[10] = 0.038103031192768;
      twoDcoeffs[11] = 0.055749810027115;
      twoDcoeffs[12] = 0.027874905013558;
      twoDcoeffs[13] = 0.055749810027115;
      twoDcoeffs[14] = 0.008957727506152;

      for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
      {
        for (CFuint iFlx = 0; iFlx < nbrSolPnts2D; ++iFlx)
        {
          m_cellAvgSolCoefs[iZta*nbrSolPnts2D+iFlx]=(twoDcoeffs[iFlx]*2.)*oneDcoeffs[iZta];
        }  
      }
    } break;

    case CFPolyOrder::ORDER5:
    {
      twoDcoeffs[0] = 0.005179687348269;
      twoDcoeffs[1] = 0.014484634686237;
      twoDcoeffs[2] = 0.023023183297967;
      twoDcoeffs[3] = 0.023023183297967;
      twoDcoeffs[4] = 0.028969269372473;
      twoDcoeffs[5] = 0.014484634686237;
      twoDcoeffs[6] = 0.037697442163369;
      twoDcoeffs[7] = 0.023023183297967;
      twoDcoeffs[8] = 0.037697442163369;
      twoDcoeffs[9] = 0.075394884326738;
      twoDcoeffs[10] = 0.048773901186621;
      twoDcoeffs[11] = 0.037697442163369;
      twoDcoeffs[12] = 0.014484634686237;
      twoDcoeffs[13] = 0.028969269372473;
      twoDcoeffs[14] = 0.028969269372473;
      twoDcoeffs[15] = 0.046046366595935;
      twoDcoeffs[16] = 0.023023183297967;
      twoDcoeffs[17] = 0.037697442163369;
      twoDcoeffs[18] = 0.014484634686237;
      twoDcoeffs[19] = 0.075394884326738;
      twoDcoeffs[20] = 0.005179687348269;

      for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
      {
        for (CFuint iFlx = 0; iFlx < nbrSolPnts2D; ++iFlx)
        {
          m_cellAvgSolCoefs[iZta*nbrSolPnts2D+iFlx]=(twoDcoeffs[iFlx]*2.)*oneDcoeffs[iZta];
        }  
      }
    } break;  
    default:
    {
      throw Common::NotImplementedException (FromHere(),"m_cellAvgSolCoefs not implemented for order "
                                    + StringOps::to_str(m_polyOrder) + ".");
    }
  }
}

//////////////////////////////////////////////////////////////////////
// !!!!!!! maybe not used
void PrismFluxReconstructionElementData::createCellCenterDerivCoefs()
{
  CFAUTOTRACE;

  // center coordinate
  vector< RealVector > centerCoord(1,RealVector(3));
  centerCoord[0][KSI] = 0.0;
  centerCoord[0][ETA] = 0.0;
  centerCoord[0][ZTA] = 0.0;

  vector< vector< vector< CFreal > > > polyDerivs =
                                            getSolPolyDerivsAtNode(centerCoord);

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // set polynomial derivatives
  m_cellCenterDerivCoefs.resize(3);
  m_cellCenterDerivCoefs[KSI].resize(nbrSolPnts);
  m_cellCenterDerivCoefs[ETA].resize(nbrSolPnts);
  m_cellCenterDerivCoefs[ZTA].resize(nbrSolPnts);
  for (CFuint iPoly = 0; iPoly < nbrSolPnts; ++iPoly)
  {
    m_cellCenterDerivCoefs[KSI][iPoly] = polyDerivs[0][KSI][iPoly];
    m_cellCenterDerivCoefs[ETA][iPoly] = polyDerivs[0][ETA][iPoly];
    m_cellCenterDerivCoefs[ZTA][iPoly] = polyDerivs[0][ZTA][iPoly];
  }
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
                                                        vector< RealVector >& nodalSet)
{
  
  const CFuint nbrPnts1D = order+1;
  const CFuint nbrSolPnts2D = (order+1)*(order+2)/2;

  // set solution point local coordinates
  nodalSet.resize(0);
  for (CFuint iZta = 0; iZta < nbrPnts1D; ++iZta)
  {
    for (CFuint iFlx = 0; iFlx < nbrSolPnts2D; ++iFlx)
    {
      RealVector node(3);
      node[KSI] = solPntsLocalCoord2D[iFlx][KSI];
      node[ETA] = solPntsLocalCoord2D[iFlx][ETA];
      node[ZTA] = m_solPntsLocalCoord1D[iZta];
      nodalSet.push_back(node);
    }  
  }
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::setCFLConvDiffRatio()
{
  CFAUTOTRACE;

  switch(m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      m_cflConvDiffRatio = 2.0;//4.0; //2.0; // check this!  2/2
    } break;
    case CFPolyOrder::ORDER1:
    {
      m_cflConvDiffRatio = 4.0;//4.0;//6.5; //100.0; //200.0; // check this! 4/6
    } break;
    case CFPolyOrder::ORDER2:
    {
      m_cflConvDiffRatio = 6.0;//6.0;//17.0; //6.0; // check this! 6/10
    } break;
    case CFPolyOrder::ORDER3:
    {
      m_cflConvDiffRatio = 8.0;//8.0;//25.0; //500.0; //8.0; // check this! 8/14
    } break;
    case CFPolyOrder::ORDER4:
    {
      m_cflConvDiffRatio = 10.0;//10.0;//50.0; //10.0; // check this! 10/18
    } break;
    case CFPolyOrder::ORDER5:
    {
      m_cflConvDiffRatio = 12.0;//22.0;//12.0;//100.0; //1200.0; // check this! 12/22
    } break;
    case CFPolyOrder::ORDER6:
    {
      m_cflConvDiffRatio = 14.0;//14.0;//200.0; //14.0; // check this! 14/26
    } break;
    case CFPolyOrder::ORDER7:
    {
      m_cflConvDiffRatio = 16.0;//16.0;//2500.0; //16.0; // check this! 16/30
    } break;
    case CFPolyOrder::ORDER8:
    {
      m_cflConvDiffRatio = 18.0;//18.0;//800.0; //18.0; // check this! 18/34
    } break;
    case CFPolyOrder::ORDER9:
    {
      m_cflConvDiffRatio = 20.0;//20.0;//5500.0; //20.0; // check this! 20/38
    } break;
    case CFPolyOrder::ORDER10:
    {
      m_cflConvDiffRatio = 22.0;//22.0;//3200.0; //22.0; // check this! 22/42
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Higher-order quadrilateral FR cell not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////
// !!!! to be checked

void PrismFluxReconstructionElementData::createFaceOutputPntCellMappedCoords()
{
  // number of points on a face
  const CFuint nbrFacePnts1D = m_polyOrder == 0 ? 2 :m_polyOrder + 1;

  // face mapped coordinates of uniform distribution of points
  vector<RealVector> faceMapCoords;
  const CFreal dKsiEta = 0 ? 2.0 : 2.0/m_polyOrder;
  CFreal ksi = -1.0;
  m_faceOutputPntFaceMappedCoords.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrFacePnts1D; ++iKsi, ksi += dKsiEta)
  {
    CFreal eta = -1.0;
    for (CFuint iEta = 0; iEta < nbrFacePnts1D; ++iEta, eta += dKsiEta)
    {
      RealVector mapCoord(2);
      mapCoord[KSI] = ksi;
      mapCoord[ETA] = eta;
      m_faceOutputPntFaceMappedCoords.push_back(mapCoord);
    }
  }
  const CFuint nbrFacePnts = m_faceOutputPntFaceMappedCoords.size();
  cf_assert(nbrFacePnts == nbrFacePnts1D*nbrFacePnts1D);
  // compute cell mapped coordinates for distribution on each face
  const CFuint nbrCellFaces = getNbrCellFaces();
  m_faceOutputPntCellMappedCoords.resize(nbrCellFaces);

  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    // current face node coordinates
    const vector<RealVector>& faceNodeCoords = m_faceNodeCoords[iFace];
    m_faceOutputPntCellMappedCoords[iFace].resize(0);
    for (CFuint iPnt = 0; iPnt < nbrFacePnts; ++iPnt)
    {
      const CFreal fun0 = 0.25*(1.0-m_faceOutputPntFaceMappedCoords[iPnt][KSI])
                              *(1.0-m_faceOutputPntFaceMappedCoords[iPnt][ETA]);
      const CFreal fun1 = 0.25*(1.0+m_faceOutputPntFaceMappedCoords[iPnt][KSI])
                              *(1.0-m_faceOutputPntFaceMappedCoords[iPnt][ETA]);
      const CFreal fun2 = 0.25*(1.0+m_faceOutputPntFaceMappedCoords[iPnt][KSI])
                              *(1.0+m_faceOutputPntFaceMappedCoords[iPnt][ETA]);
      const CFreal fun3 = 0.25*(1.0-m_faceOutputPntFaceMappedCoords[iPnt][KSI])
                              *(1.0+m_faceOutputPntFaceMappedCoords[iPnt][ETA]);
      if (iFace==0 || iFace==1){
      m_faceOutputPntCellMappedCoords[iFace].push_back(fun0*faceNodeCoords[0]+
                                                       fun1*faceNodeCoords[1]+
                                                       fun2*faceNodeCoords[2]);
      }
      else{
      m_faceOutputPntCellMappedCoords[iFace].push_back(fun0*faceNodeCoords[0]+
                                                       fun1*faceNodeCoords[1]+
                                                       fun2*faceNodeCoords[2]+
                                                       fun3*faceNodeCoords[3]);
      }

    }
  }
}

//////////////////////////////////////////////////////////////////////
// !!!! to be checked
void PrismFluxReconstructionElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;

  // number of nodes 1D
  const CFuint nbrNodes1D = m_polyOrder + 1;

  m_faceOutputPntConn.resize(0);
  for (CFuint iKsi = 0; iKsi < static_cast<CFuint>(m_polyOrder); ++iKsi)
  {
    for (CFuint iEta = 0; iEta < static_cast<CFuint>(m_polyOrder); ++iEta)
    {
      vector< CFuint > cellNodesConn(4);
      cellNodesConn[0] = (iKsi  )*nbrNodes1D + iEta  ;
      cellNodesConn[1] = (iKsi+1)*nbrNodes1D + iEta  ;
      cellNodesConn[2] = (iKsi+1)*nbrNodes1D + iEta+1;
      cellNodesConn[3] = (iKsi  )*nbrNodes1D + iEta+1;
      m_faceOutputPntConn.push_back(cellNodesConn);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createVandermondeMatrix()
{
  CFAUTOTRACE;
  
  const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  
  CFuint order;
  CFuint idx;
  
  m_vandermonde.resize(nbrSolPnts,nbrSolPnts);
  m_vandermondeInv.resize(nbrSolPnts,nbrSolPnts);
  
  if(m_polyOrder != CFPolyOrder::ORDER0)
  {
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      CFuint modalDof = 0;

      // Loop over the total order of polynomial terms for ξ and η
      for (CFuint totalOrderXY = 0; totalOrderXY <= m_polyOrder; ++totalOrderXY)
      {
        // Loop for the third dimension (Zta)
        for (CFuint iOrderZta = 0; iOrderZta <= m_polyOrder; ++iOrderZta)
        {
          // Within each total order, iterate over powers of ksi and eta in reverse for ksi
          for (CFuint iOrderKsi = totalOrderXY; iOrderKsi >= 0; --iOrderKsi)
          {
            CFuint iOrderEta = totalOrderXY - iOrderKsi;
            
            // Evaluate the polynomial basis functions at the solution point
            double a = ((2. * m_solPntsLocalCoords[iSol][KSI]) / (1. - m_solPntsLocalCoords[iSol][ETA])) - 1.;
            double b = 2. * m_solPntsLocalCoords[iSol][ETA] - 1.;
            double h1 = ComputeJacobi(iOrderKsi, 0., 0., a);
            double h2 = ComputeJacobi(iOrderEta, 2. * iOrderKsi + 1., 0., b);
            double h3 = evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA], iOrderZta); 
            
            // Compute the value of the polynomial term at the solution point and store in the Vandermonde matrix
            m_vandermonde(iSol, modalDof) = 0.25*sqrt(2.0) * h1 * h2 * h3 * pow((1. - b), iOrderKsi);

            // Increment the modal degree of freedom index
            modalDof += 1;
            // Make sure we don't go below zero
            if (iOrderKsi == 0) break;
          }
        }
      }
      cf_assert(modalDof == nbrSolPnts);
    }
    
    InvertMatrix(m_vandermonde,m_vandermondeInv);
  }
}

//////////////////////////////////////////////////////////////////////

void PrismFluxReconstructionElementData::createFlxSolDependencies()
{
  CFAUTOTRACE;

  const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  const CFuint nbrFlxPnts = (m_polyOrder+1)*(m_polyOrder+2)+3*(m_polyOrder+1)*(m_polyOrder+1); //m_flxPntsLocalCoords.size();
  const CFuint nbrSolPnts1D = m_polyOrder+1;
  const CFuint nbrSolPnts2D = (m_polyOrder+1)*(m_polyOrder+2)/2;


  m_solFlxDep.resize(nbrSolPnts);
  m_solSolDep.resize(nbrSolPnts);  // also push back the current 
  m_flxSolDep.resize(nbrFlxPnts);

  
  /////////// Flx-sol dep

  for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
  {
    for (CFuint iFlx = 0; iFlx < nbrSolPnts2D; ++iFlx)
    {
      m_flxSolDep[iFlx].push_back(iFlx+nbrSolPnts2D*iZta);
      m_flxSolDep[iFlx+nbrSolPnts2D].push_back(iFlx+nbrSolPnts2D*iZta);
    }
  }
  CFuint iFlx = 2*nbrSolPnts2D;
  for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
  {
    for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi, ++iFlx)
    {
      for (CFuint iSol = 0; iSol < nbrSolPnts2D; ++iSol)
      {
        m_flxSolDep[iFlx].push_back(iSol+iZta*nbrSolPnts2D);
        m_flxSolDep[iFlx+nbrSolPnts1D*nbrSolPnts1D].push_back(iSol+iZta*nbrSolPnts2D);
        m_flxSolDep[iFlx+2*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol+iZta*nbrSolPnts2D);
      }
    }
  }

  
  /////////// Sol-flx dep
  
  CFuint iSol = 0;
  for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts2D; ++iEta, ++iSol)
    {        
      m_solFlxDep[iSol].push_back(iEta);
      m_solFlxDep[iSol].push_back(iEta+nbrSolPnts2D);
      for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
      {
        m_solFlxDep[iSol].push_back(nbrSolPnts2D*2+iZta*nbrSolPnts1D+iKsi);
        m_solFlxDep[iSol].push_back(nbrSolPnts2D*2+iZta*nbrSolPnts1D+nbrSolPnts1D*nbrSolPnts1D+iKsi);
        m_solFlxDep[iSol].push_back(nbrSolPnts2D*2+iZta*nbrSolPnts1D+2*nbrSolPnts1D*nbrSolPnts1D+iKsi);
      }
    }
  }

  /////////// Sol-sol dep

  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    CFuint iSol2d = iSol%nbrSolPnts2D;
    CFuint iSol1d = iSol/nbrSolPnts2D;

    for(CFuint iSol1= nbrSolPnts2D*iSol1d; iSol1< nbrSolPnts2D*(iSol1d+1); ++iSol1)
    {
      if (iSol1!=iSol)
      {
        m_solSolDep[iSol].push_back(iSol1);
      }
    }
    for(CFuint iSol1= 0; iSol1< nbrSolPnts1D; ++iSol1)
    {
      m_solSolDep[iSol].push_back(iSol2d+iSol1*nbrSolPnts2D);
    }
  }


}


//////////////////////////////////////////////////////////////////////////////
CFreal PrismFluxReconstructionElementData::ComputeJacobi(CFuint N,CFuint alpha, CFuint beta,CFreal x)
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
CFreal PrismFluxReconstructionElementData::factorial(CFreal n)
{
  return (n==1. || n==0.) ? 1. : factorial(n-1.)*n;
}
//////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
