#include <fstream>

#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "MathTools/MathConsts.hh"

#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////

TetraFluxReconstructionElementData::TetraFluxReconstructionElementData() :
  FluxReconstructionElementData()
{
  m_shape = CFGeoShape::TETRA;
  m_dimensionality = DIM_3D;
}

//////////////////////////////////////////////////////////////////////

TetraFluxReconstructionElementData::TetraFluxReconstructionElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::TETRA;
  m_dimensionality = DIM_3D;
  m_polyOrder = polyOrder;
  
  m_solPntsLocalCoord1D.resize(polyOrder+1);
  m_flxPntsLocalCoord1D.resize(polyOrder+1);

  // Use a default solution and flux point distribution: Gauss Legendre.  1D
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

//////////// 2D
//m_flxPntsLocalCoord2D.resize(((polyOrder+1)*(polyOrder+2))/2);

CFuint rows = (0.5*(polyOrder+1)*(polyOrder+2)) ;
m_flxPntsLocalCoord2D.resize(rows, std::vector< CFreal > (2));
switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
                        m_flxPntsLocalCoord2D[0][0] = 0.333333333333333; 
                        m_flxPntsLocalCoord2D[0][1] = 0.333333333333333;  
      } break;
      case CFPolyOrder::ORDER1:
      {
                        m_flxPntsLocalCoord2D[0][0] = 0.1666666666667; 
                        m_flxPntsLocalCoord2D[0][1] = 0.1666666666667;
                        m_flxPntsLocalCoord2D[1][0] = 0.6666666666667; 
                        m_flxPntsLocalCoord2D[1][1] = 0.1666666666667;  
                        m_flxPntsLocalCoord2D[2][0] = 0.1666666666667; 
                        m_flxPntsLocalCoord2D[2][1] = 0.6666666666667;

    } break;

      case CFPolyOrder::ORDER2:              
      {
                        m_flxPntsLocalCoord2D[0][0] =  0.091576213509780; 
                        m_flxPntsLocalCoord2D[0][1] =  0.091576213509780; 
                        
                        m_flxPntsLocalCoord2D[1][0] =  0.445948490915964; 
                        m_flxPntsLocalCoord2D[1][1] =  0.108103018168071; 

                        m_flxPntsLocalCoord2D[2][0] =  0.816847572980440; 
                        m_flxPntsLocalCoord2D[2][1] =  0.091576213509780;  

                        m_flxPntsLocalCoord2D[3][0] =  0.108103018168071; 
                        m_flxPntsLocalCoord2D[3][1] =  0.445948490915964;
                        
                        m_flxPntsLocalCoord2D[4][0] =  0.445948490915964; 
                        m_flxPntsLocalCoord2D[4][1] =  0.445948490915964;
                        
                        m_flxPntsLocalCoord2D[5][0] =  0.091576213509780; 
                        m_flxPntsLocalCoord2D[5][1] =  0.816847572980440;
                        
      } break;

    case CFPolyOrder::ORDER3:          
    {
                      m_flxPntsLocalCoord2D[0][0] = 0.055564052669793;
                      m_flxPntsLocalCoord2D[0][1] = 0.055564052669793;

                      m_flxPntsLocalCoord2D[1][0] = 0.295533711735893; 
                      m_flxPntsLocalCoord2D[1][1] = 0.070255540518384;  
                      
                      m_flxPntsLocalCoord2D[2][0] = 0.634210747745723; 
                      m_flxPntsLocalCoord2D[2][1] = 0.070255540518384;  	

                      m_flxPntsLocalCoord2D[3][0] = 0.888871894660413; 
                      m_flxPntsLocalCoord2D[3][1] = 0.055564052669793;	

                      m_flxPntsLocalCoord2D[4][0] = 0.070255540518384; 
                      m_flxPntsLocalCoord2D[4][1] = 0.295533711735893; 	

                      m_flxPntsLocalCoord2D[5][0] = 0.333333333333333; 
                      m_flxPntsLocalCoord2D[5][1] = 0.333333333333333;

                      m_flxPntsLocalCoord2D[6][0] = 0.634210747745723; 
                      m_flxPntsLocalCoord2D[6][1] = 0.295533711735893; 

                      m_flxPntsLocalCoord2D[7][0] = 0.070255540518384; 
                      m_flxPntsLocalCoord2D[7][1] = 0.634210747745723;	
                      
                      m_flxPntsLocalCoord2D[8][0] = 0.295533711735893; 
                      m_flxPntsLocalCoord2D[8][1] = 0.634210747745723; 
                      
                      m_flxPntsLocalCoord2D[9][0] = 0.055564052669793; 
                      m_flxPntsLocalCoord2D[9][1] = 0.888871894660413; 

      } break;
      case CFPolyOrder::ORDER4:   
      {
                        m_flxPntsLocalCoord2D[0][0] = 0.035870877695734; 
                        m_flxPntsLocalCoord2D[0][1] = 0.035870877695734;

                        m_flxPntsLocalCoord2D[1][0] = 0.201503881881800; 
                        m_flxPntsLocalCoord2D[1][1] =  0.047312487011716;

                        m_flxPntsLocalCoord2D[2][0] = 0.474308787777079;
                        m_flxPntsLocalCoord2D[2][1] = 0.051382424445843; 

                        m_flxPntsLocalCoord2D[3][0] = 0.751183631106484; 
                        m_flxPntsLocalCoord2D[3][1] = 0.047312487011716;

                        m_flxPntsLocalCoord2D[4][0] = 0.928258244608533; 
                        m_flxPntsLocalCoord2D[4][1] = 0.035870877695734;


                        m_flxPntsLocalCoord2D[5][0] = 0.047312487011716;
                        m_flxPntsLocalCoord2D[5][1] = 0.201503881881800;

                        m_flxPntsLocalCoord2D[6][0] = 0.241729395767967;
                        m_flxPntsLocalCoord2D[6][1] = 0.241729395767967;

                        m_flxPntsLocalCoord2D[7][0] = 0.516541208464066; 
                        m_flxPntsLocalCoord2D[7][1] = 0.241729395767967; 

                        m_flxPntsLocalCoord2D[8][0] = 0.751183631106484; 
                        m_flxPntsLocalCoord2D[8][1] = 0.201503881881800;

                        m_flxPntsLocalCoord2D[9][0] = 0.051382424445843;
                        m_flxPntsLocalCoord2D[9][1] = 0.474308787777079;

                        m_flxPntsLocalCoord2D[10][0] = 0.241729395767967; 
                        m_flxPntsLocalCoord2D[10][1] = 0.516541208464066; 

                        m_flxPntsLocalCoord2D[11][0] = 0.474308787777079; 
                        m_flxPntsLocalCoord2D[11][1] = 0.474308787777079;

                        m_flxPntsLocalCoord2D[12][0] = 0.047312487011716; 
                        m_flxPntsLocalCoord2D[12][1] = 0.751183631106484;

                        m_flxPntsLocalCoord2D[13][0] =  0.201503881881800;
                        m_flxPntsLocalCoord2D[13][1] =  0.751183631106484;

                        m_flxPntsLocalCoord2D[14][0] = 0.035870877695734; 
                        m_flxPntsLocalCoord2D[14][1] = 0.928258244608533; 

      } break;

      case CFPolyOrder::ORDER5:
      {
                        m_flxPntsLocalCoord2D[0][0] = 0.028112952182664; 
                        m_flxPntsLocalCoord2D[0][1] = 0.028112952182664; 
 
                        m_flxPntsLocalCoord2D[1][0] = 0.148565812270887; 
                        m_flxPntsLocalCoord2D[1][1] = 0.033533207700614; 

                        m_flxPntsLocalCoord2D[2][0] = 0.357196298615681; 
                        m_flxPntsLocalCoord2D[2][1] = 0.037824789609186;

                        m_flxPntsLocalCoord2D[3][0] = 0.604978911775132; 
                        m_flxPntsLocalCoord2D[3][1] = 0.037824789609186;

                        m_flxPntsLocalCoord2D[4][0] = 0.817900980028499; 
                        m_flxPntsLocalCoord2D[4][1] = 0.033533207700614;

                        m_flxPntsLocalCoord2D[5][0] = 0.943774095634672;
                        m_flxPntsLocalCoord2D[5][1] = 0.028112952182664;

                        m_flxPntsLocalCoord2D[6][0] = 0.033533207700614;
                        m_flxPntsLocalCoord2D[6][1] = 0.148565812270887;

                        m_flxPntsLocalCoord2D[7][0] = 0.177139098469317;
                        m_flxPntsLocalCoord2D[7][1] = 0.177139098469317;

                        m_flxPntsLocalCoord2D[8][0] = 0.405508595867433;
                        m_flxPntsLocalCoord2D[8][1] = 0.188982808265134;

                        m_flxPntsLocalCoord2D[9][0] = 0.645721803061365; 
                        m_flxPntsLocalCoord2D[9][1] = 0.177139098469317;

                        m_flxPntsLocalCoord2D[10][0] = 0.817900980028499;
                        m_flxPntsLocalCoord2D[10][1] = 0.148565812270887;


                        m_flxPntsLocalCoord2D[11][0] = 0.037824789609186; 
                        m_flxPntsLocalCoord2D[11][1] = 0.357196298615681; 

                        m_flxPntsLocalCoord2D[12][0] = 0.188982808265134;
                        m_flxPntsLocalCoord2D[12][1] = 0.405508595867433;

                        m_flxPntsLocalCoord2D[13][0] = 0.405508595867433;
                        m_flxPntsLocalCoord2D[13][1] = 0.405508595867433;

                        m_flxPntsLocalCoord2D[14][0] = 0.604978911775132; 
                        m_flxPntsLocalCoord2D[14][1] = 0.357196298615681; 


                        m_flxPntsLocalCoord2D[15][0] = 0.037824789609186; 
                        m_flxPntsLocalCoord2D[15][1] = 0.604978911775132;

                        m_flxPntsLocalCoord2D[16][0] = 0.177139098469317;
                        m_flxPntsLocalCoord2D[16][1] = 0.645721803061365;

                        m_flxPntsLocalCoord2D[17][0] = 0.357196298615681;
                        m_flxPntsLocalCoord2D[17][1] = 0.604978911775132;

                        m_flxPntsLocalCoord2D[18][0] = 0.033533207700614; 
                        m_flxPntsLocalCoord2D[18][1] = 0.817900980028499;

                        m_flxPntsLocalCoord2D[19][0] = 0.148565812270887;
                        m_flxPntsLocalCoord2D[19][1] = 0.817900980028499;

                        m_flxPntsLocalCoord2D[20][0] = 0.028112952182664;
                        m_flxPntsLocalCoord2D[20][1] = 0.943774095634672;

      } break;


   case CFPolyOrder::ORDER6:
      {  
                        m_flxPntsLocalCoord2D[0][0] =  0.0000000000000;
                        m_flxPntsLocalCoord2D[0][1] =  0.9451704450174; 
                        m_flxPntsLocalCoord2D[1][0] =  0.9451704450173; 
                        m_flxPntsLocalCoord2D[1][1] =  0.0000000000000; 
                        m_flxPntsLocalCoord2D[2][0] =  0.9289002405719;
                        m_flxPntsLocalCoord2D[2][1] =  0.0685505797224;
                        m_flxPntsLocalCoord2D[3][0] =  0.0685505797224;
                        m_flxPntsLocalCoord2D[3][1] =  0.9289002405717;
                        m_flxPntsLocalCoord2D[4][0] =  0.0243268355615;
                        m_flxPntsLocalCoord2D[4][1] =  0.0243268355616; 
                        m_flxPntsLocalCoord2D[5][0] =  0.1279662835335;
                        m_flxPntsLocalCoord2D[5][1] =  0.0277838749488;
                        m_flxPntsLocalCoord2D[6][0] =  0.0277838749488; 
                        m_flxPntsLocalCoord2D[6][1] =  0.1279662835337; 
                        m_flxPntsLocalCoord2D[7][0] =  0.0287083428360; 
                        m_flxPntsLocalCoord2D[7][1] =  0.7498347588657;
                        m_flxPntsLocalCoord2D[8][0] =  0.7498347588656;
                        m_flxPntsLocalCoord2D[8][1] =  0.0287083428360;
                        m_flxPntsLocalCoord2D[9][0] =  0.7228007909707;
                        m_flxPntsLocalCoord2D[9][1] =  0.2497602062385;
                        m_flxPntsLocalCoord2D[10][0] = 0.2497602062386; 
                        m_flxPntsLocalCoord2D[10][1] = 0.7228007909707;
                        m_flxPntsLocalCoord2D[11][0] = 0.0865562992839;
                        m_flxPntsLocalCoord2D[11][1] = 0.8325513856997;
                        m_flxPntsLocalCoord2D[12][0] = 0.8325513856998;
                        m_flxPntsLocalCoord2D[12][1] = 0.0865562992839;
                        m_flxPntsLocalCoord2D[13][0] = 0.3061619157672;
                        m_flxPntsLocalCoord2D[13][1] = 0.0303526617491;
                        m_flxPntsLocalCoord2D[14][0] = 0.0303526617491;
                        m_flxPntsLocalCoord2D[14][1] = 0.3061619157675;
                        m_flxPntsLocalCoord2D[15][0] = 0.4868610595047;
                        m_flxPntsLocalCoord2D[15][1] = 0.4868610595047;
                        m_flxPntsLocalCoord2D[16][0] = 0.6657904293017;
                        m_flxPntsLocalCoord2D[16][1] = 0.1765456154219;
                        m_flxPntsLocalCoord2D[17][0] = 0.1765456154221; 
                        m_flxPntsLocalCoord2D[17][1] = 0.6657904293016;
                        m_flxPntsLocalCoord2D[18][0] = 0.0293121007360; 
                        m_flxPntsLocalCoord2D[18][1] = 0.5295657488669;
                        m_flxPntsLocalCoord2D[19][0] = 0.5295657488667;
                        m_flxPntsLocalCoord2D[19][1] = 0.0293121007360;
                        m_flxPntsLocalCoord2D[20][0] = 0.1444673824391; 
                        m_flxPntsLocalCoord2D[20][1] = 0.1444673824391;
                        m_flxPntsLocalCoord2D[21][0] = 0.3299740111411; 
                        m_flxPntsLocalCoord2D[21][1] = 0.5361815729050;
                        m_flxPntsLocalCoord2D[22][0] = 0.5361815729052; 
                        m_flxPntsLocalCoord2D[22][1] = 0.3299740111409;
                        m_flxPntsLocalCoord2D[23][0] = 0.5511507516862;
                        m_flxPntsLocalCoord2D[23][1] = 0.1437790861923;
                        m_flxPntsLocalCoord2D[24][0] = 0.1437790861923;
                        m_flxPntsLocalCoord2D[24][1] = 0.5511507516862;
                        m_flxPntsLocalCoord2D[25][0] = 0.3348066587327; 
                        m_flxPntsLocalCoord2D[25][1] = 0.1529619437161;
                        m_flxPntsLocalCoord2D[26][0] = 0.1529619437161; 
                        m_flxPntsLocalCoord2D[26][1] = 0.3348066587327;
                        m_flxPntsLocalCoord2D[27][0] = 0.3430183498147; 
                        m_flxPntsLocalCoord2D[27][1] = 0.3430183498147;
                      
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Quadrature Points not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
    }

    ///////////////////////// 3D

    //m_solPntsLocalCoord3D.resize(((polyOrder+1)*(polyOrder+2)*(polyOrder+3))/6);
    
    CFuint col = (((polyOrder+1)*(polyOrder+2)*(polyOrder+3))/6) ;
    m_solPntsLocalCoord3D.resize(col, std::vector< CFreal > (3));

    switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
        m_solPntsLocalCoord3D[0][0] = 0.25; 
        m_solPntsLocalCoord3D[0][1] = 0.25;
        m_solPntsLocalCoord3D[0][2] = 0.25;  
      } break;
      case CFPolyOrder::ORDER1:
      {
        m_solPntsLocalCoord3D[0][0]=0.1382;
        m_solPntsLocalCoord3D[0][1]=0.1382;
        m_solPntsLocalCoord3D[0][2]=0.1382;
        m_solPntsLocalCoord3D[1][0]=0.58541;
        m_solPntsLocalCoord3D[1][1]=0.1382;
        m_solPntsLocalCoord3D[1][2]=0.1382;
        m_solPntsLocalCoord3D[2][0]=0.1382;
        m_solPntsLocalCoord3D[2][1]=0.58541;
        m_solPntsLocalCoord3D[2][2]=0.1382;
        m_solPntsLocalCoord3D[3][0]=0.1382;
        m_solPntsLocalCoord3D[3][1]=0.1382;
        m_solPntsLocalCoord3D[3][2]=0.58541;
      } break;
      case CFPolyOrder::ORDER2:              
      {
        m_solPntsLocalCoord3D[0][0]=0.073835;
        m_solPntsLocalCoord3D[0][1]=0.073835;
        m_solPntsLocalCoord3D[0][2]=0.073835;
        m_solPntsLocalCoord3D[1][0]=0.7785;
        m_solPntsLocalCoord3D[1][1]=0.073835;
        m_solPntsLocalCoord3D[1][2]=0.073835;
        m_solPntsLocalCoord3D[2][0]=0.073835;
        m_solPntsLocalCoord3D[2][1]=0.7785;
        m_solPntsLocalCoord3D[2][2]=0.073835;
        m_solPntsLocalCoord3D[3][0]=0.073835;
        m_solPntsLocalCoord3D[3][1]=0.073835;
        m_solPntsLocalCoord3D[3][2]=0.7785;
        m_solPntsLocalCoord3D[4][0]=0.40624;
        m_solPntsLocalCoord3D[4][1]=0.093756;
        m_solPntsLocalCoord3D[4][2]=0.093756;
        m_solPntsLocalCoord3D[5][0]=0.093756;
        m_solPntsLocalCoord3D[5][1]=0.40624;
        m_solPntsLocalCoord3D[5][2]=0.093756;
        m_solPntsLocalCoord3D[6][0]=0.093756;
        m_solPntsLocalCoord3D[6][1]=0.093756;
        m_solPntsLocalCoord3D[6][2]=0.40624;
        m_solPntsLocalCoord3D[7][0]=0.40624;
        m_solPntsLocalCoord3D[7][1]=0.40624;
        m_solPntsLocalCoord3D[7][2]=0.093756;
        m_solPntsLocalCoord3D[8][0]=0.40624;
        m_solPntsLocalCoord3D[8][1]=0.093756;
        m_solPntsLocalCoord3D[8][2]=0.40624;
        m_solPntsLocalCoord3D[9][0]=0.093756;
        m_solPntsLocalCoord3D[9][1]=0.40624;
        m_solPntsLocalCoord3D[9][2]=0.40624;
      } break;
      case CFPolyOrder::ORDER3:              
      {
        m_solPntsLocalCoord3D[0][0]=0.032353;
        m_solPntsLocalCoord3D[0][1]=0.032353;
        m_solPntsLocalCoord3D[0][2]=0.032353;
        m_solPntsLocalCoord3D[1][0]=0.90294;
        m_solPntsLocalCoord3D[1][1]=0.032353;
        m_solPntsLocalCoord3D[1][2]=0.032353;
        m_solPntsLocalCoord3D[2][0]=0.032353;
        m_solPntsLocalCoord3D[2][1]=0.90294;
        m_solPntsLocalCoord3D[2][2]=0.032353;
        m_solPntsLocalCoord3D[3][0]=0.032353;
        m_solPntsLocalCoord3D[3][1]=0.032353;
        m_solPntsLocalCoord3D[3][2]=0.90294;
        m_solPntsLocalCoord3D[4][0]=0.6166;
        m_solPntsLocalCoord3D[4][1]=0.06036;
        m_solPntsLocalCoord3D[4][2]=0.06036;
        m_solPntsLocalCoord3D[5][0]=0.26268;
        m_solPntsLocalCoord3D[5][1]=0.06036;
        m_solPntsLocalCoord3D[5][2]=0.06036;
        m_solPntsLocalCoord3D[6][0]=0.06036;
        m_solPntsLocalCoord3D[6][1]=0.6166;
        m_solPntsLocalCoord3D[6][2]=0.06036;
        m_solPntsLocalCoord3D[7][0]=0.06036;
        m_solPntsLocalCoord3D[7][1]=0.26268;
        m_solPntsLocalCoord3D[7][2]=0.06036;
        m_solPntsLocalCoord3D[8][0]=0.06036;
        m_solPntsLocalCoord3D[8][1]=0.06036;
        m_solPntsLocalCoord3D[8][2]=0.6166;
        m_solPntsLocalCoord3D[9][0]=0.06036;
        m_solPntsLocalCoord3D[9][1]=0.06036;
        m_solPntsLocalCoord3D[9][2]=0.26268;
        m_solPntsLocalCoord3D[10][0]=0.26268;
        m_solPntsLocalCoord3D[10][1]=0.6166;
        m_solPntsLocalCoord3D[10][2]=0.06036;
        m_solPntsLocalCoord3D[11][0]=0.6166;
        m_solPntsLocalCoord3D[11][1]=0.26268;
        m_solPntsLocalCoord3D[11][2]=0.06036;
        m_solPntsLocalCoord3D[12][0]=0.26268;
        m_solPntsLocalCoord3D[12][1]=0.06036;
        m_solPntsLocalCoord3D[12][2]=0.6166;
        m_solPntsLocalCoord3D[13][0]=0.6166;
        m_solPntsLocalCoord3D[13][1]=0.06036;
        m_solPntsLocalCoord3D[13][2]=0.26268;
        m_solPntsLocalCoord3D[14][0]=0.06036;
        m_solPntsLocalCoord3D[14][1]=0.26268;
        m_solPntsLocalCoord3D[14][2]=0.6166;
        m_solPntsLocalCoord3D[15][0]=0.06036;
        m_solPntsLocalCoord3D[15][1]=0.6166;
        m_solPntsLocalCoord3D[15][2]=0.26268;
        m_solPntsLocalCoord3D[16][0]=0.30977;
        m_solPntsLocalCoord3D[16][1]=0.30977;
        m_solPntsLocalCoord3D[16][2]=0.070692;
        m_solPntsLocalCoord3D[17][0]=0.30977;
        m_solPntsLocalCoord3D[17][1]=0.070692;
        m_solPntsLocalCoord3D[17][2]=0.30977;
        m_solPntsLocalCoord3D[18][0]=0.070692;
        m_solPntsLocalCoord3D[18][1]=0.30977;
        m_solPntsLocalCoord3D[18][2]=0.30977;
        m_solPntsLocalCoord3D[19][0]=0.30977;
        m_solPntsLocalCoord3D[19][1]=0.30977;
        m_solPntsLocalCoord3D[19][2]=0.30977;
      } break;
      case CFPolyOrder::ORDER4:              
      {
        m_solPntsLocalCoord3D[0][0]=0.25;
        m_solPntsLocalCoord3D[0][1]=0.25;
        m_solPntsLocalCoord3D[0][2]=0.25;
        m_solPntsLocalCoord3D[1][0]=0.0267367755543735;
        m_solPntsLocalCoord3D[1][1]=0.0267367755543735;
        m_solPntsLocalCoord3D[1][2]=0.91978967333688;
        m_solPntsLocalCoord3D[2][0]=0.91978967333688;
        m_solPntsLocalCoord3D[2][1]=0.0267367755543735;
        m_solPntsLocalCoord3D[2][2]=0.0267367755543735;
        m_solPntsLocalCoord3D[3][0]=0.0267367755543735;
        m_solPntsLocalCoord3D[3][1]=0.91978967333688;
        m_solPntsLocalCoord3D[3][2]=0.0267367755543735;
        m_solPntsLocalCoord3D[4][0]=0.0267367755543735;
        m_solPntsLocalCoord3D[4][1]=0.0267367755543735;
        m_solPntsLocalCoord3D[4][2]=0.0267367755543735;
        m_solPntsLocalCoord3D[5][0]=0.0452454000155175;
        m_solPntsLocalCoord3D[5][1]=0.454754599984483;
        m_solPntsLocalCoord3D[5][2]=0.0452454000155175;
        m_solPntsLocalCoord3D[6][0]=0.454754599984483;
        m_solPntsLocalCoord3D[6][1]=0.0452454000155175;
        m_solPntsLocalCoord3D[6][2]=0.0452454000155175;
        m_solPntsLocalCoord3D[7][0]=0.0452454000155175;
        m_solPntsLocalCoord3D[7][1]=0.454754599984483;
        m_solPntsLocalCoord3D[7][2]=0.454754599984483;
        m_solPntsLocalCoord3D[8][0]=0.454754599984483;
        m_solPntsLocalCoord3D[8][1]=0.0452454000155175;
        m_solPntsLocalCoord3D[8][2]=0.454754599984483;
        m_solPntsLocalCoord3D[9][0]=0.454754599984483;
        m_solPntsLocalCoord3D[9][1]=0.454754599984483;
        m_solPntsLocalCoord3D[9][2]=0.0452454000155175;
        m_solPntsLocalCoord3D[10][0]=0.0452454000155175;
        m_solPntsLocalCoord3D[10][1]=0.0452454000155175;
        m_solPntsLocalCoord3D[10][2]=0.454754599984483;
        m_solPntsLocalCoord3D[11][0]=0.747759888481808;
        m_solPntsLocalCoord3D[11][1]=0.039102240635649;
        m_solPntsLocalCoord3D[11][2]=0.174035630246894;
        m_solPntsLocalCoord3D[12][0]=0.747759888481808;
        m_solPntsLocalCoord3D[12][1]=0.039102240635649;
        m_solPntsLocalCoord3D[12][2]=0.039102240635649;
        m_solPntsLocalCoord3D[13][0]=0.039102240635649;
        m_solPntsLocalCoord3D[13][1]=0.039102240635649;
        m_solPntsLocalCoord3D[13][2]=0.747759888481808;
        m_solPntsLocalCoord3D[14][0]=0.174035630246894;
        m_solPntsLocalCoord3D[14][1]=0.747759888481808;
        m_solPntsLocalCoord3D[14][2]=0.039102240635649;
        m_solPntsLocalCoord3D[15][0]=0.039102240635649;
        m_solPntsLocalCoord3D[15][1]=0.174035630246894;
        m_solPntsLocalCoord3D[15][2]=0.747759888481808;
        m_solPntsLocalCoord3D[16][0]=0.039102240635649;
        m_solPntsLocalCoord3D[16][1]=0.747759888481808;
        m_solPntsLocalCoord3D[16][2]=0.039102240635649;
        m_solPntsLocalCoord3D[17][0]=0.174035630246894;
        m_solPntsLocalCoord3D[17][1]=0.039102240635649;
        m_solPntsLocalCoord3D[17][2]=0.747759888481808;
        m_solPntsLocalCoord3D[18][0]=0.039102240635649;
        m_solPntsLocalCoord3D[18][1]=0.174035630246894;
        m_solPntsLocalCoord3D[18][2]=0.039102240635649;
        m_solPntsLocalCoord3D[19][0]=0.039102240635649;
        m_solPntsLocalCoord3D[19][1]=0.039102240635649;
        m_solPntsLocalCoord3D[19][2]=0.174035630246894;
        m_solPntsLocalCoord3D[20][0]=0.039102240635649;
        m_solPntsLocalCoord3D[20][1]=0.747759888481808;
        m_solPntsLocalCoord3D[20][2]=0.174035630246894;
        m_solPntsLocalCoord3D[21][0]=0.174035630246894;
        m_solPntsLocalCoord3D[21][1]=0.039102240635649;
        m_solPntsLocalCoord3D[21][2]=0.039102240635649;
        m_solPntsLocalCoord3D[22][0]=0.747759888481808;
        m_solPntsLocalCoord3D[22][1]=0.174035630246894;
        m_solPntsLocalCoord3D[22][2]=0.039102240635649;
        m_solPntsLocalCoord3D[23][0]=0.503118645014598;
        m_solPntsLocalCoord3D[23][1]=0.223201037962315;
        m_solPntsLocalCoord3D[23][2]=0.050479279060772;
        m_solPntsLocalCoord3D[24][0]=0.503118645014598;
        m_solPntsLocalCoord3D[24][1]=0.223201037962315;
        m_solPntsLocalCoord3D[24][2]=0.223201037962315;
        m_solPntsLocalCoord3D[25][0]=0.223201037962315;
        m_solPntsLocalCoord3D[25][1]=0.223201037962315;
        m_solPntsLocalCoord3D[25][2]=0.503118645014598;
        m_solPntsLocalCoord3D[26][0]=0.050479279060772;
        m_solPntsLocalCoord3D[26][1]=0.503118645014598;
        m_solPntsLocalCoord3D[26][2]=0.223201037962315;
        m_solPntsLocalCoord3D[27][0]=0.223201037962315;
        m_solPntsLocalCoord3D[27][1]=0.050479279060772;
        m_solPntsLocalCoord3D[27][2]=0.503118645014598;
        m_solPntsLocalCoord3D[28][0]=0.223201037962315;
        m_solPntsLocalCoord3D[28][1]=0.503118645014598;
        m_solPntsLocalCoord3D[28][2]=0.223201037962315;
        m_solPntsLocalCoord3D[29][0]=0.050479279060772;
        m_solPntsLocalCoord3D[29][1]=0.223201037962315;
        m_solPntsLocalCoord3D[29][2]=0.503118645014598;
        m_solPntsLocalCoord3D[30][0]=0.223201037962315;
        m_solPntsLocalCoord3D[30][1]=0.050479279060772;
        m_solPntsLocalCoord3D[30][2]=0.223201037962315;
        m_solPntsLocalCoord3D[31][0]=0.223201037962315;
        m_solPntsLocalCoord3D[31][1]=0.223201037962315;
        m_solPntsLocalCoord3D[31][2]=0.050479279060772;
        m_solPntsLocalCoord3D[32][0]=0.223201037962315;
        m_solPntsLocalCoord3D[32][1]=0.503118645014598;
        m_solPntsLocalCoord3D[32][2]=0.050479279060772;
        m_solPntsLocalCoord3D[33][0]=0.050479279060772;
        m_solPntsLocalCoord3D[33][1]=0.223201037962315;
        m_solPntsLocalCoord3D[33][2]=0.223201037962315;
        m_solPntsLocalCoord3D[34][0]=0.503118645014598;
        m_solPntsLocalCoord3D[34][1]=0.050479279060772;
        m_solPntsLocalCoord3D[34][2]=0.223201037962315;
      } break;
      case CFPolyOrder::ORDER5:              
      {
                        m_solPntsLocalCoord3D[0][0]=0.175924928891649;
                        m_solPntsLocalCoord3D[0][1]=0.175924928891649;
                        m_solPntsLocalCoord3D[0][2]=0.472225213325055;
                        m_solPntsLocalCoord3D[1][0]=0.175924928891649;
                        m_solPntsLocalCoord3D[1][1]=0.472225213325055;
                        m_solPntsLocalCoord3D[1][2]=0.175924928891649;
                        m_solPntsLocalCoord3D[2][0]=0.472225213325055;
                        m_solPntsLocalCoord3D[2][1]=0.175924928891649;
                        m_solPntsLocalCoord3D[2][2]=0.175924928891649;
                        m_solPntsLocalCoord3D[3][0]=0.175924928891649;
                        m_solPntsLocalCoord3D[3][1]=0.175924928891649;
                        m_solPntsLocalCoord3D[3][2]=0.175924928891649;
                        m_solPntsLocalCoord3D[4][0]=0.214336612970414;
                        m_solPntsLocalCoord3D[4][1]=0.214336612970414;
                        m_solPntsLocalCoord3D[4][2]=0.356990161088759;
                        m_solPntsLocalCoord3D[5][0]=0.214336612970414;
                        m_solPntsLocalCoord3D[5][1]=0.356990161088759;
                        m_solPntsLocalCoord3D[5][2]=0.214336612970414;
                        m_solPntsLocalCoord3D[6][0]=0.356990161088759;
                        m_solPntsLocalCoord3D[6][1]=0.214336612970414;
                        m_solPntsLocalCoord3D[6][2]=0.214336612970414;
                        m_solPntsLocalCoord3D[7][0]=0.214336612970414;
                        m_solPntsLocalCoord3D[7][1]=0.214336612970414;
                        m_solPntsLocalCoord3D[7][2]=0.214336612970414;
                        m_solPntsLocalCoord3D[8][0]=0.174463680162906;
                        m_solPntsLocalCoord3D[8][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[8][2]=0.0400417579050455;
                        m_solPntsLocalCoord3D[9][0]=0.174463680162906;
                        m_solPntsLocalCoord3D[9][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[9][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[10][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[10][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[10][2]=0.174463680162906;
                        m_solPntsLocalCoord3D[11][0]=0.0400417579050455;
                        m_solPntsLocalCoord3D[11][1]=0.174463680162906;
                        m_solPntsLocalCoord3D[11][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[12][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[12][1]=0.0400417579050455;
                        m_solPntsLocalCoord3D[12][2]=0.174463680162906;
                        m_solPntsLocalCoord3D[13][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[13][1]=0.174463680162906;
                        m_solPntsLocalCoord3D[13][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[14][0]=0.0400417579050455;
                        m_solPntsLocalCoord3D[14][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[14][2]=0.174463680162906;
                        m_solPntsLocalCoord3D[15][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[15][1]=0.0400417579050455;
                        m_solPntsLocalCoord3D[15][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[16][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[16][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[16][2]=0.0400417579050455;
                        m_solPntsLocalCoord3D[17][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[17][1]=0.174463680162906;
                        m_solPntsLocalCoord3D[17][2]=0.0400417579050455;
                        m_solPntsLocalCoord3D[18][0]=0.0400417579050455;
                        m_solPntsLocalCoord3D[18][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[18][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[19][0]=0.174463680162906;
                        m_solPntsLocalCoord3D[19][1]=0.0400417579050455;
                        m_solPntsLocalCoord3D[19][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[20][0]=0.611348339064507;
                        m_solPntsLocalCoord3D[20][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[20][2]=0.325129236353862;
                        m_solPntsLocalCoord3D[21][0]=0.611348339064507;
                        m_solPntsLocalCoord3D[21][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[21][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[22][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[22][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[22][2]=0.611348339064507;
                        m_solPntsLocalCoord3D[23][0]=0.325129236353862;
                        m_solPntsLocalCoord3D[23][1]=0.611348339064507;
                        m_solPntsLocalCoord3D[23][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[24][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[24][1]=0.325129236353862;
                        m_solPntsLocalCoord3D[24][2]=0.611348339064507;
                        m_solPntsLocalCoord3D[25][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[25][1]=0.611348339064507;
                        m_solPntsLocalCoord3D[25][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[26][0]=0.325129236353862;
                        m_solPntsLocalCoord3D[26][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[26][2]=0.611348339064507;
                        m_solPntsLocalCoord3D[27][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[27][1]=0.325129236353862;
                        m_solPntsLocalCoord3D[27][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[28][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[28][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[28][2]=0.325129236353862;
                        m_solPntsLocalCoord3D[29][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[29][1]=0.611348339064507;
                        m_solPntsLocalCoord3D[29][2]=0.325129236353862;
                        m_solPntsLocalCoord3D[30][0]=0.325129236353862;
                        m_solPntsLocalCoord3D[30][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[30][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[31][0]=0.611348339064507;
                        m_solPntsLocalCoord3D[31][1]=0.325129236353862;
                        m_solPntsLocalCoord3D[31][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[32][0]=0.856669692468503;
                        m_solPntsLocalCoord3D[32][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[32][2]=0.098311600856684;
                        m_solPntsLocalCoord3D[33][0]=0.856669692468503;
                        m_solPntsLocalCoord3D[33][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[33][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[34][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[34][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[34][2]=0.856669692468503;
                        m_solPntsLocalCoord3D[35][0]=0.098311600856684;
                        m_solPntsLocalCoord3D[35][1]=0.856669692468503;
                        m_solPntsLocalCoord3D[35][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[36][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[36][1]=0.098311600856684;
                        m_solPntsLocalCoord3D[36][2]=0.856669692468503;
                        m_solPntsLocalCoord3D[37][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[37][1]=0.856669692468503;
                        m_solPntsLocalCoord3D[37][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[38][0]=0.098311600856684;
                        m_solPntsLocalCoord3D[38][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[38][2]=0.856669692468503;
                        m_solPntsLocalCoord3D[39][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[39][1]=0.098311600856684;
                        m_solPntsLocalCoord3D[39][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[40][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[40][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[40][2]=0.098311600856684;
                        m_solPntsLocalCoord3D[41][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[41][1]=0.856669692468503;
                        m_solPntsLocalCoord3D[41][2]=0.098311600856684;
                        m_solPntsLocalCoord3D[42][0]=0.098311600856684;
                        m_solPntsLocalCoord3D[42][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[42][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[43][0]=0.856669692468503;
                        m_solPntsLocalCoord3D[43][1]=0.098311600856684;
                        m_solPntsLocalCoord3D[43][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[44][0]=0.660638331550837;
                        m_solPntsLocalCoord3D[44][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[44][2]=0.032124709578092;
                        m_solPntsLocalCoord3D[45][0]=0.660638331550837;
                        m_solPntsLocalCoord3D[45][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[45][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[46][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[46][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[46][2]=0.660638331550837;
                        m_solPntsLocalCoord3D[47][0]=0.032124709578092;
                        m_solPntsLocalCoord3D[47][1]=0.660638331550837;
                        m_solPntsLocalCoord3D[47][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[48][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[48][1]=0.032124709578092;
                        m_solPntsLocalCoord3D[48][2]=0.660638331550837;
                        m_solPntsLocalCoord3D[49][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[49][1]=0.660638331550837;
                        m_solPntsLocalCoord3D[49][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[50][0]=0.032124709578092;
                        m_solPntsLocalCoord3D[50][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[50][2]=0.660638331550837;
                        m_solPntsLocalCoord3D[51][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[51][1]=0.032124709578092;
                        m_solPntsLocalCoord3D[51][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[52][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[52][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[52][2]=0.032124709578092;
                        m_solPntsLocalCoord3D[53][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[53][1]=0.660638331550837;
                        m_solPntsLocalCoord3D[53][2]=0.032124709578092;
                        m_solPntsLocalCoord3D[54][0]=0.032124709578092;
                        m_solPntsLocalCoord3D[54][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[54][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[55][0]=0.660638331550837;
                        m_solPntsLocalCoord3D[55][1]=0.032124709578092;
                        m_solPntsLocalCoord3D[55][2]=0.153618479435536;
        /*m_solPntsLocalCoord3D[0][0]=0.0149520651530595;
        m_solPntsLocalCoord3D[0][1]=0.0149520651530595;
        m_solPntsLocalCoord3D[0][2]=0.955143804540822;
        m_solPntsLocalCoord3D[1][0]=0.955143804540822;
        m_solPntsLocalCoord3D[1][1]=0.0149520651530595;
        m_solPntsLocalCoord3D[1][2]=0.0149520651530595;
        m_solPntsLocalCoord3D[2][0]=0.0149520651530595;
        m_solPntsLocalCoord3D[2][1]=0.955143804540822;
        m_solPntsLocalCoord3D[2][2]=0.0149520651530595;
        m_solPntsLocalCoord3D[3][0]=0.0149520651530595;
        m_solPntsLocalCoord3D[3][1]=0.0149520651530595;
        m_solPntsLocalCoord3D[3][2]=0.0149520651530595;
        m_solPntsLocalCoord3D[4][0]=0.134478334792994;
        m_solPntsLocalCoord3D[4][1]=0.134478334792994;
        m_solPntsLocalCoord3D[4][2]=0.596564995621018;
        m_solPntsLocalCoord3D[5][0]=0.596564995621018;
        m_solPntsLocalCoord3D[5][1]=0.134478334792994;
        m_solPntsLocalCoord3D[5][2]=0.134478334792994;
        m_solPntsLocalCoord3D[6][0]=0.134478334792994;
        m_solPntsLocalCoord3D[6][1]=0.596564995621018;
        m_solPntsLocalCoord3D[6][2]=0.134478334792994;
        m_solPntsLocalCoord3D[7][0]=0.134478334792994;
        m_solPntsLocalCoord3D[7][1]=0.134478334792994;
        m_solPntsLocalCoord3D[7][2]=0.134478334792994;
        m_solPntsLocalCoord3D[8][0]=0.77997600844154;
        m_solPntsLocalCoord3D[8][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[8][2]=0.151831949165937;
        m_solPntsLocalCoord3D[9][0]=0.77997600844154;
        m_solPntsLocalCoord3D[9][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[9][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[10][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[10][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[10][2]=0.77997600844154;
        m_solPntsLocalCoord3D[11][0]=0.151831949165937;
        m_solPntsLocalCoord3D[11][1]=0.77997600844154;
        m_solPntsLocalCoord3D[11][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[12][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[12][1]=0.151831949165937;
        m_solPntsLocalCoord3D[12][2]=0.77997600844154;
        m_solPntsLocalCoord3D[13][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[13][1]=0.77997600844154;
        m_solPntsLocalCoord3D[13][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[14][0]=0.151831949165937;
        m_solPntsLocalCoord3D[14][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[14][2]=0.77997600844154;
        m_solPntsLocalCoord3D[15][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[15][1]=0.151831949165937;
        m_solPntsLocalCoord3D[15][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[16][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[16][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[16][2]=0.151831949165937;
        m_solPntsLocalCoord3D[17][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[17][1]=0.77997600844154;
        m_solPntsLocalCoord3D[17][2]=0.151831949165937;
        m_solPntsLocalCoord3D[18][0]=0.151831949165937;
        m_solPntsLocalCoord3D[18][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[18][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[19][0]=0.77997600844154;
        m_solPntsLocalCoord3D[19][1]=0.151831949165937;
        m_solPntsLocalCoord3D[19][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[20][0]=0.552655643106017;
        m_solPntsLocalCoord3D[20][1]=0.046205150415002;
        m_solPntsLocalCoord3D[20][2]=0.354934056063979;
        m_solPntsLocalCoord3D[21][0]=0.552655643106017;
        m_solPntsLocalCoord3D[21][1]=0.046205150415002;
        m_solPntsLocalCoord3D[21][2]=0.046205150415002;
        m_solPntsLocalCoord3D[22][0]=0.046205150415002;
        m_solPntsLocalCoord3D[22][1]=0.046205150415002;
        m_solPntsLocalCoord3D[22][2]=0.552655643106017;
        m_solPntsLocalCoord3D[23][0]=0.354934056063979;
        m_solPntsLocalCoord3D[23][1]=0.552655643106017;
        m_solPntsLocalCoord3D[23][2]=0.046205150415002;
        m_solPntsLocalCoord3D[24][0]=0.046205150415002;
        m_solPntsLocalCoord3D[24][1]=0.354934056063979;
        m_solPntsLocalCoord3D[24][2]=0.552655643106017;
        m_solPntsLocalCoord3D[25][0]=0.046205150415002;
        m_solPntsLocalCoord3D[25][1]=0.552655643106017;
        m_solPntsLocalCoord3D[25][2]=0.046205150415002;
        m_solPntsLocalCoord3D[26][0]=0.354934056063979;
        m_solPntsLocalCoord3D[26][1]=0.046205150415002;
        m_solPntsLocalCoord3D[26][2]=0.552655643106017;
        m_solPntsLocalCoord3D[27][0]=0.046205150415002;
        m_solPntsLocalCoord3D[27][1]=0.354934056063979;
        m_solPntsLocalCoord3D[27][2]=0.046205150415002;
        m_solPntsLocalCoord3D[28][0]=0.046205150415002;
        m_solPntsLocalCoord3D[28][1]=0.046205150415002;
        m_solPntsLocalCoord3D[28][2]=0.354934056063979;
        m_solPntsLocalCoord3D[29][0]=0.046205150415002;
        m_solPntsLocalCoord3D[29][1]=0.552655643106017;
        m_solPntsLocalCoord3D[29][2]=0.354934056063979;
        m_solPntsLocalCoord3D[30][0]=0.354934056063979;
        m_solPntsLocalCoord3D[30][1]=0.046205150415002;
        m_solPntsLocalCoord3D[30][2]=0.046205150415002;
        m_solPntsLocalCoord3D[31][0]=0.552655643106017;
        m_solPntsLocalCoord3D[31][1]=0.354934056063979;
        m_solPntsLocalCoord3D[31][2]=0.046205150415002;
        m_solPntsLocalCoord3D[32][0]=0.538104322888;
        m_solPntsLocalCoord3D[32][1]=0.228190461068761;
        m_solPntsLocalCoord3D[32][2]=0.00551475497447751;
        m_solPntsLocalCoord3D[33][0]=0.538104322888;
        m_solPntsLocalCoord3D[33][1]=0.228190461068761;
        m_solPntsLocalCoord3D[33][2]=0.228190461068761;
        m_solPntsLocalCoord3D[34][0]=0.228190461068761;
        m_solPntsLocalCoord3D[34][1]=0.228190461068761;
        m_solPntsLocalCoord3D[34][2]=0.538104322888;
        m_solPntsLocalCoord3D[35][0]=0.00551475497447751;
        m_solPntsLocalCoord3D[35][1]=0.538104322888;
        m_solPntsLocalCoord3D[35][2]=0.228190461068761;
        m_solPntsLocalCoord3D[36][0]=0.228190461068761;
        m_solPntsLocalCoord3D[36][1]=0.00551475497447751;
        m_solPntsLocalCoord3D[36][2]=0.538104322888;
        m_solPntsLocalCoord3D[37][0]=0.228190461068761;
        m_solPntsLocalCoord3D[37][1]=0.538104322888;
        m_solPntsLocalCoord3D[37][2]=0.228190461068761;
        m_solPntsLocalCoord3D[38][0]=0.00551475497447751;
        m_solPntsLocalCoord3D[38][1]=0.228190461068761;
        m_solPntsLocalCoord3D[38][2]=0.538104322888;
        m_solPntsLocalCoord3D[39][0]=0.228190461068761;
        m_solPntsLocalCoord3D[39][1]=0.00551475497447751;
        m_solPntsLocalCoord3D[39][2]=0.228190461068761;
        m_solPntsLocalCoord3D[40][0]=0.228190461068761;
        m_solPntsLocalCoord3D[40][1]=0.228190461068761;
        m_solPntsLocalCoord3D[40][2]=0.00551475497447751;
        m_solPntsLocalCoord3D[41][0]=0.228190461068761;
        m_solPntsLocalCoord3D[41][1]=0.538104322888;
        m_solPntsLocalCoord3D[41][2]=0.00551475497447751;
        m_solPntsLocalCoord3D[42][0]=0.00551475497447751;
        m_solPntsLocalCoord3D[42][1]=0.228190461068761;
        m_solPntsLocalCoord3D[42][2]=0.228190461068761;
        m_solPntsLocalCoord3D[43][0]=0.538104322888;
        m_solPntsLocalCoord3D[43][1]=0.00551475497447751;
        m_solPntsLocalCoord3D[43][2]=0.228190461068761;
        m_solPntsLocalCoord3D[44][0]=0.196183759574559;
        m_solPntsLocalCoord3D[44][1]=0.352305260087994;
        m_solPntsLocalCoord3D[44][2]=0.099205720249453;
        m_solPntsLocalCoord3D[45][0]=0.196183759574559;
        m_solPntsLocalCoord3D[45][1]=0.352305260087994;
        m_solPntsLocalCoord3D[45][2]=0.352305260087994;
        m_solPntsLocalCoord3D[46][0]=0.352305260087994;
        m_solPntsLocalCoord3D[46][1]=0.352305260087994;
        m_solPntsLocalCoord3D[46][2]=0.196183759574559;
        m_solPntsLocalCoord3D[47][0]=0.099205720249453;
        m_solPntsLocalCoord3D[47][1]=0.196183759574559;
        m_solPntsLocalCoord3D[47][2]=0.352305260087994;
        m_solPntsLocalCoord3D[48][0]=0.352305260087994;
        m_solPntsLocalCoord3D[48][1]=0.099205720249453;
        m_solPntsLocalCoord3D[48][2]=0.196183759574559;
        m_solPntsLocalCoord3D[49][0]=0.352305260087994;
        m_solPntsLocalCoord3D[49][1]=0.196183759574559;
        m_solPntsLocalCoord3D[49][2]=0.352305260087994;
        m_solPntsLocalCoord3D[50][0]=0.099205720249453;
        m_solPntsLocalCoord3D[50][1]=0.352305260087994;
        m_solPntsLocalCoord3D[50][2]=0.196183759574559;
        m_solPntsLocalCoord3D[51][0]=0.352305260087994;
        m_solPntsLocalCoord3D[51][1]=0.099205720249453;
        m_solPntsLocalCoord3D[51][2]=0.352305260087994;
        m_solPntsLocalCoord3D[52][0]=0.352305260087994;
        m_solPntsLocalCoord3D[52][1]=0.352305260087994;
        m_solPntsLocalCoord3D[52][2]=0.099205720249453;
        m_solPntsLocalCoord3D[53][0]=0.352305260087994;
        m_solPntsLocalCoord3D[53][1]=0.196183759574559;
        m_solPntsLocalCoord3D[53][2]=0.099205720249453;
        m_solPntsLocalCoord3D[54][0]=0.099205720249453;
        m_solPntsLocalCoord3D[54][1]=0.352305260087994;
        m_solPntsLocalCoord3D[54][2]=0.352305260087994;
        m_solPntsLocalCoord3D[55][0]=0.196183759574559;
        m_solPntsLocalCoord3D[55][1]=0.099205720249453;
        m_solPntsLocalCoord3D[55][2]=0.352305260087994;*/
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

TetraFluxReconstructionElementData::TetraFluxReconstructionElementData(CFPolyOrder::Type polyOrder, 
								     Common::SafePtr< BasePointDistribution > solPntDist, 
								     Common::SafePtr< BasePointDistribution > flxPntDist)
{
  m_shape = CFGeoShape::TETRA;
  m_dimensionality = DIM_3D;
  m_polyOrder = polyOrder;
  
  m_solPntsLocalCoord1D.resize(polyOrder+1);
  m_flxPntsLocalCoord1D.resize(polyOrder+1);

  // Use a default solution and flux point distribution: Gauss Legendre.  1D
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

//////////// 2D
//m_flxPntsLocalCoord2D.resize(((polyOrder+1)*(polyOrder+2))/2);

CFuint rows = (0.5*(polyOrder+1)*(polyOrder+2)) ;
m_flxPntsLocalCoord2D.resize(rows, std::vector< CFreal > (2));
switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
                        m_flxPntsLocalCoord2D[0][0] = 0.333333333333333; 
                        m_flxPntsLocalCoord2D[0][1] = 0.333333333333333;  
      } break;
      case CFPolyOrder::ORDER1:
      {
                        m_flxPntsLocalCoord2D[0][0] = 0.1666666666667; 
                        m_flxPntsLocalCoord2D[0][1] = 0.1666666666667;
                        m_flxPntsLocalCoord2D[1][0] = 0.6666666666667; 
                        m_flxPntsLocalCoord2D[1][1] = 0.1666666666667;  
                        m_flxPntsLocalCoord2D[2][0] = 0.1666666666667; 
                        m_flxPntsLocalCoord2D[2][1] = 0.6666666666667;

    } break;

    case CFPolyOrder::ORDER2:              
    {
                      m_flxPntsLocalCoord2D[0][0] =  0.091576213509780; 
                      m_flxPntsLocalCoord2D[0][1] =  0.091576213509780; 
                      
                      m_flxPntsLocalCoord2D[1][0] =  0.445948490915964; 
                      m_flxPntsLocalCoord2D[1][1] =  0.108103018168071; 

                      m_flxPntsLocalCoord2D[2][0] =  0.816847572980440; 
                      m_flxPntsLocalCoord2D[2][1] =  0.091576213509780;  

                      m_flxPntsLocalCoord2D[3][0] =  0.108103018168071; 
                      m_flxPntsLocalCoord2D[3][1] =  0.445948490915964;
                      
                      m_flxPntsLocalCoord2D[4][0] =  0.445948490915964; 
                      m_flxPntsLocalCoord2D[4][1] =  0.445948490915964;
                      
                      m_flxPntsLocalCoord2D[5][0] =  0.091576213509780; 
                      m_flxPntsLocalCoord2D[5][1] =  0.816847572980440;
                      
    } break;

    case CFPolyOrder::ORDER3:          
    {
                      m_flxPntsLocalCoord2D[0][0] = 0.055564052669793;
                      m_flxPntsLocalCoord2D[0][1] = 0.055564052669793;

                      m_flxPntsLocalCoord2D[1][0] = 0.295533711735893; 
                      m_flxPntsLocalCoord2D[1][1] = 0.070255540518384;  
                      
                      m_flxPntsLocalCoord2D[2][0] = 0.634210747745723; 
                      m_flxPntsLocalCoord2D[2][1] = 0.070255540518384;  	

                      m_flxPntsLocalCoord2D[3][0] = 0.888871894660413; 
                      m_flxPntsLocalCoord2D[3][1] = 0.055564052669793;	

                      m_flxPntsLocalCoord2D[4][0] = 0.070255540518384; 
                      m_flxPntsLocalCoord2D[4][1] = 0.295533711735893; 	

                      m_flxPntsLocalCoord2D[5][0] = 0.333333333333333; 
                      m_flxPntsLocalCoord2D[5][1] = 0.333333333333333;

                      m_flxPntsLocalCoord2D[6][0] = 0.634210747745723; 
                      m_flxPntsLocalCoord2D[6][1] = 0.295533711735893; 

                      m_flxPntsLocalCoord2D[7][0] = 0.070255540518384; 
                      m_flxPntsLocalCoord2D[7][1] = 0.634210747745723;	
                      
                      m_flxPntsLocalCoord2D[8][0] = 0.295533711735893; 
                      m_flxPntsLocalCoord2D[8][1] = 0.634210747745723; 
                      
                      m_flxPntsLocalCoord2D[9][0] = 0.055564052669793; 
                      m_flxPntsLocalCoord2D[9][1] = 0.888871894660413; 

      } break;
      case CFPolyOrder::ORDER4:   
      {
                        m_flxPntsLocalCoord2D[0][0] = 0.035870877695734; 
                        m_flxPntsLocalCoord2D[0][1] = 0.035870877695734;

                        m_flxPntsLocalCoord2D[1][0] = 0.201503881881800; 
                        m_flxPntsLocalCoord2D[1][1] =  0.047312487011716;

                        m_flxPntsLocalCoord2D[2][0] = 0.474308787777079;
                        m_flxPntsLocalCoord2D[2][1] = 0.051382424445843; 

                        m_flxPntsLocalCoord2D[3][0] = 0.751183631106484; 
                        m_flxPntsLocalCoord2D[3][1] = 0.047312487011716;

                        m_flxPntsLocalCoord2D[4][0] = 0.928258244608533; 
                        m_flxPntsLocalCoord2D[4][1] = 0.035870877695734;


                        m_flxPntsLocalCoord2D[5][0] = 0.047312487011716;
                        m_flxPntsLocalCoord2D[5][1] = 0.201503881881800;

                        m_flxPntsLocalCoord2D[6][0] = 0.241729395767967;
                        m_flxPntsLocalCoord2D[6][1] = 0.241729395767967;

                        m_flxPntsLocalCoord2D[7][0] = 0.516541208464066; 
                        m_flxPntsLocalCoord2D[7][1] = 0.241729395767967; 

                        m_flxPntsLocalCoord2D[8][0] = 0.751183631106484; 
                        m_flxPntsLocalCoord2D[8][1] = 0.201503881881800;

                        m_flxPntsLocalCoord2D[9][0] = 0.051382424445843;
                        m_flxPntsLocalCoord2D[9][1] = 0.474308787777079;

                        m_flxPntsLocalCoord2D[10][0] = 0.241729395767967; 
                        m_flxPntsLocalCoord2D[10][1] = 0.516541208464066; 

                        m_flxPntsLocalCoord2D[11][0] = 0.474308787777079; 
                        m_flxPntsLocalCoord2D[11][1] = 0.474308787777079;

                        m_flxPntsLocalCoord2D[12][0] = 0.047312487011716; 
                        m_flxPntsLocalCoord2D[12][1] = 0.751183631106484;

                        m_flxPntsLocalCoord2D[13][0] =  0.201503881881800;
                        m_flxPntsLocalCoord2D[13][1] =  0.751183631106484;

                        m_flxPntsLocalCoord2D[14][0] = 0.035870877695734; 
                        m_flxPntsLocalCoord2D[14][1] = 0.928258244608533; 

      } break;

      case CFPolyOrder::ORDER5:
      {
                        m_flxPntsLocalCoord2D[0][0] = 0.028112952182664; 
                        m_flxPntsLocalCoord2D[0][1] = 0.028112952182664; 
 
                        m_flxPntsLocalCoord2D[1][0] = 0.148565812270887; 
                        m_flxPntsLocalCoord2D[1][1] = 0.033533207700614; 

                        m_flxPntsLocalCoord2D[2][0] = 0.357196298615681; 
                        m_flxPntsLocalCoord2D[2][1] = 0.037824789609186;

                        m_flxPntsLocalCoord2D[3][0] = 0.604978911775132; 
                        m_flxPntsLocalCoord2D[3][1] = 0.037824789609186;

                        m_flxPntsLocalCoord2D[4][0] = 0.817900980028499; 
                        m_flxPntsLocalCoord2D[4][1] = 0.033533207700614;

                        m_flxPntsLocalCoord2D[5][0] = 0.943774095634672;
                        m_flxPntsLocalCoord2D[5][1] = 0.028112952182664;

                        m_flxPntsLocalCoord2D[6][0] = 0.033533207700614;
                        m_flxPntsLocalCoord2D[6][1] = 0.148565812270887;

                        m_flxPntsLocalCoord2D[7][0] = 0.177139098469317;
                        m_flxPntsLocalCoord2D[7][1] = 0.177139098469317;

                        m_flxPntsLocalCoord2D[8][0] = 0.405508595867433;
                        m_flxPntsLocalCoord2D[8][1] = 0.188982808265134;

                        m_flxPntsLocalCoord2D[9][0] = 0.645721803061365; 
                        m_flxPntsLocalCoord2D[9][1] = 0.177139098469317;

                        m_flxPntsLocalCoord2D[10][0] = 0.817900980028499;
                        m_flxPntsLocalCoord2D[10][1] = 0.148565812270887;


                        m_flxPntsLocalCoord2D[11][0] = 0.037824789609186; 
                        m_flxPntsLocalCoord2D[11][1] = 0.357196298615681; 

                        m_flxPntsLocalCoord2D[12][0] = 0.188982808265134;
                        m_flxPntsLocalCoord2D[12][1] = 0.405508595867433;

                        m_flxPntsLocalCoord2D[13][0] = 0.405508595867433;
                        m_flxPntsLocalCoord2D[13][1] = 0.405508595867433;

                        m_flxPntsLocalCoord2D[14][0] = 0.604978911775132; 
                        m_flxPntsLocalCoord2D[14][1] = 0.357196298615681; 


                        m_flxPntsLocalCoord2D[15][0] = 0.037824789609186; 
                        m_flxPntsLocalCoord2D[15][1] = 0.604978911775132;

                        m_flxPntsLocalCoord2D[16][0] = 0.177139098469317;
                        m_flxPntsLocalCoord2D[16][1] = 0.645721803061365;

                        m_flxPntsLocalCoord2D[17][0] = 0.357196298615681;
                        m_flxPntsLocalCoord2D[17][1] = 0.604978911775132;

                        m_flxPntsLocalCoord2D[18][0] = 0.033533207700614; 
                        m_flxPntsLocalCoord2D[18][1] = 0.817900980028499;

                        m_flxPntsLocalCoord2D[19][0] = 0.148565812270887;
                        m_flxPntsLocalCoord2D[19][1] = 0.817900980028499;

                        m_flxPntsLocalCoord2D[20][0] = 0.028112952182664;
                        m_flxPntsLocalCoord2D[20][1] = 0.943774095634672;

      } break;


   case CFPolyOrder::ORDER6:
      {  
                        m_flxPntsLocalCoord2D[0][0] =  0.0000000000000;
                        m_flxPntsLocalCoord2D[0][1] =  0.9451704450174; 
                        m_flxPntsLocalCoord2D[1][0] =  0.9451704450173; 
                        m_flxPntsLocalCoord2D[1][1] =  0.0000000000000; 
                        m_flxPntsLocalCoord2D[2][0] =  0.9289002405719;
                        m_flxPntsLocalCoord2D[2][1] =  0.0685505797224;
                        m_flxPntsLocalCoord2D[3][0] =  0.0685505797224;
                        m_flxPntsLocalCoord2D[3][1] =  0.9289002405717;
                        m_flxPntsLocalCoord2D[4][0] =  0.0243268355615;
                        m_flxPntsLocalCoord2D[4][1] =  0.0243268355616; 
                        m_flxPntsLocalCoord2D[5][0] =  0.1279662835335;
                        m_flxPntsLocalCoord2D[5][1] =  0.0277838749488;
                        m_flxPntsLocalCoord2D[6][0] =  0.0277838749488; 
                        m_flxPntsLocalCoord2D[6][1] =  0.1279662835337; 
                        m_flxPntsLocalCoord2D[7][0] =  0.0287083428360; 
                        m_flxPntsLocalCoord2D[7][1] =  0.7498347588657;
                        m_flxPntsLocalCoord2D[8][0] =  0.7498347588656;
                        m_flxPntsLocalCoord2D[8][1] =  0.0287083428360;
                        m_flxPntsLocalCoord2D[9][0] =  0.7228007909707;
                        m_flxPntsLocalCoord2D[9][1] =  0.2497602062385;
                        m_flxPntsLocalCoord2D[10][0] = 0.2497602062386; 
                        m_flxPntsLocalCoord2D[10][1] = 0.7228007909707;
                        m_flxPntsLocalCoord2D[11][0] = 0.0865562992839;
                        m_flxPntsLocalCoord2D[11][1] = 0.8325513856997;
                        m_flxPntsLocalCoord2D[12][0] = 0.8325513856998;
                        m_flxPntsLocalCoord2D[12][1] = 0.0865562992839;
                        m_flxPntsLocalCoord2D[13][0] = 0.3061619157672;
                        m_flxPntsLocalCoord2D[13][1] = 0.0303526617491;
                        m_flxPntsLocalCoord2D[14][0] = 0.0303526617491;
                        m_flxPntsLocalCoord2D[14][1] = 0.3061619157675;
                        m_flxPntsLocalCoord2D[15][0] = 0.4868610595047;
                        m_flxPntsLocalCoord2D[15][1] = 0.4868610595047;
                        m_flxPntsLocalCoord2D[16][0] = 0.6657904293017;
                        m_flxPntsLocalCoord2D[16][1] = 0.1765456154219;
                        m_flxPntsLocalCoord2D[17][0] = 0.1765456154221; 
                        m_flxPntsLocalCoord2D[17][1] = 0.6657904293016;
                        m_flxPntsLocalCoord2D[18][0] = 0.0293121007360; 
                        m_flxPntsLocalCoord2D[18][1] = 0.5295657488669;
                        m_flxPntsLocalCoord2D[19][0] = 0.5295657488667;
                        m_flxPntsLocalCoord2D[19][1] = 0.0293121007360;
                        m_flxPntsLocalCoord2D[20][0] = 0.1444673824391; 
                        m_flxPntsLocalCoord2D[20][1] = 0.1444673824391;
                        m_flxPntsLocalCoord2D[21][0] = 0.3299740111411; 
                        m_flxPntsLocalCoord2D[21][1] = 0.5361815729050;
                        m_flxPntsLocalCoord2D[22][0] = 0.5361815729052; 
                        m_flxPntsLocalCoord2D[22][1] = 0.3299740111409;
                        m_flxPntsLocalCoord2D[23][0] = 0.5511507516862;
                        m_flxPntsLocalCoord2D[23][1] = 0.1437790861923;
                        m_flxPntsLocalCoord2D[24][0] = 0.1437790861923;
                        m_flxPntsLocalCoord2D[24][1] = 0.5511507516862;
                        m_flxPntsLocalCoord2D[25][0] = 0.3348066587327; 
                        m_flxPntsLocalCoord2D[25][1] = 0.1529619437161;
                        m_flxPntsLocalCoord2D[26][0] = 0.1529619437161; 
                        m_flxPntsLocalCoord2D[26][1] = 0.3348066587327;
                        m_flxPntsLocalCoord2D[27][0] = 0.3430183498147; 
                        m_flxPntsLocalCoord2D[27][1] = 0.3430183498147;
                      
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Quadrature Points not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
    }

    ///////////////////////// 3D

    //m_solPntsLocalCoord3D.resize(((polyOrder+1)*(polyOrder+2)*(polyOrder+3))/6);
    
    CFuint col = (((polyOrder+1)*(polyOrder+2)*(polyOrder+3))/6) ;
    m_solPntsLocalCoord3D.resize(col, std::vector< CFreal > (3));

    switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
        m_solPntsLocalCoord3D[0][0] = 0.25; 
        m_solPntsLocalCoord3D[0][1] = 0.25;
        m_solPntsLocalCoord3D[0][2] = 0.25;  
      } break;
      case CFPolyOrder::ORDER1:
      {
        m_solPntsLocalCoord3D[0][0]=0.1382;
        m_solPntsLocalCoord3D[0][1]=0.1382;
        m_solPntsLocalCoord3D[0][2]=0.1382;
        m_solPntsLocalCoord3D[1][0]=0.58541;
        m_solPntsLocalCoord3D[1][1]=0.1382;
        m_solPntsLocalCoord3D[1][2]=0.1382;
        m_solPntsLocalCoord3D[2][0]=0.1382;
        m_solPntsLocalCoord3D[2][1]=0.58541;
        m_solPntsLocalCoord3D[2][2]=0.1382;
        m_solPntsLocalCoord3D[3][0]=0.1382;
        m_solPntsLocalCoord3D[3][1]=0.1382;
        m_solPntsLocalCoord3D[3][2]=0.58541;
      } break;
      case CFPolyOrder::ORDER2:              
      {
        m_solPntsLocalCoord3D[0][0]=0.073835;
        m_solPntsLocalCoord3D[0][1]=0.073835;
        m_solPntsLocalCoord3D[0][2]=0.073835;
        m_solPntsLocalCoord3D[1][0]=0.7785;
        m_solPntsLocalCoord3D[1][1]=0.073835;
        m_solPntsLocalCoord3D[1][2]=0.073835;
        m_solPntsLocalCoord3D[2][0]=0.073835;
        m_solPntsLocalCoord3D[2][1]=0.7785;
        m_solPntsLocalCoord3D[2][2]=0.073835;
        m_solPntsLocalCoord3D[3][0]=0.073835;
        m_solPntsLocalCoord3D[3][1]=0.073835;
        m_solPntsLocalCoord3D[3][2]=0.7785;
        m_solPntsLocalCoord3D[4][0]=0.40624;
        m_solPntsLocalCoord3D[4][1]=0.093756;
        m_solPntsLocalCoord3D[4][2]=0.093756;
        m_solPntsLocalCoord3D[5][0]=0.093756;
        m_solPntsLocalCoord3D[5][1]=0.40624;
        m_solPntsLocalCoord3D[5][2]=0.093756;
        m_solPntsLocalCoord3D[6][0]=0.093756;
        m_solPntsLocalCoord3D[6][1]=0.093756;
        m_solPntsLocalCoord3D[6][2]=0.40624;
        m_solPntsLocalCoord3D[7][0]=0.40624;
        m_solPntsLocalCoord3D[7][1]=0.40624;
        m_solPntsLocalCoord3D[7][2]=0.093756;
        m_solPntsLocalCoord3D[8][0]=0.40624;
        m_solPntsLocalCoord3D[8][1]=0.093756;
        m_solPntsLocalCoord3D[8][2]=0.40624;
        m_solPntsLocalCoord3D[9][0]=0.093756;
        m_solPntsLocalCoord3D[9][1]=0.40624;
        m_solPntsLocalCoord3D[9][2]=0.40624;
      } break;
      case CFPolyOrder::ORDER3:              
      {
        m_solPntsLocalCoord3D[0][0]=0.032353;
        m_solPntsLocalCoord3D[0][1]=0.032353;
        m_solPntsLocalCoord3D[0][2]=0.032353;
        m_solPntsLocalCoord3D[1][0]=0.90294;
        m_solPntsLocalCoord3D[1][1]=0.032353;
        m_solPntsLocalCoord3D[1][2]=0.032353;
        m_solPntsLocalCoord3D[2][0]=0.032353;
        m_solPntsLocalCoord3D[2][1]=0.90294;
        m_solPntsLocalCoord3D[2][2]=0.032353;
        m_solPntsLocalCoord3D[3][0]=0.032353;
        m_solPntsLocalCoord3D[3][1]=0.032353;
        m_solPntsLocalCoord3D[3][2]=0.90294;
        m_solPntsLocalCoord3D[4][0]=0.6166;
        m_solPntsLocalCoord3D[4][1]=0.06036;
        m_solPntsLocalCoord3D[4][2]=0.06036;
        m_solPntsLocalCoord3D[5][0]=0.26268;
        m_solPntsLocalCoord3D[5][1]=0.06036;
        m_solPntsLocalCoord3D[5][2]=0.06036;
        m_solPntsLocalCoord3D[6][0]=0.06036;
        m_solPntsLocalCoord3D[6][1]=0.6166;
        m_solPntsLocalCoord3D[6][2]=0.06036;
        m_solPntsLocalCoord3D[7][0]=0.06036;
        m_solPntsLocalCoord3D[7][1]=0.26268;
        m_solPntsLocalCoord3D[7][2]=0.06036;
        m_solPntsLocalCoord3D[8][0]=0.06036;
        m_solPntsLocalCoord3D[8][1]=0.06036;
        m_solPntsLocalCoord3D[8][2]=0.6166;
        m_solPntsLocalCoord3D[9][0]=0.06036;
        m_solPntsLocalCoord3D[9][1]=0.06036;
        m_solPntsLocalCoord3D[9][2]=0.26268;
        m_solPntsLocalCoord3D[10][0]=0.26268;
        m_solPntsLocalCoord3D[10][1]=0.6166;
        m_solPntsLocalCoord3D[10][2]=0.06036;
        m_solPntsLocalCoord3D[11][0]=0.6166;
        m_solPntsLocalCoord3D[11][1]=0.26268;
        m_solPntsLocalCoord3D[11][2]=0.06036;
        m_solPntsLocalCoord3D[12][0]=0.26268;
        m_solPntsLocalCoord3D[12][1]=0.06036;
        m_solPntsLocalCoord3D[12][2]=0.6166;
        m_solPntsLocalCoord3D[13][0]=0.6166;
        m_solPntsLocalCoord3D[13][1]=0.06036;
        m_solPntsLocalCoord3D[13][2]=0.26268;
        m_solPntsLocalCoord3D[14][0]=0.06036;
        m_solPntsLocalCoord3D[14][1]=0.26268;
        m_solPntsLocalCoord3D[14][2]=0.6166;
        m_solPntsLocalCoord3D[15][0]=0.06036;
        m_solPntsLocalCoord3D[15][1]=0.6166;
        m_solPntsLocalCoord3D[15][2]=0.26268;
        m_solPntsLocalCoord3D[16][0]=0.30977;
        m_solPntsLocalCoord3D[16][1]=0.30977;
        m_solPntsLocalCoord3D[16][2]=0.070692;
        m_solPntsLocalCoord3D[17][0]=0.30977;
        m_solPntsLocalCoord3D[17][1]=0.070692;
        m_solPntsLocalCoord3D[17][2]=0.30977;
        m_solPntsLocalCoord3D[18][0]=0.070692;
        m_solPntsLocalCoord3D[18][1]=0.30977;
        m_solPntsLocalCoord3D[18][2]=0.30977;
        m_solPntsLocalCoord3D[19][0]=0.30977;
        m_solPntsLocalCoord3D[19][1]=0.30977;
        m_solPntsLocalCoord3D[19][2]=0.30977;
      } break;
      case CFPolyOrder::ORDER4:              
      {
        m_solPntsLocalCoord3D[0][0]=0.25;
        m_solPntsLocalCoord3D[0][1]=0.25;
        m_solPntsLocalCoord3D[0][2]=0.25;
        m_solPntsLocalCoord3D[1][0]=0.0267367755543735;
        m_solPntsLocalCoord3D[1][1]=0.0267367755543735;
        m_solPntsLocalCoord3D[1][2]=0.91978967333688;
        m_solPntsLocalCoord3D[2][0]=0.91978967333688;
        m_solPntsLocalCoord3D[2][1]=0.0267367755543735;
        m_solPntsLocalCoord3D[2][2]=0.0267367755543735;
        m_solPntsLocalCoord3D[3][0]=0.0267367755543735;
        m_solPntsLocalCoord3D[3][1]=0.91978967333688;
        m_solPntsLocalCoord3D[3][2]=0.0267367755543735;
        m_solPntsLocalCoord3D[4][0]=0.0267367755543735;
        m_solPntsLocalCoord3D[4][1]=0.0267367755543735;
        m_solPntsLocalCoord3D[4][2]=0.0267367755543735;
        m_solPntsLocalCoord3D[5][0]=0.0452454000155175;
        m_solPntsLocalCoord3D[5][1]=0.454754599984483;
        m_solPntsLocalCoord3D[5][2]=0.0452454000155175;
        m_solPntsLocalCoord3D[6][0]=0.454754599984483;
        m_solPntsLocalCoord3D[6][1]=0.0452454000155175;
        m_solPntsLocalCoord3D[6][2]=0.0452454000155175;
        m_solPntsLocalCoord3D[7][0]=0.0452454000155175;
        m_solPntsLocalCoord3D[7][1]=0.454754599984483;
        m_solPntsLocalCoord3D[7][2]=0.454754599984483;
        m_solPntsLocalCoord3D[8][0]=0.454754599984483;
        m_solPntsLocalCoord3D[8][1]=0.0452454000155175;
        m_solPntsLocalCoord3D[8][2]=0.454754599984483;
        m_solPntsLocalCoord3D[9][0]=0.454754599984483;
        m_solPntsLocalCoord3D[9][1]=0.454754599984483;
        m_solPntsLocalCoord3D[9][2]=0.0452454000155175;
        m_solPntsLocalCoord3D[10][0]=0.0452454000155175;
        m_solPntsLocalCoord3D[10][1]=0.0452454000155175;
        m_solPntsLocalCoord3D[10][2]=0.454754599984483;
        m_solPntsLocalCoord3D[11][0]=0.747759888481808;
        m_solPntsLocalCoord3D[11][1]=0.039102240635649;
        m_solPntsLocalCoord3D[11][2]=0.174035630246894;
        m_solPntsLocalCoord3D[12][0]=0.747759888481808;
        m_solPntsLocalCoord3D[12][1]=0.039102240635649;
        m_solPntsLocalCoord3D[12][2]=0.039102240635649;
        m_solPntsLocalCoord3D[13][0]=0.039102240635649;
        m_solPntsLocalCoord3D[13][1]=0.039102240635649;
        m_solPntsLocalCoord3D[13][2]=0.747759888481808;
        m_solPntsLocalCoord3D[14][0]=0.174035630246894;
        m_solPntsLocalCoord3D[14][1]=0.747759888481808;
        m_solPntsLocalCoord3D[14][2]=0.039102240635649;
        m_solPntsLocalCoord3D[15][0]=0.039102240635649;
        m_solPntsLocalCoord3D[15][1]=0.174035630246894;
        m_solPntsLocalCoord3D[15][2]=0.747759888481808;
        m_solPntsLocalCoord3D[16][0]=0.039102240635649;
        m_solPntsLocalCoord3D[16][1]=0.747759888481808;
        m_solPntsLocalCoord3D[16][2]=0.039102240635649;
        m_solPntsLocalCoord3D[17][0]=0.174035630246894;
        m_solPntsLocalCoord3D[17][1]=0.039102240635649;
        m_solPntsLocalCoord3D[17][2]=0.747759888481808;
        m_solPntsLocalCoord3D[18][0]=0.039102240635649;
        m_solPntsLocalCoord3D[18][1]=0.174035630246894;
        m_solPntsLocalCoord3D[18][2]=0.039102240635649;
        m_solPntsLocalCoord3D[19][0]=0.039102240635649;
        m_solPntsLocalCoord3D[19][1]=0.039102240635649;
        m_solPntsLocalCoord3D[19][2]=0.174035630246894;
        m_solPntsLocalCoord3D[20][0]=0.039102240635649;
        m_solPntsLocalCoord3D[20][1]=0.747759888481808;
        m_solPntsLocalCoord3D[20][2]=0.174035630246894;
        m_solPntsLocalCoord3D[21][0]=0.174035630246894;
        m_solPntsLocalCoord3D[21][1]=0.039102240635649;
        m_solPntsLocalCoord3D[21][2]=0.039102240635649;
        m_solPntsLocalCoord3D[22][0]=0.747759888481808;
        m_solPntsLocalCoord3D[22][1]=0.174035630246894;
        m_solPntsLocalCoord3D[22][2]=0.039102240635649;
        m_solPntsLocalCoord3D[23][0]=0.503118645014598;
        m_solPntsLocalCoord3D[23][1]=0.223201037962315;
        m_solPntsLocalCoord3D[23][2]=0.050479279060772;
        m_solPntsLocalCoord3D[24][0]=0.503118645014598;
        m_solPntsLocalCoord3D[24][1]=0.223201037962315;
        m_solPntsLocalCoord3D[24][2]=0.223201037962315;
        m_solPntsLocalCoord3D[25][0]=0.223201037962315;
        m_solPntsLocalCoord3D[25][1]=0.223201037962315;
        m_solPntsLocalCoord3D[25][2]=0.503118645014598;
        m_solPntsLocalCoord3D[26][0]=0.050479279060772;
        m_solPntsLocalCoord3D[26][1]=0.503118645014598;
        m_solPntsLocalCoord3D[26][2]=0.223201037962315;
        m_solPntsLocalCoord3D[27][0]=0.223201037962315;
        m_solPntsLocalCoord3D[27][1]=0.050479279060772;
        m_solPntsLocalCoord3D[27][2]=0.503118645014598;
        m_solPntsLocalCoord3D[28][0]=0.223201037962315;
        m_solPntsLocalCoord3D[28][1]=0.503118645014598;
        m_solPntsLocalCoord3D[28][2]=0.223201037962315;
        m_solPntsLocalCoord3D[29][0]=0.050479279060772;
        m_solPntsLocalCoord3D[29][1]=0.223201037962315;
        m_solPntsLocalCoord3D[29][2]=0.503118645014598;
        m_solPntsLocalCoord3D[30][0]=0.223201037962315;
        m_solPntsLocalCoord3D[30][1]=0.050479279060772;
        m_solPntsLocalCoord3D[30][2]=0.223201037962315;
        m_solPntsLocalCoord3D[31][0]=0.223201037962315;
        m_solPntsLocalCoord3D[31][1]=0.223201037962315;
        m_solPntsLocalCoord3D[31][2]=0.050479279060772;
        m_solPntsLocalCoord3D[32][0]=0.223201037962315;
        m_solPntsLocalCoord3D[32][1]=0.503118645014598;
        m_solPntsLocalCoord3D[32][2]=0.050479279060772;
        m_solPntsLocalCoord3D[33][0]=0.050479279060772;
        m_solPntsLocalCoord3D[33][1]=0.223201037962315;
        m_solPntsLocalCoord3D[33][2]=0.223201037962315;
        m_solPntsLocalCoord3D[34][0]=0.503118645014598;
        m_solPntsLocalCoord3D[34][1]=0.050479279060772;
        m_solPntsLocalCoord3D[34][2]=0.223201037962315;
      } break;
      case CFPolyOrder::ORDER5:              
      {
                        m_solPntsLocalCoord3D[0][0]=0.175924928891649;
                        m_solPntsLocalCoord3D[0][1]=0.175924928891649;
                        m_solPntsLocalCoord3D[0][2]=0.472225213325055;
                        m_solPntsLocalCoord3D[1][0]=0.175924928891649;
                        m_solPntsLocalCoord3D[1][1]=0.472225213325055;
                        m_solPntsLocalCoord3D[1][2]=0.175924928891649;
                        m_solPntsLocalCoord3D[2][0]=0.472225213325055;
                        m_solPntsLocalCoord3D[2][1]=0.175924928891649;
                        m_solPntsLocalCoord3D[2][2]=0.175924928891649;
                        m_solPntsLocalCoord3D[3][0]=0.175924928891649;
                        m_solPntsLocalCoord3D[3][1]=0.175924928891649;
                        m_solPntsLocalCoord3D[3][2]=0.175924928891649;
                        m_solPntsLocalCoord3D[4][0]=0.214336612970414;
                        m_solPntsLocalCoord3D[4][1]=0.214336612970414;
                        m_solPntsLocalCoord3D[4][2]=0.356990161088759;
                        m_solPntsLocalCoord3D[5][0]=0.214336612970414;
                        m_solPntsLocalCoord3D[5][1]=0.356990161088759;
                        m_solPntsLocalCoord3D[5][2]=0.214336612970414;
                        m_solPntsLocalCoord3D[6][0]=0.356990161088759;
                        m_solPntsLocalCoord3D[6][1]=0.214336612970414;
                        m_solPntsLocalCoord3D[6][2]=0.214336612970414;
                        m_solPntsLocalCoord3D[7][0]=0.214336612970414;
                        m_solPntsLocalCoord3D[7][1]=0.214336612970414;
                        m_solPntsLocalCoord3D[7][2]=0.214336612970414;
                        m_solPntsLocalCoord3D[8][0]=0.174463680162906;
                        m_solPntsLocalCoord3D[8][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[8][2]=0.0400417579050455;
                        m_solPntsLocalCoord3D[9][0]=0.174463680162906;
                        m_solPntsLocalCoord3D[9][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[9][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[10][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[10][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[10][2]=0.174463680162906;
                        m_solPntsLocalCoord3D[11][0]=0.0400417579050455;
                        m_solPntsLocalCoord3D[11][1]=0.174463680162906;
                        m_solPntsLocalCoord3D[11][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[12][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[12][1]=0.0400417579050455;
                        m_solPntsLocalCoord3D[12][2]=0.174463680162906;
                        m_solPntsLocalCoord3D[13][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[13][1]=0.174463680162906;
                        m_solPntsLocalCoord3D[13][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[14][0]=0.0400417579050455;
                        m_solPntsLocalCoord3D[14][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[14][2]=0.174463680162906;
                        m_solPntsLocalCoord3D[15][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[15][1]=0.0400417579050455;
                        m_solPntsLocalCoord3D[15][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[16][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[16][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[16][2]=0.0400417579050455;
                        m_solPntsLocalCoord3D[17][0]=0.392747280966025;
                        m_solPntsLocalCoord3D[17][1]=0.174463680162906;
                        m_solPntsLocalCoord3D[17][2]=0.0400417579050455;
                        m_solPntsLocalCoord3D[18][0]=0.0400417579050455;
                        m_solPntsLocalCoord3D[18][1]=0.392747280966025;
                        m_solPntsLocalCoord3D[18][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[19][0]=0.174463680162906;
                        m_solPntsLocalCoord3D[19][1]=0.0400417579050455;
                        m_solPntsLocalCoord3D[19][2]=0.392747280966025;
                        m_solPntsLocalCoord3D[20][0]=0.611348339064507;
                        m_solPntsLocalCoord3D[20][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[20][2]=0.325129236353862;
                        m_solPntsLocalCoord3D[21][0]=0.611348339064507;
                        m_solPntsLocalCoord3D[21][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[21][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[22][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[22][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[22][2]=0.611348339064507;
                        m_solPntsLocalCoord3D[23][0]=0.325129236353862;
                        m_solPntsLocalCoord3D[23][1]=0.611348339064507;
                        m_solPntsLocalCoord3D[23][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[24][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[24][1]=0.325129236353862;
                        m_solPntsLocalCoord3D[24][2]=0.611348339064507;
                        m_solPntsLocalCoord3D[25][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[25][1]=0.611348339064507;
                        m_solPntsLocalCoord3D[25][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[26][0]=0.325129236353862;
                        m_solPntsLocalCoord3D[26][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[26][2]=0.611348339064507;
                        m_solPntsLocalCoord3D[27][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[27][1]=0.325129236353862;
                        m_solPntsLocalCoord3D[27][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[28][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[28][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[28][2]=0.325129236353862;
                        m_solPntsLocalCoord3D[29][0]=0.0317612122908155;
                        m_solPntsLocalCoord3D[29][1]=0.611348339064507;
                        m_solPntsLocalCoord3D[29][2]=0.325129236353862;
                        m_solPntsLocalCoord3D[30][0]=0.325129236353862;
                        m_solPntsLocalCoord3D[30][1]=0.0317612122908155;
                        m_solPntsLocalCoord3D[30][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[31][0]=0.611348339064507;
                        m_solPntsLocalCoord3D[31][1]=0.325129236353862;
                        m_solPntsLocalCoord3D[31][2]=0.0317612122908155;
                        m_solPntsLocalCoord3D[32][0]=0.856669692468503;
                        m_solPntsLocalCoord3D[32][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[32][2]=0.098311600856684;
                        m_solPntsLocalCoord3D[33][0]=0.856669692468503;
                        m_solPntsLocalCoord3D[33][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[33][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[34][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[34][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[34][2]=0.856669692468503;
                        m_solPntsLocalCoord3D[35][0]=0.098311600856684;
                        m_solPntsLocalCoord3D[35][1]=0.856669692468503;
                        m_solPntsLocalCoord3D[35][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[36][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[36][1]=0.098311600856684;
                        m_solPntsLocalCoord3D[36][2]=0.856669692468503;
                        m_solPntsLocalCoord3D[37][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[37][1]=0.856669692468503;
                        m_solPntsLocalCoord3D[37][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[38][0]=0.098311600856684;
                        m_solPntsLocalCoord3D[38][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[38][2]=0.856669692468503;
                        m_solPntsLocalCoord3D[39][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[39][1]=0.098311600856684;
                        m_solPntsLocalCoord3D[39][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[40][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[40][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[40][2]=0.098311600856684;
                        m_solPntsLocalCoord3D[41][0]=0.0225093533374065;
                        m_solPntsLocalCoord3D[41][1]=0.856669692468503;
                        m_solPntsLocalCoord3D[41][2]=0.098311600856684;
                        m_solPntsLocalCoord3D[42][0]=0.098311600856684;
                        m_solPntsLocalCoord3D[42][1]=0.0225093533374065;
                        m_solPntsLocalCoord3D[42][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[43][0]=0.856669692468503;
                        m_solPntsLocalCoord3D[43][1]=0.098311600856684;
                        m_solPntsLocalCoord3D[43][2]=0.0225093533374065;
                        m_solPntsLocalCoord3D[44][0]=0.660638331550837;
                        m_solPntsLocalCoord3D[44][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[44][2]=0.032124709578092;
                        m_solPntsLocalCoord3D[45][0]=0.660638331550837;
                        m_solPntsLocalCoord3D[45][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[45][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[46][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[46][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[46][2]=0.660638331550837;
                        m_solPntsLocalCoord3D[47][0]=0.032124709578092;
                        m_solPntsLocalCoord3D[47][1]=0.660638331550837;
                        m_solPntsLocalCoord3D[47][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[48][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[48][1]=0.032124709578092;
                        m_solPntsLocalCoord3D[48][2]=0.660638331550837;
                        m_solPntsLocalCoord3D[49][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[49][1]=0.660638331550837;
                        m_solPntsLocalCoord3D[49][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[50][0]=0.032124709578092;
                        m_solPntsLocalCoord3D[50][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[50][2]=0.660638331550837;
                        m_solPntsLocalCoord3D[51][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[51][1]=0.032124709578092;
                        m_solPntsLocalCoord3D[51][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[52][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[52][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[52][2]=0.032124709578092;
                        m_solPntsLocalCoord3D[53][0]=0.153618479435536;
                        m_solPntsLocalCoord3D[53][1]=0.660638331550837;
                        m_solPntsLocalCoord3D[53][2]=0.032124709578092;
                        m_solPntsLocalCoord3D[54][0]=0.032124709578092;
                        m_solPntsLocalCoord3D[54][1]=0.153618479435536;
                        m_solPntsLocalCoord3D[54][2]=0.153618479435536;
                        m_solPntsLocalCoord3D[55][0]=0.660638331550837;
                        m_solPntsLocalCoord3D[55][1]=0.032124709578092;
                        m_solPntsLocalCoord3D[55][2]=0.153618479435536;
      /*  m_solPntsLocalCoord3D[0][0]=0.0149520651530595;
        m_solPntsLocalCoord3D[0][1]=0.0149520651530595;
        m_solPntsLocalCoord3D[0][2]=0.955143804540822;
        m_solPntsLocalCoord3D[1][0]=0.955143804540822;
        m_solPntsLocalCoord3D[1][1]=0.0149520651530595;
        m_solPntsLocalCoord3D[1][2]=0.0149520651530595;
        m_solPntsLocalCoord3D[2][0]=0.0149520651530595;
        m_solPntsLocalCoord3D[2][1]=0.955143804540822;
        m_solPntsLocalCoord3D[2][2]=0.0149520651530595;
        m_solPntsLocalCoord3D[3][0]=0.0149520651530595;
        m_solPntsLocalCoord3D[3][1]=0.0149520651530595;
        m_solPntsLocalCoord3D[3][2]=0.0149520651530595;
        m_solPntsLocalCoord3D[4][0]=0.134478334792994;
        m_solPntsLocalCoord3D[4][1]=0.134478334792994;
        m_solPntsLocalCoord3D[4][2]=0.596564995621018;
        m_solPntsLocalCoord3D[5][0]=0.596564995621018;
        m_solPntsLocalCoord3D[5][1]=0.134478334792994;
        m_solPntsLocalCoord3D[5][2]=0.134478334792994;
        m_solPntsLocalCoord3D[6][0]=0.134478334792994;
        m_solPntsLocalCoord3D[6][1]=0.596564995621018;
        m_solPntsLocalCoord3D[6][2]=0.134478334792994;
        m_solPntsLocalCoord3D[7][0]=0.134478334792994;
        m_solPntsLocalCoord3D[7][1]=0.134478334792994;
        m_solPntsLocalCoord3D[7][2]=0.134478334792994;
        m_solPntsLocalCoord3D[8][0]=0.77997600844154;
        m_solPntsLocalCoord3D[8][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[8][2]=0.151831949165937;
        m_solPntsLocalCoord3D[9][0]=0.77997600844154;
        m_solPntsLocalCoord3D[9][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[9][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[10][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[10][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[10][2]=0.77997600844154;
        m_solPntsLocalCoord3D[11][0]=0.151831949165937;
        m_solPntsLocalCoord3D[11][1]=0.77997600844154;
        m_solPntsLocalCoord3D[11][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[12][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[12][1]=0.151831949165937;
        m_solPntsLocalCoord3D[12][2]=0.77997600844154;
        m_solPntsLocalCoord3D[13][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[13][1]=0.77997600844154;
        m_solPntsLocalCoord3D[13][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[14][0]=0.151831949165937;
        m_solPntsLocalCoord3D[14][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[14][2]=0.77997600844154;
        m_solPntsLocalCoord3D[15][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[15][1]=0.151831949165937;
        m_solPntsLocalCoord3D[15][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[16][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[16][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[16][2]=0.151831949165937;
        m_solPntsLocalCoord3D[17][0]=0.0340960211962615;
        m_solPntsLocalCoord3D[17][1]=0.77997600844154;
        m_solPntsLocalCoord3D[17][2]=0.151831949165937;
        m_solPntsLocalCoord3D[18][0]=0.151831949165937;
        m_solPntsLocalCoord3D[18][1]=0.0340960211962615;
        m_solPntsLocalCoord3D[18][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[19][0]=0.77997600844154;
        m_solPntsLocalCoord3D[19][1]=0.151831949165937;
        m_solPntsLocalCoord3D[19][2]=0.0340960211962615;
        m_solPntsLocalCoord3D[20][0]=0.552655643106017;
        m_solPntsLocalCoord3D[20][1]=0.046205150415002;
        m_solPntsLocalCoord3D[20][2]=0.354934056063979;
        m_solPntsLocalCoord3D[21][0]=0.552655643106017;
        m_solPntsLocalCoord3D[21][1]=0.046205150415002;
        m_solPntsLocalCoord3D[21][2]=0.046205150415002;
        m_solPntsLocalCoord3D[22][0]=0.046205150415002;
        m_solPntsLocalCoord3D[22][1]=0.046205150415002;
        m_solPntsLocalCoord3D[22][2]=0.552655643106017;
        m_solPntsLocalCoord3D[23][0]=0.354934056063979;
        m_solPntsLocalCoord3D[23][1]=0.552655643106017;
        m_solPntsLocalCoord3D[23][2]=0.046205150415002;
        m_solPntsLocalCoord3D[24][0]=0.046205150415002;
        m_solPntsLocalCoord3D[24][1]=0.354934056063979;
        m_solPntsLocalCoord3D[24][2]=0.552655643106017;
        m_solPntsLocalCoord3D[25][0]=0.046205150415002;
        m_solPntsLocalCoord3D[25][1]=0.552655643106017;
        m_solPntsLocalCoord3D[25][2]=0.046205150415002;
        m_solPntsLocalCoord3D[26][0]=0.354934056063979;
        m_solPntsLocalCoord3D[26][1]=0.046205150415002;
        m_solPntsLocalCoord3D[26][2]=0.552655643106017;
        m_solPntsLocalCoord3D[27][0]=0.046205150415002;
        m_solPntsLocalCoord3D[27][1]=0.354934056063979;
        m_solPntsLocalCoord3D[27][2]=0.046205150415002;
        m_solPntsLocalCoord3D[28][0]=0.046205150415002;
        m_solPntsLocalCoord3D[28][1]=0.046205150415002;
        m_solPntsLocalCoord3D[28][2]=0.354934056063979;
        m_solPntsLocalCoord3D[29][0]=0.046205150415002;
        m_solPntsLocalCoord3D[29][1]=0.552655643106017;
        m_solPntsLocalCoord3D[29][2]=0.354934056063979;
        m_solPntsLocalCoord3D[30][0]=0.354934056063979;
        m_solPntsLocalCoord3D[30][1]=0.046205150415002;
        m_solPntsLocalCoord3D[30][2]=0.046205150415002;
        m_solPntsLocalCoord3D[31][0]=0.552655643106017;
        m_solPntsLocalCoord3D[31][1]=0.354934056063979;
        m_solPntsLocalCoord3D[31][2]=0.046205150415002;
        m_solPntsLocalCoord3D[32][0]=0.538104322888;
        m_solPntsLocalCoord3D[32][1]=0.228190461068761;
        m_solPntsLocalCoord3D[32][2]=0.00551475497447751;
        m_solPntsLocalCoord3D[33][0]=0.538104322888;
        m_solPntsLocalCoord3D[33][1]=0.228190461068761;
        m_solPntsLocalCoord3D[33][2]=0.228190461068761;
        m_solPntsLocalCoord3D[34][0]=0.228190461068761;
        m_solPntsLocalCoord3D[34][1]=0.228190461068761;
        m_solPntsLocalCoord3D[34][2]=0.538104322888;
        m_solPntsLocalCoord3D[35][0]=0.00551475497447751;
        m_solPntsLocalCoord3D[35][1]=0.538104322888;
        m_solPntsLocalCoord3D[35][2]=0.228190461068761;
        m_solPntsLocalCoord3D[36][0]=0.228190461068761;
        m_solPntsLocalCoord3D[36][1]=0.00551475497447751;
        m_solPntsLocalCoord3D[36][2]=0.538104322888;
        m_solPntsLocalCoord3D[37][0]=0.228190461068761;
        m_solPntsLocalCoord3D[37][1]=0.538104322888;
        m_solPntsLocalCoord3D[37][2]=0.228190461068761;
        m_solPntsLocalCoord3D[38][0]=0.00551475497447751;
        m_solPntsLocalCoord3D[38][1]=0.228190461068761;
        m_solPntsLocalCoord3D[38][2]=0.538104322888;
        m_solPntsLocalCoord3D[39][0]=0.228190461068761;
        m_solPntsLocalCoord3D[39][1]=0.00551475497447751;
        m_solPntsLocalCoord3D[39][2]=0.228190461068761;
        m_solPntsLocalCoord3D[40][0]=0.228190461068761;
        m_solPntsLocalCoord3D[40][1]=0.228190461068761;
        m_solPntsLocalCoord3D[40][2]=0.00551475497447751;
        m_solPntsLocalCoord3D[41][0]=0.228190461068761;
        m_solPntsLocalCoord3D[41][1]=0.538104322888;
        m_solPntsLocalCoord3D[41][2]=0.00551475497447751;
        m_solPntsLocalCoord3D[42][0]=0.00551475497447751;
        m_solPntsLocalCoord3D[42][1]=0.228190461068761;
        m_solPntsLocalCoord3D[42][2]=0.228190461068761;
        m_solPntsLocalCoord3D[43][0]=0.538104322888;
        m_solPntsLocalCoord3D[43][1]=0.00551475497447751;
        m_solPntsLocalCoord3D[43][2]=0.228190461068761;
        m_solPntsLocalCoord3D[44][0]=0.196183759574559;
        m_solPntsLocalCoord3D[44][1]=0.352305260087994;
        m_solPntsLocalCoord3D[44][2]=0.099205720249453;
        m_solPntsLocalCoord3D[45][0]=0.196183759574559;
        m_solPntsLocalCoord3D[45][1]=0.352305260087994;
        m_solPntsLocalCoord3D[45][2]=0.352305260087994;
        m_solPntsLocalCoord3D[46][0]=0.352305260087994;
        m_solPntsLocalCoord3D[46][1]=0.352305260087994;
        m_solPntsLocalCoord3D[46][2]=0.196183759574559;
        m_solPntsLocalCoord3D[47][0]=0.099205720249453;
        m_solPntsLocalCoord3D[47][1]=0.196183759574559;
        m_solPntsLocalCoord3D[47][2]=0.352305260087994;
        m_solPntsLocalCoord3D[48][0]=0.352305260087994;
        m_solPntsLocalCoord3D[48][1]=0.099205720249453;
        m_solPntsLocalCoord3D[48][2]=0.196183759574559;
        m_solPntsLocalCoord3D[49][0]=0.352305260087994;
        m_solPntsLocalCoord3D[49][1]=0.196183759574559;
        m_solPntsLocalCoord3D[49][2]=0.352305260087994;
        m_solPntsLocalCoord3D[50][0]=0.099205720249453;
        m_solPntsLocalCoord3D[50][1]=0.352305260087994;
        m_solPntsLocalCoord3D[50][2]=0.196183759574559;
        m_solPntsLocalCoord3D[51][0]=0.352305260087994;
        m_solPntsLocalCoord3D[51][1]=0.099205720249453;
        m_solPntsLocalCoord3D[51][2]=0.352305260087994;
        m_solPntsLocalCoord3D[52][0]=0.352305260087994;
        m_solPntsLocalCoord3D[52][1]=0.352305260087994;
        m_solPntsLocalCoord3D[52][2]=0.099205720249453;
        m_solPntsLocalCoord3D[53][0]=0.352305260087994;
        m_solPntsLocalCoord3D[53][1]=0.196183759574559;
        m_solPntsLocalCoord3D[53][2]=0.099205720249453;
        m_solPntsLocalCoord3D[54][0]=0.099205720249453;
        m_solPntsLocalCoord3D[54][1]=0.352305260087994;
        m_solPntsLocalCoord3D[54][2]=0.352305260087994;
        m_solPntsLocalCoord3D[55][0]=0.196183759574559;
        m_solPntsLocalCoord3D[55][1]=0.099205720249453;
        m_solPntsLocalCoord3D[55][2]=0.352305260087994;*/
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Quadrature Points not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
    }

  resetFluxReconstructionElementData();
 /* m_shape = CFGeoShape::TETRA;
  m_dimensionality = DIM_3D;
  m_polyOrder = polyOrder;
  m_solPntsLocalCoord1D = solPntDist->getLocalCoords1D(polyOrder);
  m_flxPntsLocalCoord1D = flxPntDist->getLocalCoords1D(polyOrder);
  m_solPntDistribution = solPntDist;
  m_flxPntDistribution = flxPntDist;

  resetFluxReconstructionElementData();*/
}

//////////////////////////////////////////////////////////////////////

TetraFluxReconstructionElementData::~TetraFluxReconstructionElementData()
{
}

//////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////// NOT USED
std::vector<CFreal> TetraFluxReconstructionElementData::getPercentage(CFPolyOrder::Type solOrder){
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

void TetraFluxReconstructionElementData::createFlxPntsLocalCoords()
{
  CFAUTOTRACE;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = m_flxPntsLocalCoord1D.size();
  const CFuint nbrFlxPnts2DTriag = (m_polyOrder+1)*(m_polyOrder+2)/2;
  // set flux point local coordinates
  m_flxPntsLocalCoords.resize(0);

  RealVector flxCoords(3);
  //face 0
  flxCoords[ZTA] = 0;
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2DTriag; ++iFlx)
  {
    flxCoords[KSI] = m_flxPntsLocalCoord2D[iFlx][KSI];
    flxCoords[ETA] = m_flxPntsLocalCoord2D[iFlx][ETA];
    m_flxPntsLocalCoords.push_back(flxCoords);
  }

  //face 1
  flxCoords[ETA] = 0;
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2DTriag; ++iFlx)
  {
    flxCoords[KSI] = m_flxPntsLocalCoord2D[iFlx][ETA];
    flxCoords[ZTA] = m_flxPntsLocalCoord2D[iFlx][KSI];
    m_flxPntsLocalCoords.push_back(flxCoords);
  }

  //face 2 (old 3) aka oblique face
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2DTriag; ++iFlx)
  {
    flxCoords[KSI] = 1.-m_flxPntsLocalCoord2D[iFlx][KSI]-m_flxPntsLocalCoord2D[iFlx][ETA];
    flxCoords[ETA] = m_flxPntsLocalCoord2D[iFlx][ETA];
    flxCoords[ZTA] = m_flxPntsLocalCoord2D[iFlx][KSI];
    m_flxPntsLocalCoords.push_back(flxCoords);
  }

  //face 3 (old 2)
  flxCoords[KSI] = 0;
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2DTriag; ++iFlx)
  {
    flxCoords[ETA] = m_flxPntsLocalCoord2D[iFlx][KSI];
    flxCoords[ZTA] = m_flxPntsLocalCoord2D[iFlx][ETA];
    m_flxPntsLocalCoords.push_back(flxCoords);
  }



  cf_assert(m_flxPntsLocalCoords.size() == (2*(m_polyOrder+1)*(m_polyOrder+2)));
}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createSolPntsLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts3D = m_solPntsLocalCoord3D.size();
  // set solution point local coordinates
  m_solPntsLocalCoords.resize(0);
  for (CFuint iSol = 0; iSol < nbrSolPnts3D; ++iSol)
  {
    RealVector solCoords(3);
    solCoords[KSI] = m_solPntsLocalCoord3D[iSol][KSI];
    solCoords[ETA] = m_solPntsLocalCoord3D[iSol][ETA];
    solCoords[ZTA] = m_solPntsLocalCoord3D[iSol][ZTA];
    m_solPntsLocalCoords.push_back(solCoords);
  }
    cf_assert(m_solPntsLocalCoord3D.size() == (((m_polyOrder+1)*(m_polyOrder+2)*(m_polyOrder+3))/6)); 

}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createFaceFlxPntsFaceLocalCoords()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrFlxPnts2D = m_flxPntsLocalCoord2D.size();

  // set face flux point face local coordinates
  m_faceFlxPntsFaceLocalCoords.resize(0);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx)
  {
    RealVector flxCoord(2);
    flxCoord[KSI] = m_flxPntsLocalCoord2D[iFlx][KSI];
    flxCoord[ETA] = m_flxPntsLocalCoord2D[iFlx][ETA];
    m_faceFlxPntsFaceLocalCoords.push_back(flxCoord);
  }
}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createSolPolyExponents() 
{
  CFAUTOTRACE;
  
  // number of solution points in 2D
  const CFuint nbrSolPnts2D = m_polyOrder+1;//m_flxPntsLocalCoord2D.size();

  // define exponents
  m_solPolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrSolPnts2D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts2D-iKsi; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrSolPnts2D-iKsi-iEta; ++iZta)
      {
        vector< CFint > solPolyExps(3);
        solPolyExps[KSI] = iKsi;
        solPolyExps[ETA] = iEta;
        solPolyExps[ZTA] = iZta;
        m_solPolyExponents.push_back(solPolyExps);
      }
    }
  }/*
  const CFuint polyOrderP1 = m_polyOrder+1;
  for (CFuint iXYZ = 0; iXYZ < polyOrderP1; ++iXYZ)
  {
    for (CFuint iX = 0; iX < iXYZ+1; ++iX)
    {
      for (CFuint iY = 0; iY < iXYZ-iX+1; ++iY)
      {
        vector< CFint > solPolyExps(3);
        solPolyExps[KSI] = iX;
        solPolyExps[ETA] = iY;
        solPolyExps[ZTA] = iXYZ-iX-iY;
        m_solPolyExponents.push_back(solPolyExps);
      }
    }
  }*/
}



//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createNodePolyExponents()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrNodes3D = 4;

  // define exponents
  m_nodePolyExponents.resize(0);

  vector< CFint > nodePolyExps(3);
        nodePolyExps[KSI] = 0.;
        nodePolyExps[ETA] = 0.;
        nodePolyExps[ZTA] = 0.;
        m_nodePolyExponents.push_back(nodePolyExps);

        nodePolyExps[KSI] = 1.;
        nodePolyExps[ETA] = 0.;
        nodePolyExps[ZTA] = 0.;
        m_nodePolyExponents.push_back(nodePolyExps);

        nodePolyExps[KSI] = 0.;
        nodePolyExps[ETA] = 1.;
        nodePolyExps[ZTA] = 0.;
        m_nodePolyExponents.push_back(nodePolyExps);

        nodePolyExps[KSI] = 0.;
        nodePolyExps[ETA] = 0.;
        nodePolyExps[ZTA] = 1.;
        m_nodePolyExponents.push_back(nodePolyExps);            
}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createFaceFluxPntsConn()
{
  CFAUTOTRACE;

  // number of flux points in 2D
  const CFuint nbrFlxPnts2D = (m_polyOrder+1)*(m_polyOrder+2)/2; //on triag face

  // resize m_faceFlxPntConn
  m_faceFlxPntConn.resize(4);

  // variable holding the face index
  CFuint faceIdx = 0;
  CFuint iFlxg = 0; 
  
  // zeroth face (zta=0)
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx, ++iFlxg)
  {
    m_faceFlxPntConn[faceIdx].push_back(iFlxg); 
  }
  ++faceIdx;
  
  // first face (eta=0)
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx, ++iFlxg)
  {
    m_faceFlxPntConn[faceIdx].push_back(iFlxg);
  }
  ++faceIdx;

  // second face (oblique)
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx, ++iFlxg)
  {
    m_faceFlxPntConn[faceIdx].push_back(iFlxg); 
  }
  ++faceIdx;

  // third face (ksi=0)
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx, ++iFlxg)
  {
    m_faceFlxPntConn[faceIdx].push_back(iFlxg); 
  }
  



}

//////////////////////////////////////////////////////////////////////
// !!!! to be checked
void TetraFluxReconstructionElementData::createFaceFluxPntsConnPerOrient()
{
    // number of faces
  const CFuint nbrFaces = 4;

  // number of solution points in 3D
  const CFuint nbrSolPnts3D = m_solPntsLocalCoord3D.size();

  const CFuint nbrFlxPnts2D = m_flxPntsLocalCoord2D.size();

  // number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();

  // flux point indexes for inverted face
  vector< CFuint > invFlxIdxs;
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts2D; ++iKsi)
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

  // number of possible orientations
  const CFuint nbrOrient = 30;

  // create data structure
  m_faceFlxPntConnPerOrient.resize(nbrOrient);
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < nbrFaces; ++iFaceL)
  {
    vector< CFuint > faceFlxConnL = m_faceFlxPntConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrFaces; ++iFaceR)
    {
      vector< CFuint > faceFlxConnR = m_faceFlxPntConn[iFaceR];
      for (CFuint iRot = 0; iRot < 3; ++iRot, ++iOrient)
      {
        m_faceFlxPntConnPerOrient[iOrient].resize(2);
        for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
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
  
  cf_assert(iOrient == nbrOrient);

}



//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createCellNodeCoords()
{
  CFAUTOTRACE;

  m_cellNodeCoords.resize(4);

  // first node
  m_cellNodeCoords[0].resize(3);
  m_cellNodeCoords[0][KSI] = +0.;
  m_cellNodeCoords[0][ETA] = +0.;
  m_cellNodeCoords[0][ZTA] = +0.;

  // second node
  m_cellNodeCoords[1].resize(3);
  m_cellNodeCoords[1][KSI] = +1.0;
  m_cellNodeCoords[1][ETA] = +0.0;
  m_cellNodeCoords[1][ZTA] = +0.0;

  // third node
  m_cellNodeCoords[2].resize(3);
  m_cellNodeCoords[2][KSI] = +0.0;
  m_cellNodeCoords[2][ETA] = +1.0;
  m_cellNodeCoords[2][ZTA] = +0.0;

  // fourth node
  m_cellNodeCoords[3].resize(3);
  m_cellNodeCoords[3][KSI] = +0.0;
  m_cellNodeCoords[3][ETA] = +0.0;
  m_cellNodeCoords[3][ZTA] = +1.0;

}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createFaceNodeConnectivity()
{
  CFAUTOTRACE;

  m_faceNodeConn.resize(4);

  m_faceNodeConn[0].resize(3);
  m_faceNodeConn[0][0] = 0;
  m_faceNodeConn[0][1] = 1;
  m_faceNodeConn[0][2] = 2;

  m_faceNodeConn[1].resize(3);
  m_faceNodeConn[1][0] = 0;
  m_faceNodeConn[1][1] = 3;
  m_faceNodeConn[1][2] = 1;

  m_faceNodeConn[2].resize(3);
  m_faceNodeConn[2][0] = 1;
  m_faceNodeConn[2][1] = 3;
  m_faceNodeConn[2][2] = 2;

  m_faceNodeConn[3].resize(3);
  m_faceNodeConn[3][0] = 0;
  m_faceNodeConn[3][1] = 2;
  m_faceNodeConn[3][2] = 3;

}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createFaceMappedCoordDir()
{
  CFAUTOTRACE;

  m_faceMappedCoordDir.resize(4);

  m_faceMappedCoordDir[0] = 1;
  m_faceMappedCoordDir[1] = 1;
  m_faceMappedCoordDir[2] = 1;
  m_faceMappedCoordDir[3] = 1;
}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createFluxPntFluxDim()
{
  CFAUTOTRACE;

  m_flxPntFlxDim.resize(4 * (m_flxPntsLocalCoord2D.size()));
  
  for (CFuint iFace = 0; iFace < m_faceFlxPntConn.size(); ++iFace)
  {
    for (CFuint iFlx = 0; iFlx < m_faceFlxPntConn[iFace].size(); ++iFlx)
    {
     if (iFace == 0)
      {
        m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 0;
      }
      else if (iFace == 1)
      {
	      m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 1;
      }
      else if (iFace == 2)
      {
	      m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 2;
      }
      else if (iFace == 3)
      {
	      m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 3;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createFaceNormals()
{
  CFAUTOTRACE;

  m_faceNormals.resize(4);
  for (CFuint iFace = 0; iFace < 4; ++iFace)
  {
    m_faceNormals[iFace].resize(3);
  }

  m_faceNormals[0][0] = 0.;
  m_faceNormals[0][1] = 0.;
  m_faceNormals[0][2] = -1.;

  m_faceNormals[1][0] = 0.;
  m_faceNormals[1][1] = -1.;
  m_faceNormals[1][2] = 0.;

  m_faceNormals[3][0] = -1.;
  m_faceNormals[3][1] = 0.;
  m_faceNormals[3][2] = 0.;

  m_faceNormals[2][0] = 1./sqrt(3.);
  m_faceNormals[2][1] = 1./sqrt(3.);
  m_faceNormals[2][2] = 1./sqrt(3.);


}

//////////////////////////////////////////////////////////////////////
// !!!! to be checked

void TetraFluxReconstructionElementData::createFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;

  // number of faces
  const CFuint nbrFaces = 4;


  // number of possible orientations
  const CFuint nbrOrient = 30;


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
  
  for (CFuint iFaceL = 0; iFaceL < nbrFaces; ++iFaceL)
  {
    vector< CFuint > faceNodesL = m_faceNodeConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrFaces; ++iFaceR) //iFaceL
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
  cf_assert(iOrient == nbrOrient);

}


//////////////////////////////////////////////////////////////////////
// !!!! to be checked

void TetraFluxReconstructionElementData::createFaceIntegrationCoefs()
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
// !!!! to be updated for higher order (P>3)

void TetraFluxReconstructionElementData::createCellAvgSolCoefs()
{
  CFAUTOTRACE;

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // resize m_cellAvgSolCoefs
  m_cellAvgSolCoefs.resize(nbrSolPnts);

  switch(m_polyOrder)
  {
    case CFPolyOrder::ORDER0:
    {
      m_cellAvgSolCoefs[0] = 1.333333333333333;
    } break;
    case CFPolyOrder::ORDER1:
    {
      m_cellAvgSolCoefs[0] = 0.333333333333333;
      m_cellAvgSolCoefs[1] = 0.333333333333333;
      m_cellAvgSolCoefs[2] = 0.333333333333333;
      m_cellAvgSolCoefs[3] = 0.333333333333333;

    } break;
    case CFPolyOrder::ORDER2:
    {
      m_cellAvgSolCoefs[0] = 0.063510846457611867;
      m_cellAvgSolCoefs[1] = 0.063510846457611867;
      m_cellAvgSolCoefs[2] = 0.063510846457611867;
      m_cellAvgSolCoefs[3] = 0.063510846457611867;
      m_cellAvgSolCoefs[4] = 0.179881657917148;
      m_cellAvgSolCoefs[5] = 0.179881657917148;
      m_cellAvgSolCoefs[6] = 0.179881657917148;
      m_cellAvgSolCoefs[7] = 0.179881657917148;
      m_cellAvgSolCoefs[8] = 0.179881657917148;
      m_cellAvgSolCoefs[9] = 0.179881657917148;

    } break;
    case CFPolyOrder::ORDER3:
    {
      m_cellAvgSolCoefs[0] = 0.009422766392626;
      m_cellAvgSolCoefs[1] = 0.009422766392626;
      m_cellAvgSolCoefs[2] = 0.009422766392626;
      m_cellAvgSolCoefs[3] = 0.009422766392626;

      m_cellAvgSolCoefs[4] = 0.062664891962516933;
      m_cellAvgSolCoefs[5] = 0.062664891962516933;
      m_cellAvgSolCoefs[6] = 0.062664891962516933;
      m_cellAvgSolCoefs[7] = 0.062664891962516933;
      m_cellAvgSolCoefs[8] = 0.062664891962516933;
      m_cellAvgSolCoefs[9] = 0.062664891962516933;
      m_cellAvgSolCoefs[10] = 0.062664891962516933;
      m_cellAvgSolCoefs[11] = 0.062664891962516933;
      m_cellAvgSolCoefs[12] = 0.062664891962516933;
      m_cellAvgSolCoefs[13] = 0.062664891962516933;
      m_cellAvgSolCoefs[14] = 0.062664891962516933;
      m_cellAvgSolCoefs[15] = 0.062664891962516933;

      m_cellAvgSolCoefs[16] = 0.13591589105315733;
      m_cellAvgSolCoefs[17] = 0.13591589105315733;
      m_cellAvgSolCoefs[18] = 0.13591589105315733;
      m_cellAvgSolCoefs[19] = 0.13591589105315733;
      
    } break;
        case CFPolyOrder::ORDER4:
    {
      m_cellAvgSolCoefs[0]=0.124232764159378;
      m_cellAvgSolCoefs[1]=0.00292006186205173;
      m_cellAvgSolCoefs[2]=0.00292006186205173;
      m_cellAvgSolCoefs[3]=0.00292006186205173;
      m_cellAvgSolCoefs[4]=0.00292006186205173;
      m_cellAvgSolCoefs[5]=0.0333740527582328;
      m_cellAvgSolCoefs[6]=0.0333740527582328;
      m_cellAvgSolCoefs[7]=0.0333740527582328;
      m_cellAvgSolCoefs[8]=0.0333740527582328;
      m_cellAvgSolCoefs[9]=0.0333740527582328;
      m_cellAvgSolCoefs[10]=0.0333740527582328;
      m_cellAvgSolCoefs[11]=0.0191194226903553;
      m_cellAvgSolCoefs[12]=0.0191194226903553;
      m_cellAvgSolCoefs[13]=0.0191194226903553;
      m_cellAvgSolCoefs[14]=0.0191194226903553;
      m_cellAvgSolCoefs[15]=0.0191194226903553;
      m_cellAvgSolCoefs[16]=0.0191194226903553;
      m_cellAvgSolCoefs[17]=0.0191194226903553;
      m_cellAvgSolCoefs[18]=0.0191194226903553;
      m_cellAvgSolCoefs[19]=0.0191194226903553;
      m_cellAvgSolCoefs[20]=0.0191194226903553;
      m_cellAvgSolCoefs[21]=0.0191194226903553;
      m_cellAvgSolCoefs[22]=0.0191194226903553;
      m_cellAvgSolCoefs[23]=0.0639785777410072;
      m_cellAvgSolCoefs[24]=0.0639785777410072;
      m_cellAvgSolCoefs[25]=0.0639785777410072;
      m_cellAvgSolCoefs[26]=0.0639785777410072;
      m_cellAvgSolCoefs[27]=0.0639785777410072;
      m_cellAvgSolCoefs[28]=0.0639785777410072;
      m_cellAvgSolCoefs[29]=0.0639785777410072;
      m_cellAvgSolCoefs[30]=0.0639785777410072;
      m_cellAvgSolCoefs[31]=0.0639785777410072;
      m_cellAvgSolCoefs[32]=0.0639785777410072;
      m_cellAvgSolCoefs[33]=0.0639785777410072;
      m_cellAvgSolCoefs[34]=0.0639785777410072; 
    } break;
        case CFPolyOrder::ORDER5:
    {
      m_cellAvgSolCoefs[0]=0.00138308164481866;
      m_cellAvgSolCoefs[1]=0.00138308164481866;
      m_cellAvgSolCoefs[2]=0.00138308164481866;
      m_cellAvgSolCoefs[3]=0.00138308164481866;
      m_cellAvgSolCoefs[4]=0.0488388488540144;
      m_cellAvgSolCoefs[5]=0.0488388488540144;
      m_cellAvgSolCoefs[6]=0.0488388488540144;
      m_cellAvgSolCoefs[7]=0.0488388488540144;
      m_cellAvgSolCoefs[8]=0.0128022193865973;
      m_cellAvgSolCoefs[9]=0.0128022193865973;
      m_cellAvgSolCoefs[10]=0.0128022193865973;
      m_cellAvgSolCoefs[11]=0.0128022193865973;
      m_cellAvgSolCoefs[12]=0.0128022193865973;
      m_cellAvgSolCoefs[13]=0.0128022193865973;
      m_cellAvgSolCoefs[14]=0.0128022193865973;
      m_cellAvgSolCoefs[15]=0.0128022193865973;
      m_cellAvgSolCoefs[16]=0.0128022193865973;
      m_cellAvgSolCoefs[17]=0.0128022193865973;
      m_cellAvgSolCoefs[18]=0.0128022193865973;
      m_cellAvgSolCoefs[19]=0.0128022193865973;
      m_cellAvgSolCoefs[20]=0.0219325302397642;
      m_cellAvgSolCoefs[21]=0.0219325302397642;
      m_cellAvgSolCoefs[22]=0.0219325302397642;
      m_cellAvgSolCoefs[23]=0.0219325302397642;
      m_cellAvgSolCoefs[24]=0.0219325302397642;
      m_cellAvgSolCoefs[25]=0.0219325302397642;
      m_cellAvgSolCoefs[26]=0.0219325302397642;
      m_cellAvgSolCoefs[27]=0.0219325302397642;
      m_cellAvgSolCoefs[28]=0.0219325302397642;
      m_cellAvgSolCoefs[29]=0.0219325302397642;
      m_cellAvgSolCoefs[30]=0.0219325302397642;
      m_cellAvgSolCoefs[31]=0.0219325302397642;
      m_cellAvgSolCoefs[32]=0.0204997022017746;
      m_cellAvgSolCoefs[33]=0.0204997022017746;
      m_cellAvgSolCoefs[34]=0.0204997022017746;
      m_cellAvgSolCoefs[35]=0.0204997022017746;
      m_cellAvgSolCoefs[36]=0.0204997022017746;
      m_cellAvgSolCoefs[37]=0.0204997022017746;
      m_cellAvgSolCoefs[38]=0.0204997022017746;
      m_cellAvgSolCoefs[39]=0.0204997022017746;
      m_cellAvgSolCoefs[40]=0.0204997022017746;
      m_cellAvgSolCoefs[41]=0.0204997022017746;
      m_cellAvgSolCoefs[42]=0.0204997022017746;
      m_cellAvgSolCoefs[43]=0.0204997022017746;
      m_cellAvgSolCoefs[44]=0.039136015783364;
      m_cellAvgSolCoefs[45]=0.039136015783364;
      m_cellAvgSolCoefs[46]=0.039136015783364;
      m_cellAvgSolCoefs[47]=0.039136015783364;
      m_cellAvgSolCoefs[48]=0.039136015783364;
      m_cellAvgSolCoefs[49]=0.039136015783364;
      m_cellAvgSolCoefs[50]=0.039136015783364;
      m_cellAvgSolCoefs[51]=0.039136015783364;
      m_cellAvgSolCoefs[52]=0.039136015783364;
      m_cellAvgSolCoefs[53]=0.039136015783364;
      m_cellAvgSolCoefs[54]=0.039136015783364;
      m_cellAvgSolCoefs[55]=0.039136015783364;
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Quadrature Weights not implemented for order "
                                    + StringOps::to_str(m_polyOrder) + ".");
    }
  }
}

//////////////////////////////////////////////////////////////////////
// !!!!!!! maybe not used
void TetraFluxReconstructionElementData::createCellCenterDerivCoefs()
{
  CFAUTOTRACE;

  // center coordinate
  vector< RealVector > centerCoord(1,RealVector(3));
  centerCoord[0][KSI] = 0.33;
  centerCoord[0][ETA] = 0.33;
  centerCoord[0][ZTA] = 0.33;

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

void TetraFluxReconstructionElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
                                                        vector< RealVector >& nodalSet)
{
  
  const CFuint nbrSolPnts3D = ((order+1)*(order+2)*(order+3))/6;

  // set solution point local coordinates
  nodalSet.resize(0);
  for (CFuint iSol = 0; iSol < nbrSolPnts3D; ++iSol)
  {
    RealVector node(3);
    node[KSI] = m_solPntsLocalCoord3D[iSol][KSI];
    node[ETA] = m_solPntsLocalCoord3D[iSol][ETA];
    node[ZTA] = m_solPntsLocalCoord3D[iSol][ZTA];
    nodalSet.push_back(node);
  }
}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::setCFLConvDiffRatio()
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
// !!!! NOT USED!! !!!!

void TetraFluxReconstructionElementData::createFaceOutputPntCellMappedCoords()
{
  // number of points on a face
  const CFuint nbrFacePnts2D = m_polyOrder == 0 ? 2 :(m_polyOrder+ 1);

  // face mapped coordinates of uniform distribution of points
  vector<RealVector> faceMapCoords;
  const CFreal dKsiEta = 0 ? 2.0 : 2.0/m_polyOrder;
  CFreal ksi = 0.0;
  m_faceOutputPntFaceMappedCoords.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrFacePnts2D; ++iKsi, ksi += dKsiEta)
  {
    CFreal eta = 0.0;
    for (CFuint iEta = 0; iEta < nbrFacePnts2D; ++iEta, eta += dKsiEta)
    {
      RealVector mapCoord(2);
      mapCoord[KSI] = ksi;
      mapCoord[ETA] = eta;
      m_faceOutputPntFaceMappedCoords.push_back(mapCoord);
    }
  }
  const CFuint nbrFacePnts = m_faceOutputPntFaceMappedCoords.size();

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


      m_faceOutputPntCellMappedCoords[iFace].push_back(fun0*faceNodeCoords[0]+
                                                       fun1*faceNodeCoords[1]+
                                                       fun2*faceNodeCoords[2]);
      
    }
  }
  
 //cout << "createFaceOutputPntCellMappedCoords -- not done yet"<<endl;
}

//////////////////////////////////////////////////////////////////////
// !!!! maybe not used
void TetraFluxReconstructionElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;

  // number of nodes 1D
  /*const CFuint nbrNodes1D = m_polyOrder + 1;

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
  }*/
  //cout << "createFaceOutputPntConn -- not done yet"<<endl;
}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createVandermondeMatrix()
{
  CFAUTOTRACE;
  
  const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  
  CFuint order;
  CFuint idx;
  
  m_vandermonde.resize(nbrSolPnts,nbrSolPnts);
  m_vandermondeInv.resize(nbrSolPnts,nbrSolPnts);
  
  if(m_polyOrder != CFPolyOrder::ORDER0 && m_polyOrder != CFPolyOrder::ORDER1)
  {
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
     CFuint modalDof = 0;
    
        for (CFuint iOrderKsi = 0; iOrderKsi < (m_polyOrder+1); ++iOrderKsi)
        {
          for (CFuint iOrderEta = 0; iOrderEta < (m_polyOrder+1)-iOrderKsi; ++iOrderEta)
          {
            for (CFuint iOrderZta = 0; iOrderZta < (m_polyOrder+1)-iOrderKsi-iOrderEta; ++iOrderZta)
            {
              double a = ((2.*m_solPntsLocalCoords[iSol][KSI])/(1.-m_solPntsLocalCoords[iSol][ETA]-m_solPntsLocalCoords[iSol][ZTA]))-1.;
              double b = ((2.*m_solPntsLocalCoords[iSol][ETA])/(1.-m_solPntsLocalCoords[iSol][ZTA]))-1.;    
              double c = (2.*m_solPntsLocalCoords[iSol][ZTA])-1.;

              double h1 = ComputeJacobi(iOrderKsi, 0., 0., a);
              double h2 = ComputeJacobi(iOrderEta, 2.*iOrderKsi+1., 0., b);
              double h3 = ComputeJacobi(iOrderZta, 2.*iOrderKsi+2.*iOrderEta+2., 0., c);
              m_vandermonde(iSol,modalDof)=sqrt(8.0)*0.25*h1*h2*h3*pow((1.-b),iOrderKsi)*pow((1.-c),iOrderKsi+iOrderEta);
              modalDof+=1;
            }
          } 
        }
    }
    InvertMatrix(m_vandermonde,m_vandermondeInv);
  }
}

//////////////////////////////////////////////////////////////////////

void TetraFluxReconstructionElementData::createFlxSolDependencies()
{
  CFAUTOTRACE;

  const CFuint nbrSolPnts = m_solPntsLocalCoord3D.size();
  const CFuint nbrFlxPnts = 4*(m_flxPntsLocalCoord2D.size());

  m_solFlxDep.resize(nbrSolPnts);
  m_solSolDep.resize(nbrSolPnts);  // also push back the current 
  m_flxSolDep.resize(nbrFlxPnts);

  CFuint iSol = 0;
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx){
      for(CFuint iSol= 0; iSol< nbrSolPnts; ++iSol){
        m_flxSolDep[iFlx].push_back(iSol);
      }
    }
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol){
      for(CFuint iFlx= 0; iFlx< nbrFlxPnts; ++iFlx){
        m_solFlxDep[iSol].push_back(iFlx);
      }
    }
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol){
      for(CFuint iSol1= 0; iSol1< nbrSolPnts; ++iSol1){
        m_solSolDep[iSol].push_back(iSol1);
      }
    }
}

//////////////////////////////////////////////////////////////////////////////
CFreal TetraFluxReconstructionElementData::ComputeJacobi(CFuint N,CFuint alpha, CFuint beta,CFreal x)
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
CFreal TetraFluxReconstructionElementData::factorial(CFreal n)
{
  return (n==1. || n==0.) ? 1. : factorial(n-1.)*n;
}

//////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
