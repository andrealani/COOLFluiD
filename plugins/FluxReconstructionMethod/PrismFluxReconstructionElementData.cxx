#include <fstream>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/PrismFluxReconstructionElementData.hh"
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
                        
                        solPntsLocalCoord2D[1][0] =  0.816847572980440; 
                        solPntsLocalCoord2D[1][1] =  0.091576213509780;  
                        
                        solPntsLocalCoord2D[2][0] =  0.091576213509780; 
                        solPntsLocalCoord2D[2][1] =  0.816847572980440;
                        
                        solPntsLocalCoord2D[3][0] =  0.445948490915964; 
                        solPntsLocalCoord2D[3][1] =  0.108103018168071; 
                        
                        solPntsLocalCoord2D[4][0] =  0.108103018168071; 
                        solPntsLocalCoord2D[4][1] =  0.445948490915964;
                        
                        solPntsLocalCoord2D[5][0] =  0.445948490915964; 
                        solPntsLocalCoord2D[5][1] =  0.445948490915964;


      } break;

      case CFPolyOrder::ORDER3:          
      {
                        solPntsLocalCoord2D[0][0] = 0.055564052669793;
                        solPntsLocalCoord2D[0][1] = 0.055564052669793;
                        
                        solPntsLocalCoord2D[1][0] = 0.888871894660413; 
                        solPntsLocalCoord2D[1][1] = 0.055564052669793;	
                        
                        solPntsLocalCoord2D[2][0] = 0.055564052669793; 
                        solPntsLocalCoord2D[2][1] = 0.888871894660413;  	

                        
                        solPntsLocalCoord2D[3][0] = 0.634210747745723; 
                        solPntsLocalCoord2D[3][1] = 0.070255540518384;  	

                        
                        solPntsLocalCoord2D[4][0] = 0.070255540518384; 
                        solPntsLocalCoord2D[4][1] = 0.634210747745723;	
                        
                        solPntsLocalCoord2D[5][0] = 0.295533711735893; 
                        solPntsLocalCoord2D[5][1] = 0.634210747745723; 	

                        
                        solPntsLocalCoord2D[6][0] = 0.295533711735893; 
                        solPntsLocalCoord2D[6][1] = 0.070255540518384;  	
                        
                        solPntsLocalCoord2D[7][0] = 0.070255540518384; 
                        solPntsLocalCoord2D[7][1] = 0.295533711735893; 	
                        
                        solPntsLocalCoord2D[8][0] = 0.634210747745723; 
                        solPntsLocalCoord2D[8][1] = 0.295533711735893; 	
                        
                        solPntsLocalCoord2D[9][0] = 0.333333333333333; 
                        solPntsLocalCoord2D[9][1] = 0.333333333333333;


      } break;
      case CFPolyOrder::ORDER4:   
      {
                        solPntsLocalCoord2D[0][0] = 0.035870877695734; 
                        solPntsLocalCoord2D[0][1] = 0.035870877695734;

                        solPntsLocalCoord2D[1][0] = 0.928258244608533; 
                        solPntsLocalCoord2D[1][1] = 0.035870877695734;

                        solPntsLocalCoord2D[2][0] = 0.035870877695734; 
                        solPntsLocalCoord2D[2][1] = 0.928258244608533; 

                        solPntsLocalCoord2D[3][0] = 0.241729395767967;
                        solPntsLocalCoord2D[3][1] = 0.241729395767967;

                        solPntsLocalCoord2D[4][0] = 0.516541208464066; 
                        solPntsLocalCoord2D[4][1] = 0.241729395767967; 

                        solPntsLocalCoord2D[5][0] = 0.241729395767967; 
                        solPntsLocalCoord2D[5][1] = 0.516541208464066; 

                        solPntsLocalCoord2D[6][0] = 0.474308787777079;
                        solPntsLocalCoord2D[6][1] = 0.051382424445843; 

                        solPntsLocalCoord2D[7][0] = 0.051382424445843;
                        solPntsLocalCoord2D[7][1] = 0.474308787777079;

                        solPntsLocalCoord2D[8][0] = 0.474308787777079; 
                        solPntsLocalCoord2D[8][1] = 0.474308787777079;

                        solPntsLocalCoord2D[9][0] = 0.751183631106484; 
                        solPntsLocalCoord2D[9][1] = 0.047312487011716;

                        solPntsLocalCoord2D[10][0] = 0.047312487011716; 
                        solPntsLocalCoord2D[10][1] = 0.751183631106484;

                        solPntsLocalCoord2D[11][0] =  0.201503881881800;
                        solPntsLocalCoord2D[11][1] =  0.751183631106484;

                        solPntsLocalCoord2D[12][0] = 0.201503881881800; 
                        solPntsLocalCoord2D[12][1] =  0.047312487011716;

                        solPntsLocalCoord2D[13][0] = 0.047312487011716;
                        solPntsLocalCoord2D[13][1] = 0.201503881881800;

                        solPntsLocalCoord2D[14][0] = 0.751183631106484; 
                        solPntsLocalCoord2D[14][1] = 0.201503881881800;

      } break;

      case CFPolyOrder::ORDER5:
      {
                        solPntsLocalCoord2D[0][0] = 0.028112952182664; 
                        solPntsLocalCoord2D[0][1] = 0.028112952182664; 

                        solPntsLocalCoord2D[1][0] = 0.943774095634672;
                        solPntsLocalCoord2D[1][1] = 0.028112952182664; 

                        solPntsLocalCoord2D[2][0] = 0.028112952182664;
                        solPntsLocalCoord2D[2][1] = 0.943774095634672;

                        solPntsLocalCoord2D[3][0] = 0.177139098469317;
                        solPntsLocalCoord2D[3][1] = 0.177139098469317;

                        solPntsLocalCoord2D[4][0] = 0.645721803061365; 
                        solPntsLocalCoord2D[4][1] = 0.177139098469317;

                        solPntsLocalCoord2D[5][0] = 0.177139098469317;
                        solPntsLocalCoord2D[5][1] = 0.645721803061365;

                        solPntsLocalCoord2D[6][0] = 0.405508595867433;
                        solPntsLocalCoord2D[6][1] = 0.188982808265134;

                        solPntsLocalCoord2D[7][0] = 0.188982808265134;
                        solPntsLocalCoord2D[7][1] = 0.405508595867433;

                        solPntsLocalCoord2D[8][0] = 0.405508595867433;
                        solPntsLocalCoord2D[8][1] = 0.405508595867433;

                        solPntsLocalCoord2D[9][0] = 0.817900980028499; 
                        solPntsLocalCoord2D[9][1] = 0.033533207700614;

                        solPntsLocalCoord2D[10][0] = 0.033533207700614; 
                        solPntsLocalCoord2D[10][1] = 0.817900980028499;

                        solPntsLocalCoord2D[11][0] = 0.148565812270887;
                        solPntsLocalCoord2D[11][1] = 0.817900980028499;

                        solPntsLocalCoord2D[12][0] = 0.148565812270887; 
                        solPntsLocalCoord2D[12][1] = 0.033533207700614; 

                        solPntsLocalCoord2D[13][0] = 0.033533207700614;
                        solPntsLocalCoord2D[13][1] = 0.148565812270887;

                        solPntsLocalCoord2D[14][0] = 0.817900980028499;
                        solPntsLocalCoord2D[14][1] = 0.148565812270887;

                        solPntsLocalCoord2D[15][0] = 0.604978911775132; 
                        solPntsLocalCoord2D[15][1] = 0.037824789609186;

                        solPntsLocalCoord2D[16][0] = 0.037824789609186; 
                        solPntsLocalCoord2D[16][1] = 0.604978911775132;

                        solPntsLocalCoord2D[17][0] = 0.357196298615681;
                        solPntsLocalCoord2D[17][1] = 0.604978911775132;

                        solPntsLocalCoord2D[18][0] = 0.357196298615681; 
                        solPntsLocalCoord2D[18][1] = 0.037824789609186;

                        solPntsLocalCoord2D[19][0] = 0.037824789609186; 
                        solPntsLocalCoord2D[19][1] = 0.357196298615681; 

                        solPntsLocalCoord2D[20][0] = 0.604978911775132; 
                        solPntsLocalCoord2D[20][1] = 0.357196298615681; 

      } break;


   case CFPolyOrder::ORDER6:
      {  
                        solPntsLocalCoord2D[0][0] =  0.0000000000000;
                        solPntsLocalCoord2D[0][1] =  0.9451704450174; 
                        solPntsLocalCoord2D[1][0] =  0.9451704450173; 
                        solPntsLocalCoord2D[1][1] =  0.0000000000000; 
                        solPntsLocalCoord2D[2][0] =  0.9289002405719;
                        solPntsLocalCoord2D[2][1] =  0.0685505797224;
                        solPntsLocalCoord2D[3][0] =  0.0685505797224;
                        solPntsLocalCoord2D[3][1] =  0.9289002405717;
                        solPntsLocalCoord2D[4][0] =  0.0243268355615;
                        solPntsLocalCoord2D[4][1] =  0.0243268355616; 
                        solPntsLocalCoord2D[5][0] =  0.1279662835335;
                        solPntsLocalCoord2D[5][1] =  0.0277838749488;
                        solPntsLocalCoord2D[6][0] =  0.0277838749488; 
                        solPntsLocalCoord2D[6][1] =  0.1279662835337; 
                        solPntsLocalCoord2D[7][0] =  0.0287083428360; 
                        solPntsLocalCoord2D[7][1] =  0.7498347588657;
                        solPntsLocalCoord2D[8][0] =  0.7498347588656;
                        solPntsLocalCoord2D[8][1] =  0.0287083428360;
                        solPntsLocalCoord2D[9][0] =  0.7228007909707;
                        solPntsLocalCoord2D[9][1] =  0.2497602062385;
                        solPntsLocalCoord2D[10][0] = 0.2497602062386; 
                        solPntsLocalCoord2D[10][1] = 0.7228007909707;
                        solPntsLocalCoord2D[11][0] = 0.0865562992839;
                        solPntsLocalCoord2D[11][1] = 0.8325513856997;
                        solPntsLocalCoord2D[12][0] = 0.8325513856998;
                        solPntsLocalCoord2D[12][1] = 0.0865562992839;
                        solPntsLocalCoord2D[13][0] = 0.3061619157672;
                        solPntsLocalCoord2D[13][1] = 0.0303526617491;
                        solPntsLocalCoord2D[14][0] = 0.0303526617491;
                        solPntsLocalCoord2D[14][1] = 0.3061619157675;
                        solPntsLocalCoord2D[15][0] = 0.4868610595047;
                        solPntsLocalCoord2D[15][1] = 0.4868610595047;
                        solPntsLocalCoord2D[16][0] = 0.6657904293017;
                        solPntsLocalCoord2D[16][1] = 0.1765456154219;
                        solPntsLocalCoord2D[17][0] = 0.1765456154221; 
                        solPntsLocalCoord2D[17][1] = 0.6657904293016;
                        solPntsLocalCoord2D[18][0] = 0.0293121007360; 
                        solPntsLocalCoord2D[18][1] = 0.5295657488669;
                        solPntsLocalCoord2D[19][0] = 0.5295657488667;
                        solPntsLocalCoord2D[19][1] = 0.0293121007360;
                        solPntsLocalCoord2D[20][0] = 0.1444673824391; 
                        solPntsLocalCoord2D[20][1] = 0.1444673824391;
                        solPntsLocalCoord2D[21][0] = 0.3299740111411; 
                        solPntsLocalCoord2D[21][1] = 0.5361815729050;
                        solPntsLocalCoord2D[22][0] = 0.5361815729052; 
                        solPntsLocalCoord2D[22][1] = 0.3299740111409;
                        solPntsLocalCoord2D[23][0] = 0.5511507516862;
                        solPntsLocalCoord2D[23][1] = 0.1437790861923;
                        solPntsLocalCoord2D[24][0] = 0.1437790861923;
                        solPntsLocalCoord2D[24][1] = 0.5511507516862;
                        solPntsLocalCoord2D[25][0] = 0.3348066587327; 
                        solPntsLocalCoord2D[25][1] = 0.1529619437161;
                        solPntsLocalCoord2D[26][0] = 0.1529619437161; 
                        solPntsLocalCoord2D[26][1] = 0.3348066587327;
                        solPntsLocalCoord2D[27][0] = 0.3430183498147; 
                        solPntsLocalCoord2D[27][1] = 0.3430183498147;
                      
      } break;


      case CFPolyOrder::ORDER7:
      {
                        solPntsLocalCoord2D[0][0] = 0.0242935351590;
                        solPntsLocalCoord2D[0][1] = 0.9493059293846; 
                        solPntsLocalCoord2D[1][0] = 0.0265193427722;
                        solPntsLocalCoord2D[1][1] = 0.0242695130640; 
                        solPntsLocalCoord2D[2][0] = 0.9492126023551;
                        solPntsLocalCoord2D[2][1] = 0.0265067966437;
                        solPntsLocalCoord2D[3][0] = 0.0033775763749;
                        solPntsLocalCoord2D[3][1] = 0.4767316412363;
                        solPntsLocalCoord2D[4][0] = 0.4757672298101;
                        solPntsLocalCoord2D[4][1] = 0.5198921829102; 
                        solPntsLocalCoord2D[5][0] = 0.5190783193471;
                        solPntsLocalCoord2D[5][1] = 0.0055912706202;
                        solPntsLocalCoord2D[6][0] = 0.8616839745321; 
                        solPntsLocalCoord2D[6][1] = 0.0133996048618; 
                        solPntsLocalCoord2D[7][0] = 0.1249209759926;
                        solPntsLocalCoord2D[7][1] = 0.8613054321334;
                        solPntsLocalCoord2D[8][0] = 0.0138565453861;
                        solPntsLocalCoord2D[8][1] = 0.1247733717358;
                        solPntsLocalCoord2D[9][0] = 0.0211887064222;
                        solPntsLocalCoord2D[9][1] = 0.8438438351223;
                        solPntsLocalCoord2D[10][0]  =0.8432296787219; 
                        solPntsLocalCoord2D[10][1] = 0.1354563645830;
                        solPntsLocalCoord2D[11][0] = 0.1354231797865;
                        solPntsLocalCoord2D[11][1] = 0.0213482820656;
                        solPntsLocalCoord2D[12][0] = 0.3088853510679;
                        solPntsLocalCoord2D[12][1] = 0.0221919663014;
                        solPntsLocalCoord2D[13][0] = 0.6685057595169;
                        solPntsLocalCoord2D[13][1] = 0.3089012879389;
                        solPntsLocalCoord2D[14][0] = 0.0226545012557;
                        solPntsLocalCoord2D[14][1] = 0.6691709943321;
                        solPntsLocalCoord2D[15][0] = 0.2808515408772;
                        solPntsLocalCoord2D[15][1] = 0.6924718155106;
                        solPntsLocalCoord2D[16][0] = 0.6922446749051;
                        solPntsLocalCoord2D[16][1] = 0.0268723345026;
                        solPntsLocalCoord2D[17][0] = 0.0268617447119; 
                        solPntsLocalCoord2D[17][1] = 0.2810093973222;
                        solPntsLocalCoord2D[18][0] = 0.1141778485470; 
                        solPntsLocalCoord2D[18][1] = 0.7973581413586;
                        solPntsLocalCoord2D[19][0] = 0.7974807922061;
                        solPntsLocalCoord2D[19][1] = 0.0879806508791;
                        solPntsLocalCoord2D[20][0] = 0.0892807293894; 
                        solPntsLocalCoord2D[20][1] = 0.1145020561128;
                        solPntsLocalCoord2D[21][0] = 0.1052487892455; 
                        solPntsLocalCoord2D[21][1] = 0.6686904119922;
                        solPntsLocalCoord2D[22][0] = 0.6663022280740; 
                        solPntsLocalCoord2D[22][1] = 0.2275051631832;
                        solPntsLocalCoord2D[23][0] = 0.2307803737547;
                        solPntsLocalCoord2D[23][1] = 0.1054572561221;
                        solPntsLocalCoord2D[24][0] = 0.1705059157540;
                        solPntsLocalCoord2D[24][1] = 0.5174064398658;
                        solPntsLocalCoord2D[25][0] = 0.5086593973043; 
                        solPntsLocalCoord2D[25][1] = 0.3170523855209;
                        solPntsLocalCoord2D[26][0] = 0.3141823862281; 
                        solPntsLocalCoord2D[26][1] = 0.1810706361659;
                        solPntsLocalCoord2D[27][0] = 0.4617460817864; 
                        solPntsLocalCoord2D[27][1] = 0.4678594539804;
                        solPntsLocalCoord2D[28][0] = 0.0693087496081; 
                        solPntsLocalCoord2D[28][1] = 0.4622856042085;
                        solPntsLocalCoord2D[29][0] = 0.4622856042085; 
                        solPntsLocalCoord2D[29][1] = 0.0724357805669;
                        solPntsLocalCoord2D[30][0] = 0.2578625857893;
                        solPntsLocalCoord2D[30][1] = 0.6131395039177;
                        solPntsLocalCoord2D[31][0] = 0.6112627766779;
                        solPntsLocalCoord2D[31][1] = 0.1300360834609;
                        solPntsLocalCoord2D[32][0] = 0.1305182135934;
                        solPntsLocalCoord2D[32][1] = 0.2581713828884;
                        solPntsLocalCoord2D[33][0] = 0.4281437991828; 
                        solPntsLocalCoord2D[33][1] = 0.2362005969817;
                        solPntsLocalCoord2D[34][0] = 0.3356995783730; 
                        solPntsLocalCoord2D[34][1] = 0.4311026308588;
                        solPntsLocalCoord2D[35][0] = 0.2305424298836; 
                        solPntsLocalCoord2D[35][1] = 0.3456013949376;

                      
      } break;

        case CFPolyOrder::ORDER8:
      {
                        solPntsLocalCoord2D[0][0] = 0.0000000000000;
                        solPntsLocalCoord2D[0][1] = 1.0000000000000; 
                        solPntsLocalCoord2D[1][0] = 1.0000000000000;
                        solPntsLocalCoord2D[1][1] = 0.0000000000000; 
                        solPntsLocalCoord2D[2][0] = 0.0000000000000;
                        solPntsLocalCoord2D[2][1] = 0.0000000000000;
                        solPntsLocalCoord2D[3][0] = 0.0573330873026;
                        solPntsLocalCoord2D[3][1] = 0.0151382269814;
                        solPntsLocalCoord2D[4][0] = 0.0573330873026;
                        solPntsLocalCoord2D[4][1] = 0.9275286857160; 
                        solPntsLocalCoord2D[5][0] = 0.9275286857160;
                        solPntsLocalCoord2D[5][1] = 0.0573330873026;
                        solPntsLocalCoord2D[6][0] = 0.0151382269814; 
                        solPntsLocalCoord2D[6][1] = 0.0573330873026; 
                        solPntsLocalCoord2D[7][0] = 0.9275286857160;
                        solPntsLocalCoord2D[7][1] = 0.0151382269814;
                        solPntsLocalCoord2D[8][0] = 0.0151382269814;
                        solPntsLocalCoord2D[8][1] = 0.9275286857160;
                        solPntsLocalCoord2D[9][0] = 0.8159625040711;
                        solPntsLocalCoord2D[9][1] = 0.1659719969565;
                        solPntsLocalCoord2D[10][0]  =0.8159625040711; 
                        solPntsLocalCoord2D[10][1] = 0.0180654989724;
                        solPntsLocalCoord2D[11][0] = 0.1659719969565;
                        solPntsLocalCoord2D[11][1] = 0.8159625040711;
                        solPntsLocalCoord2D[12][0] = 0.0180654989724;
                        solPntsLocalCoord2D[12][1] = 0.8159625040711;
                        solPntsLocalCoord2D[13][0] = 0.1659719969565;
                        solPntsLocalCoord2D[13][1] = 0.0180654989724;
                        solPntsLocalCoord2D[14][0] = 0.0180654989724;
                        solPntsLocalCoord2D[14][1] = 0.1659719969565;
                        solPntsLocalCoord2D[15][0] = 0.3165475556378;
                        solPntsLocalCoord2D[15][1] = 0.0186886898773;
                        solPntsLocalCoord2D[16][0] = 0.6647637544849;
                        solPntsLocalCoord2D[16][1] = 0.0186886898773;
                        solPntsLocalCoord2D[17][0] = 0.0186886898773; 
                        solPntsLocalCoord2D[17][1] = 0.6647637544849;
                        solPntsLocalCoord2D[18][0] = 0.0186886898773; 
                        solPntsLocalCoord2D[18][1] = 0.3165475556378;
                        solPntsLocalCoord2D[19][0] = 0.3165475556378;
                        solPntsLocalCoord2D[19][1] = 0.6647637544849;
                        solPntsLocalCoord2D[20][0] = 0.6647637544849; 
                        solPntsLocalCoord2D[20][1] = 0.3165475556378;
                        solPntsLocalCoord2D[21][0] = 0.0192662192492; 
                        solPntsLocalCoord2D[21][1] = 0.4903668903754;
                        solPntsLocalCoord2D[22][0] = 0.4903668903754; 
                        solPntsLocalCoord2D[22][1] = 0.0192662192492;
                        solPntsLocalCoord2D[23][0] = 0.4903668903754;
                        solPntsLocalCoord2D[23][1] = 0.4903668903754;
                        solPntsLocalCoord2D[24][0] = 0.0875134669581;
                        solPntsLocalCoord2D[24][1] = 0.8249730660837;
                        solPntsLocalCoord2D[25][0] = 0.0875134669581; 
                        solPntsLocalCoord2D[25][1] = 0.0875134669581;
                        solPntsLocalCoord2D[26][0] = 0.8249730660837; 
                        solPntsLocalCoord2D[26][1] = 0.0875134669581;
                        solPntsLocalCoord2D[27][0] = 0.0935526036219; 
                        solPntsLocalCoord2D[27][1] = 0.2079865423167;
                        solPntsLocalCoord2D[28][0] = 0.0935526036219; 
                        solPntsLocalCoord2D[28][1] = 0.6984608540613;
                        solPntsLocalCoord2D[29][0] = 0.2079865423167; 
                        solPntsLocalCoord2D[29][1] = 0.0935526036219;
                        solPntsLocalCoord2D[30][0] = 0.6984608540613;
                        solPntsLocalCoord2D[30][1] = 0.0935526036219;
                        solPntsLocalCoord2D[31][0] = 0.6984608540613;
                        solPntsLocalCoord2D[31][1] = 0.2079865423167;
                        solPntsLocalCoord2D[32][0] = 0.2079865423167;
                        solPntsLocalCoord2D[32][1] = 0.6984608540613;
                        solPntsLocalCoord2D[33][0] = 0.0974892983467; 
                        solPntsLocalCoord2D[33][1] = 0.5380088595149;
                        solPntsLocalCoord2D[34][0] = 0.3645018421383; 
                        solPntsLocalCoord2D[34][1] = 0.0974892983467;
                        solPntsLocalCoord2D[35][0] = 0.5380088595149;
                        solPntsLocalCoord2D[35][1] = 0.0974892983467;
                        solPntsLocalCoord2D[36][0] = 0.5380088595149; 
                        solPntsLocalCoord2D[36][1] = 0.3645018421383; 
                        solPntsLocalCoord2D[37][0] = 0.3645018421383; 
                        solPntsLocalCoord2D[37][1] = 0.5380088595149; 
                        solPntsLocalCoord2D[38][0] = 0.0974892983467;
                        solPntsLocalCoord2D[38][1] = 0.3645018421383;
                        solPntsLocalCoord2D[39][0] = 0.2217145894873;
                        solPntsLocalCoord2D[39][1] = 0.5565708210253;
                        solPntsLocalCoord2D[40][0] = 0.5565708210253;
                        solPntsLocalCoord2D[40][1] = 0.2217145894873; 
                        solPntsLocalCoord2D[41][0] = 0.2217145894873;
                        solPntsLocalCoord2D[41][1] = 0.2217145894873;
                        solPntsLocalCoord2D[42][0] = 0.3860471669296; 
                        solPntsLocalCoord2D[42][1] = 0.2279056661408; 
                        solPntsLocalCoord2D[43][0] = 0.2279056661408; 
                        solPntsLocalCoord2D[43][1] = 0.3860471669296;
                        solPntsLocalCoord2D[44][0] = 0.3860471669296;
                        solPntsLocalCoord2D[44][1] = 0.3860471669296;                      
      } break;
          case CFPolyOrder::ORDER9:
      {
                        solPntsLocalCoord2D[0][0] = 1.0000000000000;
                        solPntsLocalCoord2D[0][1] = 0.0000000000000; 
                        solPntsLocalCoord2D[1][0] = 0.0000000000000;
                        solPntsLocalCoord2D[1][1] = 1.0000000000000; 
                        solPntsLocalCoord2D[2][0] = 0.0000000000000;
                        solPntsLocalCoord2D[2][1] = 0.0000000000000;
                        solPntsLocalCoord2D[3][0] = 0.9398863583577;
                        solPntsLocalCoord2D[3][1] = 0.0049848744634;
                        solPntsLocalCoord2D[4][0] =0.0543806683058;
                        solPntsLocalCoord2D[4][1] = 0.9386405618617; 
                        solPntsLocalCoord2D[5][0] = 0.0093940049164;
                        solPntsLocalCoord2D[5][1] = 0.0526424462697;
                        solPntsLocalCoord2D[6][0] = 0.0164345086362; 
                        solPntsLocalCoord2D[6][1] = 0.9469035517351; 
                        solPntsLocalCoord2D[7][0] = 0.9469487269862;
                        solPntsLocalCoord2D[7][1] = 0.0363373677167;
                        solPntsLocalCoord2D[8][0] = 0.0426604005768;
                        solPntsLocalCoord2D[8][1] = 0.0151224541799;
                        solPntsLocalCoord2D[9][0] = 0.0122269495439;
                        solPntsLocalCoord2D[9][1] = 0.8693773510664;
                        solPntsLocalCoord2D[10][0] = 0.8673696521047; 
                        solPntsLocalCoord2D[10][1] = 0.1204917285774;
                        solPntsLocalCoord2D[11][0] = 0.8456744021389;
                        solPntsLocalCoord2D[11][1] = 0.0157763967870;
                        solPntsLocalCoord2D[12][0] = 0.1395759632103;
                        solPntsLocalCoord2D[12][1] = 0.8448120870375;
                        solPntsLocalCoord2D[13][0] = 0.1317821743231;
                        solPntsLocalCoord2D[13][1] = 0.0135009605584;
                        solPntsLocalCoord2D[14][0] = 0.0157955126300;
                        solPntsLocalCoord2D[14][1] = 0.1455274938536;
                        solPntsLocalCoord2D[15][0] = 0.7365462884436;
                        solPntsLocalCoord2D[15][1] = 0.0155697540908;
                        solPntsLocalCoord2D[16][0] = 0.0139688430330;
                        solPntsLocalCoord2D[16][1] = 0.7379836894450;
                        solPntsLocalCoord2D[17][0] = 0.2547895186039; 
                        solPntsLocalCoord2D[17][1] = 0.7297615689771;
                        solPntsLocalCoord2D[18][0] = 0.7316386522555;
                        solPntsLocalCoord2D[18][1] = 0.2543076683315;
                        solPntsLocalCoord2D[19][0] = 0.0157253728951;
                        solPntsLocalCoord2D[19][1] = 0.2696239795791;
                        solPntsLocalCoord2D[20][0] = 0.2662302843647; 
                        solPntsLocalCoord2D[20][1] = 0.0144783956308;
                        solPntsLocalCoord2D[21][0] = 0.8673504065214; 
                        solPntsLocalCoord2D[21][1] = 0.0591679410400;
                        solPntsLocalCoord2D[22][0] = 0.0741493666957; 
                        solPntsLocalCoord2D[22][1] = 0.8634782575061;
                        solPntsLocalCoord2D[23][0] = 0.0159285948360;
                        solPntsLocalCoord2D[23][1] = 0.4191238955238;
                        solPntsLocalCoord2D[24][0] = 0.0156061028068;
                        solPntsLocalCoord2D[24][1] = 0.5809222921146;
                        solPntsLocalCoord2D[25][0] = 0.5910094817484;
                        solPntsLocalCoord2D[25][1] = 0.0159251452651;
                        solPntsLocalCoord2D[26][0] = 0.4034771496889; 
                        solPntsLocalCoord2D[26][1] = 0.5806700368104;
                        solPntsLocalCoord2D[27][0] = 0.5694745628526; 
                        solPntsLocalCoord2D[27][1] = 0.4149495146302;
                        solPntsLocalCoord2D[28][0] = 0.0678493700650;
                        solPntsLocalCoord2D[28][1] = 0.0761218678591;
                        solPntsLocalCoord2D[29][0] = 0.4265968590272; 
                        solPntsLocalCoord2D[29][1] = 0.0157509692312;
                        solPntsLocalCoord2D[30][0] = 0.0670982507890;
                        solPntsLocalCoord2D[30][1] = 0.7741898312421;
                        solPntsLocalCoord2D[31][0] = 0.7528310231480;
                        solPntsLocalCoord2D[31][1] = 0.0819119495639;
                        solPntsLocalCoord2D[32][0] = 0.7753727783557;
                        solPntsLocalCoord2D[32][1] = 0.1577128457292;
                        solPntsLocalCoord2D[33][0] = 0.1689073157787; 
                        solPntsLocalCoord2D[33][1] = 0.7503943099742;
                        solPntsLocalCoord2D[34][0] = 0.1687335832919; 
                        solPntsLocalCoord2D[34][1] = 0.0708311507268;
                        solPntsLocalCoord2D[35][0] = 0.0821244708436;
                        solPntsLocalCoord2D[35][1] = 0.1762996626771;
                        solPntsLocalCoord2D[36][0] = 0.6288705363345; 
                        solPntsLocalCoord2D[36][1] = 0.0807744953317; 
                        solPntsLocalCoord2D[37][0] = 0.0811413015266; 
                        solPntsLocalCoord2D[37][1] = 0.3054373589776; 
                        solPntsLocalCoord2D[38][0] = 0.2969112065080;
                        solPntsLocalCoord2D[38][1] = 0.6227485988871;
                        solPntsLocalCoord2D[39][0] = 0.0767542314171;
                        solPntsLocalCoord2D[39][1] = 0.6247247149546;
                        solPntsLocalCoord2D[40][0] = 0.6223022333845;
                        solPntsLocalCoord2D[40][1] = 0.3011485821166; 
                        solPntsLocalCoord2D[41][0] = 0.3103786288051;
                        solPntsLocalCoord2D[41][1] = 0.0779098365079;
                        solPntsLocalCoord2D[42][0] = 0.0819218215187; 
                        solPntsLocalCoord2D[42][1] = 0.4603633038351; 
                        solPntsLocalCoord2D[43][0] = 0.4717022665013; 
                        solPntsLocalCoord2D[43][1] = 0.0821554006797;
                        solPntsLocalCoord2D[44][0] = 0.4546603415250;
                        solPntsLocalCoord2D[44][1] = 0.4637565033890;
                        solPntsLocalCoord2D[45][0] = 0.1701091339237;
                        solPntsLocalCoord2D[45][1] = 0.6422277808188; 
                        solPntsLocalCoord2D[46][0] = 0.6406004329487;
                        solPntsLocalCoord2D[46][1] = 0.1898293537256; 
                        solPntsLocalCoord2D[47][0] = 0.1912267583717;
                        solPntsLocalCoord2D[47][1] = 0.1739955685343;
                        solPntsLocalCoord2D[48][0] = 0.1885315767070;                        
                        solPntsLocalCoord2D[48][1] = 0.4798914070406;
                        solPntsLocalCoord2D[49][0] = 0.4772929957691;
                        solPntsLocalCoord2D[49][1] = 0.3348356598119; 
                        solPntsLocalCoord2D[50][0] = 0.3126974621760;
                        solPntsLocalCoord2D[50][1] = 0.4957972197259;
                        solPntsLocalCoord2D[51][0] = 0.4961225945946; 
                        solPntsLocalCoord2D[51][1] = 0.1927553668904; 
                        solPntsLocalCoord2D[52][0] = 0.1928805312867;
                        solPntsLocalCoord2D[52][1] = 0.3161015807261;
                        solPntsLocalCoord2D[53][0] = 0.3360041453816;
                        solPntsLocalCoord2D[53][1] = 0.1894892801290;
                        solPntsLocalCoord2D[54][0] = 0.3337280550848;
                        solPntsLocalCoord2D[54][1] = 0.3343571021811;
                      
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
    //faces ZTA=1 and -1 --> triag faces
    flxCoords[ZTA] = 1;
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts2DTriag; ++iFlx)
    {
      flxCoords[KSI] = solPntsLocalCoord2D[iFlx][KSI];
      flxCoords[ETA] = solPntsLocalCoord2D[iFlx][ETA];
      m_flxPntsLocalCoords.push_back(flxCoords);
    }
    flxCoords[ZTA] = -1;
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts2DTriag; ++iFlx)
    {
      flxCoords[KSI] = solPntsLocalCoord2D[iFlx][KSI];
      flxCoords[ETA] = solPntsLocalCoord2D[iFlx][ETA];
      m_flxPntsLocalCoords.push_back(flxCoords);
    }
    // 3 rectangular faces 

    //face at eta=0
    flxCoords[ETA] = 0;
    for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
    {
        for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
        {
            flxCoords[KSI] = (m_flxPntsLocalCoord1D[iKsi]+1.)/2.;
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }
    // (Oblique) face at ksi+eta=1
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
        for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
        {
            flxCoords[KSI] = 1-alpha[iEta];
            flxCoords[ETA] = 0+alpha[iEta];
            flxCoords[ZTA] = m_flxPntsLocalCoord1D[iZta];
            m_flxPntsLocalCoords.push_back(flxCoords);
        }
    }

    //face at ksi=0
    flxCoords[KSI] = 0;
    for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
    {
        for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
        {
            flxCoords[ETA] = (m_flxPntsLocalCoord1D[iEta]+1.)/2.;
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

void PrismFluxReconstructionElementData::createSolPolyExponents()
{
  CFAUTOTRACE;

  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  // define exponents
  m_solPolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
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
    for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
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
  CFuint iFlxg = 0; 
  // zeroth face (triag face zta=1)
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx, ++iFlxg)
  {
    m_faceFlxPntConn[faceIdx].push_back(iFlxg); 
  }
  ++faceIdx;

  // first face (triag face zta=-1)
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx, ++iFlxg)
  {
    m_faceFlxPntConn[faceIdx].push_back(iFlxg); 
  }
  ++faceIdx;

  // second face (eta=0)
  for (CFuint iKsi = 0; iKsi < nbrFlxPnts1D; ++iKsi)
  {
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta, ++iFlxg)
    {
      m_faceFlxPntConn[faceIdx].push_back(iFlxg);
    }
  }
  ++faceIdx;

  // third face (ksi+eta=1)
  for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
  {
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta)
    {
      m_faceFlxPntConn[faceIdx].push_back(iFlxg);
    }
  }
  ++faceIdx;
  
  // fourth face (ksi=0)
  for (CFuint iEta = 0; iEta < nbrFlxPnts1D; ++iEta)
  {
    for (CFuint iZta = 0; iZta < nbrFlxPnts1D; ++iZta, ++iFlxg)
    {
      m_faceFlxPntConn[faceIdx].push_back(iFlxg);
    }
  }

}

//////////////////////////////////////////////////////////////////////
// !!!! to be checked

void PrismFluxReconstructionElementData::createFaceFluxPntsConnPerOrient()
{
  CFAUTOTRACE;
  // number of faces
  const CFuint nbrFaces = m_faceNodeConn.size();
  CFLog(VERBOSE,"nbrFaces:" << nbrFaces << "\n");
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
  CFLog(VERBOSE,"nbrSolPnts1D:" << nbrSolPnts1D << "\n");
  // number of face flux points
  const CFuint nbrFaceFlxPnts = getNbrOfFaceFlxPnts();
  CFLog(VERBOSE,"nbrFaceFlxPnts:" << nbrFaceFlxPnts << "\n");
  // flux point indexes for inverted face
  vector< CFuint > invFlxIdxs;
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    const CFuint idxKsi = nbrSolPnts1D - iKsi - 1;
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
    {
      invFlxIdxs.push_back(nbrSolPnts1D*idxKsi+iEta);
    }
  }
  // number of rotatable flux point groups
  const CFuint nbrRotFlxGroups = nbrFaceFlxPnts/4;
  CFLog(VERBOSE,"nbrRotFlxGroups:" << nbrRotFlxGroups << "\n");
  // maximum number of flux points in a line of flux points
  const CFuint maxNbrFlxPntsInLine = (nbrSolPnts1D+1)/2;
  CFLog(VERBOSE,"maxNbrFlxPntsInLine:" << maxNbrFlxPntsInLine << "\n");
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

  // number of possible orientations
  const CFuint nbrOrient = 33;

  // create data structure
  m_faceFlxPntConnPerOrient.resize(nbrOrient);
  CFuint iOrient = 0;
  for (CFuint iFaceL = 2; iFaceL < nbrFaces; ++iFaceL)
  {
    vector< CFuint > faceFlxConnL = m_faceFlxPntConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < nbrFaces; ++iFaceR)
    {
      vector< CFuint > faceFlxConnR = m_faceFlxPntConn[iFaceR];
      for (CFuint iRot = 0; iRot < 4; ++iRot, ++iOrient)
      {
        m_faceFlxPntConnPerOrient[iOrient].resize(2);
        for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
        {
          m_faceFlxPntConnPerOrient[iOrient][LEFT ]
              .push_back(faceFlxConnL[iFlx            ]);
          m_faceFlxPntConnPerOrient[iOrient][RIGHT]
              .push_back(faceFlxConnR[invFlxIdxs[iFlx]]);
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
          const CFuint swap = invFlxIdxs[flx3Idx];
          invFlxIdxs[flx3Idx] = invFlxIdxs[flx2Idx];
          invFlxIdxs[flx2Idx] = invFlxIdxs[flx1Idx];
          invFlxIdxs[flx1Idx] = invFlxIdxs[flx0Idx];
          invFlxIdxs[flx0Idx] = swap;
        }
      }
    }
  }

// flux point indexes for inverted face
    CFuint nbrSolPntsTriag=  (m_polyOrder+1)*(m_polyOrder+2)/2;

    for (CFuint iEta = 0; iEta < nbrSolPntsTriag; ++iEta)
    {
      invFlxIdxs.push_back(nbrSolPntsTriag-1-iEta);
    }
  

    for (CFuint iFaceL = 0; iFaceL < 2; ++iFaceL)
  {
    vector< CFuint > faceFlxConnL = m_faceFlxPntConn[iFaceL];
    for (CFuint iFaceR = iFaceL; iFaceR < 2; ++iFaceR)
    {
      vector< CFuint > faceFlxConnR = m_faceFlxPntConn[iFaceR];
      for (CFuint iRot = 0; iRot < 3; ++iRot, ++iOrient)
      {
        m_faceFlxPntConnPerOrient[iOrient].resize(2);
        for (CFuint iFlx = 0; iFlx < nbrSolPntsTriag ; ++iFlx)
        {
          m_faceFlxPntConnPerOrient[iOrient][LEFT ]
              .push_back(faceFlxConnL[iFlx            ]);
          m_faceFlxPntConnPerOrient[iOrient][RIGHT]
              .push_back(faceFlxConnR[invFlxIdxs[iFlx]]);
        }
        // rotate the right face
        for (CFuint iRotGroup = 0; iRotGroup < nbrRotFlxGroups; ++iRotGroup)
        {
          // indexes of flux points to be rotated
          const CFuint flx0Idx = rotFlxIdxs[iRotGroup][0];
          const CFuint flx1Idx = rotFlxIdxs[iRotGroup][1];
          const CFuint flx2Idx = rotFlxIdxs[iRotGroup][2];
          //CFLog(VERBOSE,"flxrotIdx" << flx0Idx << flx1Idx << flx2Idx << flx3Idx << "\n");

          // rotate flux points
          const CFuint swap = invFlxIdxs[flx2Idx];
          invFlxIdxs[flx2Idx] = invFlxIdxs[flx1Idx];
          invFlxIdxs[flx1Idx] = invFlxIdxs[flx0Idx];
          invFlxIdxs[flx0Idx] = swap;
        }
      }
    }
  }
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
  m_cellNodeCoords[0][ZTA] = +1.0;

  // second node
  m_cellNodeCoords[1].resize(3);
  m_cellNodeCoords[1][KSI] = +1.0;
  m_cellNodeCoords[1][ETA] = +0.0;
  m_cellNodeCoords[1][ZTA] = +1.0;

  // third node
  m_cellNodeCoords[2].resize(3);
  m_cellNodeCoords[2][KSI] = +0.0;
  m_cellNodeCoords[2][ETA] = +1.0;
  m_cellNodeCoords[2][ZTA] = +1.0;

  // fourth node
  m_cellNodeCoords[3].resize(3);
  m_cellNodeCoords[3][KSI] = +0.0;
  m_cellNodeCoords[3][ETA] = +0.0;
  m_cellNodeCoords[3][ZTA] = -1.0;

  // fifth node
  m_cellNodeCoords[4].resize(3);
  m_cellNodeCoords[4][KSI] = +1.0;
  m_cellNodeCoords[4][ETA] = +0.0;
  m_cellNodeCoords[4][ZTA] = -1.0;

  // sixth node
  m_cellNodeCoords[5].resize(3);
  m_cellNodeCoords[5][KSI] = +0.0;
  m_cellNodeCoords[5][ETA] = +1.0;
  m_cellNodeCoords[5][ZTA] = -1.0;

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

  m_faceMappedCoordDir.resize(6);

  m_faceMappedCoordDir[0] = 1;
  m_faceMappedCoordDir[1] = 1;
  m_faceMappedCoordDir[2] = 1;
  m_faceMappedCoordDir[3] = 1;
  m_faceMappedCoordDir[4] = 1;
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
      else
      {
	      m_flxPntFlxDim[m_faceFlxPntConn[iFace][iFlx]] = 4;
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
  m_faceNormals[0][2] = 1.;

  m_faceNormals[1][0] = 0.;
  m_faceNormals[1][1] = 0.;
  m_faceNormals[1][2] = -1.;

  m_faceNormals[2][0] = 0.;
  m_faceNormals[2][1] = -1.;
  m_faceNormals[2][2] = 0.;

  m_faceNormals[3][0] = 0.5;
  m_faceNormals[3][1] = 0.5;
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
            for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_faceNodeConnPerOrient[iOrient][iSide].resize(3);
    }
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
  
cout << "createFaceNodeConnectivityPerOrient -- Nb Orientation=  "<<iOrient<<endl;
  //cf_assert(iOrient == nbrOrient);

/*////////////////////////
  // resize the variable
  m_faceNodeConnPerOrient.resize(nbrOrient);
  m_faceConnPerOrient.resize(nbrOrient);
  m_faceMappedCoordDirPerOrient.resize(nbrOrient);
  for (CFuint iOrient = 0; iOrient < nbrOrient; ++iOrient)
  {
    m_faceNodeConnPerOrient[iOrient].resize(2);
    m_faceConnPerOrient[iOrient].resize(2);
    m_faceMappedCoordDirPerOrient[iOrient].resize(2);
    for (CFuint iCell = 0; iCell < 2; ++iCell)
      m_faceNodeConnPerOrient[iOrient][iCell].resize(2);
  }

  // fill the variable
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < nbrFaces; ++iFaceL)
  {
    for (CFuint iFaceR = iFaceL; iFaceR < nbrFaces; ++iFaceR, ++iOrient)
    {
      m_faceConnPerOrient[iOrient][LEFT ] = iFaceL;
      m_faceConnPerOrient[iOrient][RIGHT] = iFaceR;

      m_faceMappedCoordDirPerOrient[iOrient][LEFT ] = 2.;
      m_faceMappedCoordDirPerOrient[iOrient][RIGHT] = -2.;

      for (CFuint iNode = 0; iNode < 2; ++iNode)
      {
        m_faceNodeConnPerOrient[iOrient][LEFT ][iNode] = m_faceNodeConn[iFaceL][iNode  ];
        m_faceNodeConnPerOrient[iOrient][RIGHT][iNode] = m_faceNodeConn[iFaceR][1-iNode];
      }
    }
  }
////////////////////////*/

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

  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tpIntegrator(DIM_3D,m_polyOrder);

  // create cell node local coordinates
  vector< RealVector > nodeCoord(8);
  nodeCoord[0].resize(3);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[0][ETA] = -1.0;
  nodeCoord[0][ZTA] = -1.0;
  nodeCoord[1].resize(3);
  nodeCoord[1][KSI] = +1.0;
  nodeCoord[1][ETA] = -1.0;
  nodeCoord[1][ZTA] = -1.0;
  nodeCoord[2].resize(3);
  nodeCoord[2][KSI] = +1.0;
  nodeCoord[2][ETA] = +1.0;
  nodeCoord[2][ZTA] = -1.0;
  nodeCoord[3].resize(3);
  nodeCoord[3][KSI] = -1.0;
  nodeCoord[3][ETA] = +1.0;
  nodeCoord[3][ZTA] = -1.0;
  nodeCoord[4].resize(3);
  nodeCoord[4][KSI] = -1.0;
  nodeCoord[4][ETA] = -1.0;
  nodeCoord[4][ZTA] = +1.0;
  nodeCoord[5].resize(3);
  nodeCoord[5][KSI] = +1.0;
  nodeCoord[5][ETA] = -1.0;
  nodeCoord[5][ZTA] = +1.0;
  nodeCoord[6].resize(3);
  nodeCoord[6][KSI] = +1.0;
  nodeCoord[6][ETA] = +1.0;
  nodeCoord[6][ZTA] = +1.0;
  nodeCoord[7].resize(3);
  nodeCoord[7][KSI] = -1.0;
  nodeCoord[7][ETA] = +1.0;
  nodeCoord[7][ZTA] = +1.0;

  // get quadrature point coordinates and wheights
  vector< RealVector > quadPntCoords   = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  vector< CFreal     > quadPntWheights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWheights.size() == nbrQPnts);

  // get the solution polynomial values at the quadrature points
  vector< vector< CFreal > > quadPntPolyVals = getSolPolyValsAtNode(quadPntCoords);

  // compute the coefficients for integration over a face
  // loop over solution points
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_cellAvgSolCoefs[iSol] = 0.0;
    for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
    {
      m_cellAvgSolCoefs[iSol] += quadPntWheights[iQPnt]*quadPntPolyVals[iQPnt][iSol];
    }
    m_cellAvgSolCoefs[iSol] *= 0.125;
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
// !!!! to be checked

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
  
  if(m_polyOrder != CFPolyOrder::ORDER0 && m_polyOrder != CFPolyOrder::ORDER1)
  {
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      CFuint modalDof = 0;
      RealVector orderIdx(nbrSolPnts1D);
      orderIdx = 0.0;
      
      for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
        {
          for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++modalDof)
          {
            order = max(iKsi,iEta);
            order = max(order,iZta);
            idx = order*order*order+orderIdx[order];
            orderIdx[order] += 1;
              
            m_vandermonde(iSol,idx) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iEta)*evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA],iZta);
          }
        }
      }
//       for (CFuint iOrder = 0; iOrder < nbrSolPnts1D; ++iOrder)
//       {
//         for (CFuint iOrderKsi = 0; iOrderKsi < iOrder; ++iOrderKsi, ++modalDof)
//         {
// 	  CFuint iOrderEta = iOrder;
// 	  CFuint iOrderZta = iOrder;
//           m_vandermonde(iSol,modalDof) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iOrderEta)*evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA],iOrderZta);
//         }
//         for (CFuint iOrderEta = 0; iOrderEta < iOrder+1; ++iOrderEta, ++modalDof)
//         {
// 	  CFuint iOrderKsi = iOrder;
// 	  CFuint iOrderZta = iOrder;
//           m_vandermonde(iSol,modalDof) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iOrderEta)*evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA],iOrderZta);
//         }
//         for (CFuint iOrderZta = 0; iOrderZta < iOrder+2; ++iOrderZta, ++modalDof)
//         {
// 	  CFuint iOrderKsi = iOrder;
// 	  CFuint iOrderEta = iOrder;
//           m_vandermonde(iSol,modalDof) = evaluateLegendre(m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(m_solPntsLocalCoords[iSol][ETA],iOrderEta)*evaluateLegendre(m_solPntsLocalCoords[iSol][ZTA],iOrderZta);
//         }
//       }
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
  const CFuint nbrFlxPnts = m_flxPntsLocalCoords.size();

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

 /* const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  const CFuint nbrFlxPnts = m_flxPntsLocalCoords.size();
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();

  m_solFlxDep.resize(nbrSolPnts);
  m_solSolDep.resize(nbrSolPnts);
  m_flxSolDep.resize(nbrFlxPnts);
  m_closestSolToFlxIdx.resize(nbrFlxPnts);

  CFuint iSol = 0;
  
  for (CFuint iKsi = 0; iKsi < nbrSolPnts1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbrSolPnts1D; ++iEta)
    {
      for (CFuint iZta = 0; iZta < nbrSolPnts1D; ++iZta, ++iSol)
      {
        m_solFlxDep[iSol].push_back(iEta + iKsi*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iEta + iKsi*nbrSolPnts1D + nbrSolPnts1D*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iZta + iEta*nbrSolPnts1D + 2*nbrSolPnts1D*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iZta + iEta*nbrSolPnts1D + 3*nbrSolPnts1D*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iKsi + iZta*nbrSolPnts1D + 4*nbrSolPnts1D*nbrSolPnts1D);
        m_solFlxDep[iSol].push_back(iKsi + iZta*nbrSolPnts1D + 5*nbrSolPnts1D*nbrSolPnts1D);        

        m_flxSolDep[iEta + iKsi*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iEta + iKsi*nbrSolPnts1D + nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iZta + iEta*nbrSolPnts1D + 2*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iZta + iEta*nbrSolPnts1D + 3*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iKsi + iZta*nbrSolPnts1D + 4*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);
        m_flxSolDep[iKsi + iZta*nbrSolPnts1D + 5*nbrSolPnts1D*nbrSolPnts1D].push_back(iSol);

        for (CFuint jSol = 0; jSol < nbrSolPnts1D; ++jSol)
        {
          m_solSolDep[iSol].push_back(iKsi*nbrSolPnts1D*nbrSolPnts1D + iEta*nbrSolPnts1D + jSol);
          if (iSol != iZta + iKsi*nbrSolPnts1D*nbrSolPnts1D + jSol*nbrSolPnts1D) m_solSolDep[iSol].push_back(iZta + iKsi*nbrSolPnts1D*nbrSolPnts1D + jSol*nbrSolPnts1D);
          if (iSol != iZta + jSol*nbrSolPnts1D*nbrSolPnts1D + iEta*nbrSolPnts1D) m_solSolDep[iSol].push_back(iZta + jSol*nbrSolPnts1D*nbrSolPnts1D + iEta*nbrSolPnts1D);
        }
      }
    }
  }

  for (CFuint i = 0; i < nbrSolPnts1D; ++i)
  {
    for (CFuint j = 0; j < nbrSolPnts1D; ++j)
    {
      // face 0 and 1
      m_closestSolToFlxIdx[j+i*nbrSolPnts1D] = j*nbrSolPnts1D+i*nbrSolPnts1D*nbrSolPnts1D; 
      m_closestSolToFlxIdx[j+i*nbrSolPnts1D + nbrSolPnts1D*nbrSolPnts1D] = j*nbrSolPnts1D+i*nbrSolPnts1D*nbrSolPnts1D + nbrSolPnts1D-1; 
      
      // face 2 and 3
      m_closestSolToFlxIdx[j+i*nbrSolPnts1D + 2*nbrSolPnts1D*nbrSolPnts1D] = j+i*nbrSolPnts1D; 
      m_closestSolToFlxIdx[j+i*nbrSolPnts1D + 3*nbrSolPnts1D*nbrSolPnts1D] = j+i*nbrSolPnts1D + nbrSolPnts1D*nbrSolPnts1D*(nbrSolPnts1D-1); 
        
      // face 4 and 5
      m_closestSolToFlxIdx[j+i*nbrSolPnts1D + 4*nbrSolPnts1D*nbrSolPnts1D] = i+j*nbrSolPnts1D*nbrSolPnts1D; 
      m_closestSolToFlxIdx[j+i*nbrSolPnts1D + 5*nbrSolPnts1D*nbrSolPnts1D] = i+j*nbrSolPnts1D*nbrSolPnts1D + nbrSolPnts1D*(nbrSolPnts1D-1); 
    }
  }*/
}

//////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
