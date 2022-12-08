#include <fstream>
#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh"
#include "MathTools/MathConsts.hh"
#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;


namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////

TriagFluxReconstructionElementData::TriagFluxReconstructionElementData() :
  FluxReconstructionElementData()
{
  m_shape = CFGeoShape::TRIAG;
  m_dimensionality = DIM_2D;
}

//////////////////////////////////////////////////////////////////////

TriagFluxReconstructionElementData::TriagFluxReconstructionElementData(CFPolyOrder::Type polyOrder)
{
  m_shape = CFGeoShape::TRIAG;
  m_dimensionality = DIM_2D;
  m_polyOrder = polyOrder;
  m_solPntsLocalCoord1D.resize(polyOrder+1);
  m_flxPntsLocalCoord1D.resize(polyOrder+1);

  CFuint rows = (0.5*(polyOrder+1)*(polyOrder+2)) ;
  std::vector< std::vector< CFreal > > coords( 2 , std::vector<CFreal> (rows)); 
  coords.resize(rows, std::vector< CFreal > (2));
  solPntsLocalCoord2D.resize(rows, std::vector< CFreal > (2));
  CFuint nbfaces = 3;
  flxPntsLocalCoord2D.resize(nbfaces, std::vector< std::vector<CFreal> > (polyOrder+1));
	for (CFuint iFace=0; iFace <nbfaces; ++iFace){
    for (CFuint iFlx=0; iFlx<polyOrder+1; ++iFlx){
		  (flxPntsLocalCoord2D[iFace][iFlx]).resize(0);
	  }
  }
  
  for (CFuint iFace=0; iFace <nbfaces; ++iFace){
    for (CFuint iFlx=0; iFlx<polyOrder+1; ++iFlx){
		  flxPntsLocalCoord2D[iFace][iFlx].push_back(0.);
		  flxPntsLocalCoord2D[iFace][iFlx].push_back(0.);
	  }
  }

  switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
                        coords[0][0] = 0.333333333333333; 
                        coords[0][1] = 0.333333333333333;  
      } break;
      case CFPolyOrder::ORDER1:
      {
                        coords[0][0] = 0.1666666666667; 
                        coords[0][1] = 0.1666666666667;
                        coords[1][0] = 0.6666666666667; 
                        coords[1][1] = 0.1666666666667;  
                        coords[2][0] = 0.1666666666667; 
                        coords[2][1] = 0.6666666666667;

    } break;

      case CFPolyOrder::ORDER2:              
      {
                        coords[0][0] =  0.091576213509780; 
                        coords[0][1] =  0.091576213509780; 
                        
                        coords[1][0] =  0.816847572980440; 
                        coords[1][1] =  0.091576213509780;  
                        
                        coords[2][0] =  0.091576213509780; 
                        coords[2][1] =  0.816847572980440;
                        
                        coords[3][0] =  0.445948490915964; 
                        coords[3][1] =  0.108103018168071; 
                        
                        coords[4][0] =  0.108103018168071; 
                        coords[4][1] =  0.445948490915964;
                        
                        coords[5][0] =  0.445948490915964; 
                        coords[5][1] =  0.445948490915964;


      } break;

      case CFPolyOrder::ORDER3:          
      {
                        coords[0][0] = 0.055564052669793;
                        coords[0][1] = 0.055564052669793;
                        
                        coords[1][0] = 0.888871894660413; 
                        coords[1][1] = 0.055564052669793;	
                        
                        coords[2][0] = 0.055564052669793; 
                        coords[2][1] = 0.888871894660413;  	

                        
                        coords[3][0] = 0.634210747745723; 
                        coords[3][1] = 0.070255540518384;  	

                        
                        coords[4][0] = 0.070255540518384; 
                        coords[4][1] = 0.634210747745723;	
                        
                        coords[5][0] = 0.295533711735893; 
                        coords[5][1] = 0.634210747745723; 	

                        
                        coords[6][0] = 0.295533711735893; 
                        coords[6][1] = 0.070255540518384;  	
                        
                        coords[7][0] = 0.070255540518384; 
                        coords[7][1] = 0.295533711735893; 	
                        
                        coords[8][0] = 0.634210747745723; 
                        coords[8][1] = 0.295533711735893; 	
                        
                        coords[9][0] = 0.333333333333333; 
                        coords[9][1] = 0.333333333333333;


      } break;
      case CFPolyOrder::ORDER4:   
      {
                        coords[0][0] = 0.035870877695734; 
                        coords[0][1] = 0.035870877695734;

                        coords[1][0] = 0.928258244608533; 
                        coords[1][1] = 0.035870877695734;

                        coords[2][0] = 0.035870877695734; 
                        coords[2][1] = 0.928258244608533; 

                        coords[3][0] = 0.241729395767967;
                        coords[3][1] = 0.241729395767967;

                        coords[4][0] = 0.516541208464066; 
                        coords[4][1] = 0.241729395767967; 

                        coords[5][0] = 0.241729395767967; 
                        coords[5][1] = 0.516541208464066; 

                        coords[6][0] = 0.474308787777079;
                        coords[6][1] = 0.051382424445843; 

                        coords[7][0] = 0.051382424445843;
                        coords[7][1] = 0.474308787777079;

                        coords[8][0] = 0.474308787777079; 
                        coords[8][1] = 0.474308787777079;

                        coords[9][0] = 0.751183631106484; 
                        coords[9][1] = 0.047312487011716;

                        coords[10][0] = 0.047312487011716; 
                        coords[10][1] = 0.751183631106484;

                        coords[11][0] =  0.201503881881800;
                        coords[11][1] =  0.751183631106484;

                        coords[12][0] = 0.201503881881800; 
                        coords[12][1] =  0.047312487011716;

                        coords[13][0] = 0.047312487011716;
                        coords[13][1] = 0.201503881881800;

                        coords[14][0] = 0.751183631106484; 
                        coords[14][1] = 0.201503881881800;

      } break;

      case CFPolyOrder::ORDER5:
      {
                        coords[0][0] = 0.028112952182664; 
                        coords[0][1] = 0.028112952182664; 

                        coords[1][0] = 0.943774095634672;
                        coords[1][1] = 0.028112952182664; 

                        coords[2][0] = 0.028112952182664;
                        coords[2][1] = 0.943774095634672;

                        coords[3][0] = 0.177139098469317;
                        coords[3][1] = 0.177139098469317;

                        coords[4][0] = 0.645721803061365; 
                        coords[4][1] = 0.177139098469317;

                        coords[5][0] = 0.177139098469317;
                        coords[5][1] = 0.645721803061365;

                        coords[6][0] = 0.405508595867433;
                        coords[6][1] = 0.188982808265134;

                        coords[7][0] = 0.188982808265134;
                        coords[7][1] = 0.405508595867433;

                        coords[8][0] = 0.405508595867433;
                        coords[8][1] = 0.405508595867433;

                        coords[9][0] = 0.817900980028499; 
                        coords[9][1] = 0.033533207700614;

                        coords[10][0] = 0.033533207700614; 
                        coords[10][1] = 0.817900980028499;

                        coords[11][0] = 0.148565812270887;
                        coords[11][1] = 0.817900980028499;

                        coords[12][0] = 0.148565812270887; 
                        coords[12][1] = 0.033533207700614; 

                        coords[13][0] = 0.033533207700614;
                        coords[13][1] = 0.148565812270887;

                        coords[14][0] = 0.817900980028499;
                        coords[14][1] = 0.148565812270887;

                        coords[15][0] = 0.604978911775132; 
                        coords[15][1] = 0.037824789609186;

                        coords[16][0] = 0.037824789609186; 
                        coords[16][1] = 0.604978911775132;

                        coords[17][0] = 0.357196298615681;
                        coords[17][1] = 0.604978911775132;

                        coords[18][0] = 0.357196298615681; 
                        coords[18][1] = 0.037824789609186;

                        coords[19][0] = 0.037824789609186; 
                        coords[19][1] = 0.357196298615681; 

                        coords[20][0] = 0.604978911775132; 
                        coords[20][1] = 0.357196298615681; 

      } break;


   case CFPolyOrder::ORDER6:
      {  
                        coords[0][0] =  0.0000000000000;
                        coords[0][1] =  0.9451704450174; 
                        coords[1][0] =  0.9451704450173; 
                        coords[1][1] =  0.0000000000000; 
                        coords[2][0] =  0.9289002405719;
                        coords[2][1] =  0.0685505797224;
                        coords[3][0] =  0.0685505797224;
                        coords[3][1] =  0.9289002405717;
                        coords[4][0] =  0.0243268355615;
                        coords[4][1] =  0.0243268355616; 
                        coords[5][0] =  0.1279662835335;
                        coords[5][1] =  0.0277838749488;
                        coords[6][0] =  0.0277838749488; 
                        coords[6][1] =  0.1279662835337; 
                        coords[7][0] =  0.0287083428360; 
                        coords[7][1] =  0.7498347588657;
                        coords[8][0] =  0.7498347588656;
                        coords[8][1] =  0.0287083428360;
                        coords[9][0] =  0.7228007909707;
                        coords[9][1] =  0.2497602062385;
                        coords[10][0] = 0.2497602062386; 
                        coords[10][1] = 0.7228007909707;
                        coords[11][0] = 0.0865562992839;
                        coords[11][1] = 0.8325513856997;
                        coords[12][0] = 0.8325513856998;
                        coords[12][1] = 0.0865562992839;
                        coords[13][0] = 0.3061619157672;
                        coords[13][1] = 0.0303526617491;
                        coords[14][0] = 0.0303526617491;
                        coords[14][1] = 0.3061619157675;
                        coords[15][0] = 0.4868610595047;
                        coords[15][1] = 0.4868610595047;
                        coords[16][0] = 0.6657904293017;
                        coords[16][1] = 0.1765456154219;
                        coords[17][0] = 0.1765456154221; 
                        coords[17][1] = 0.6657904293016;
                        coords[18][0] = 0.0293121007360; 
                        coords[18][1] = 0.5295657488669;
                        coords[19][0] = 0.5295657488667;
                        coords[19][1] = 0.0293121007360;
                        coords[20][0] = 0.1444673824391; 
                        coords[20][1] = 0.1444673824391;
                        coords[21][0] = 0.3299740111411; 
                        coords[21][1] = 0.5361815729050;
                        coords[22][0] = 0.5361815729052; 
                        coords[22][1] = 0.3299740111409;
                        coords[23][0] = 0.5511507516862;
                        coords[23][1] = 0.1437790861923;
                        coords[24][0] = 0.1437790861923;
                        coords[24][1] = 0.5511507516862;
                        coords[25][0] = 0.3348066587327; 
                        coords[25][1] = 0.1529619437161;
                        coords[26][0] = 0.1529619437161; 
                        coords[26][1] = 0.3348066587327;
                        coords[27][0] = 0.3430183498147; 
                        coords[27][1] = 0.3430183498147;
                      
      } break;


      case CFPolyOrder::ORDER7:
      {
                        coords[0][0] = 0.0242935351590;
                        coords[0][1] = 0.9493059293846; 
                        coords[1][0] = 0.0265193427722;
                        coords[1][1] = 0.0242695130640; 
                        coords[2][0] = 0.9492126023551;
                        coords[2][1] = 0.0265067966437;
                        coords[3][0] = 0.0033775763749;
                        coords[3][1] = 0.4767316412363;
                        coords[4][0] = 0.4757672298101;
                        coords[4][1] = 0.5198921829102; 
                        coords[5][0] = 0.5190783193471;
                        coords[5][1] = 0.0055912706202;
                        coords[6][0] = 0.8616839745321; 
                        coords[6][1] = 0.0133996048618; 
                        coords[7][0] = 0.1249209759926;
                        coords[7][1] = 0.8613054321334;
                        coords[8][0] = 0.0138565453861;
                        coords[8][1] = 0.1247733717358;
                        coords[9][0] = 0.0211887064222;
                        coords[9][1] = 0.8438438351223;
                        coords[10][0]  =0.8432296787219; 
                        coords[10][1] = 0.1354563645830;
                        coords[11][0] = 0.1354231797865;
                        coords[11][1] = 0.0213482820656;
                        coords[12][0] = 0.3088853510679;
                        coords[12][1] = 0.0221919663014;
                        coords[13][0] = 0.6685057595169;
                        coords[13][1] = 0.3089012879389;
                        coords[14][0] = 0.0226545012557;
                        coords[14][1] = 0.6691709943321;
                        coords[15][0] = 0.2808515408772;
                        coords[15][1] = 0.6924718155106;
                        coords[16][0] = 0.6922446749051;
                        coords[16][1] = 0.0268723345026;
                        coords[17][0] = 0.0268617447119; 
                        coords[17][1] = 0.2810093973222;
                        coords[18][0] = 0.1141778485470; 
                        coords[18][1] = 0.7973581413586;
                        coords[19][0] = 0.7974807922061;
                        coords[19][1] = 0.0879806508791;
                        coords[20][0] = 0.0892807293894; 
                        coords[20][1] = 0.1145020561128;
                        coords[21][0] = 0.1052487892455; 
                        coords[21][1] = 0.6686904119922;
                        coords[22][0] = 0.6663022280740; 
                        coords[22][1] = 0.2275051631832;
                        coords[23][0] = 0.2307803737547;
                        coords[23][1] = 0.1054572561221;
                        coords[24][0] = 0.1705059157540;
                        coords[24][1] = 0.5174064398658;
                        coords[25][0] = 0.5086593973043; 
                        coords[25][1] = 0.3170523855209;
                        coords[26][0] = 0.3141823862281; 
                        coords[26][1] = 0.1810706361659;
                        coords[27][0] = 0.4617460817864; 
                        coords[27][1] = 0.4678594539804;
                        coords[28][0] = 0.0693087496081; 
                        coords[28][1] = 0.4622856042085;
                        coords[29][0] = 0.4622856042085; 
                        coords[29][1] = 0.0724357805669;
                        coords[30][0] = 0.2578625857893;
                        coords[30][1] = 0.6131395039177;
                        coords[31][0] = 0.6112627766779;
                        coords[31][1] = 0.1300360834609;
                        coords[32][0] = 0.1305182135934;
                        coords[32][1] = 0.2581713828884;
                        coords[33][0] = 0.4281437991828; 
                        coords[33][1] = 0.2362005969817;
                        coords[34][0] = 0.3356995783730; 
                        coords[34][1] = 0.4311026308588;
                        coords[35][0] = 0.2305424298836; 
                        coords[35][1] = 0.3456013949376;

                      
      } break;

        case CFPolyOrder::ORDER8:
      {
                        coords[0][0] = 0.0000000000001;
                        coords[0][1] = 0.9999999999999; 
                        coords[1][0] = 0.9999999999999;
                        coords[1][1] = 0.0000000000001; 
                        coords[2][0] = 0.0000000000001;
                        coords[2][1] = 0.0000000000001;
                        coords[3][0] = 0.0573330873026;
                        coords[3][1] = 0.0151382269814;
                        coords[4][0] = 0.0573330873026;
                        coords[4][1] = 0.9275286857160; 
                        coords[5][0] = 0.9275286857160;
                        coords[5][1] = 0.0573330873026;
                        coords[6][0] = 0.0151382269814; 
                        coords[6][1] = 0.0573330873026; 
                        coords[7][0] = 0.9275286857160;
                        coords[7][1] = 0.0151382269814;
                        coords[8][0] = 0.0151382269814;
                        coords[8][1] = 0.9275286857160;
                        coords[9][0] = 0.8159625040711;
                        coords[9][1] = 0.1659719969565;
                        coords[10][0]  =0.8159625040711; 
                        coords[10][1] = 0.0180654989724;
                        coords[11][0] = 0.1659719969565;
                        coords[11][1] = 0.8159625040711;
                        coords[12][0] = 0.0180654989724;
                        coords[12][1] = 0.8159625040711;
                        coords[13][0] = 0.1659719969565;
                        coords[13][1] = 0.0180654989724;
                        coords[14][0] = 0.0180654989724;
                        coords[14][1] = 0.1659719969565;
                        coords[15][0] = 0.3165475556378;
                        coords[15][1] = 0.0186886898773;
                        coords[16][0] = 0.6647637544849;
                        coords[16][1] = 0.0186886898773;
                        coords[17][0] = 0.0186886898773; 
                        coords[17][1] = 0.6647637544849;
                        coords[18][0] = 0.0186886898773; 
                        coords[18][1] = 0.3165475556378;
                        coords[19][0] = 0.3165475556378;
                        coords[19][1] = 0.6647637544849;
                        coords[20][0] = 0.6647637544849; 
                        coords[20][1] = 0.3165475556378;
                        coords[21][0] = 0.0192662192492; 
                        coords[21][1] = 0.4903668903754;
                        coords[22][0] = 0.4903668903754; 
                        coords[22][1] = 0.0192662192492;
                        coords[23][0] = 0.4903668903754;
                        coords[23][1] = 0.4903668903754;
                        coords[24][0] = 0.0875134669581;
                        coords[24][1] = 0.8249730660837;
                        coords[25][0] = 0.0875134669581; 
                        coords[25][1] = 0.0875134669581;
                        coords[26][0] = 0.8249730660837; 
                        coords[26][1] = 0.0875134669581;
                        coords[27][0] = 0.0935526036219; 
                        coords[27][1] = 0.2079865423167;
                        coords[28][0] = 0.0935526036219; 
                        coords[28][1] = 0.6984608540613;
                        coords[29][0] = 0.2079865423167; 
                        coords[29][1] = 0.0935526036219;
                        coords[30][0] = 0.6984608540613;
                        coords[30][1] = 0.0935526036219;
                        coords[31][0] = 0.6984608540613;
                        coords[31][1] = 0.2079865423167;
                        coords[32][0] = 0.2079865423167;
                        coords[32][1] = 0.6984608540613;
                        coords[33][0] = 0.0974892983467; 
                        coords[33][1] = 0.5380088595149;
                        coords[34][0] = 0.3645018421383; 
                        coords[34][1] = 0.0974892983467;
                        coords[35][0] = 0.5380088595149;
                        coords[35][1] = 0.0974892983467;
                        coords[36][0] = 0.5380088595149; 
                        coords[36][1] = 0.3645018421383; 
                        coords[37][0] = 0.3645018421383; 
                        coords[37][1] = 0.5380088595149; 
                        coords[38][0] = 0.0974892983467;
                        coords[38][1] = 0.3645018421383;
                        coords[39][0] = 0.2217145894873;
                        coords[39][1] = 0.5565708210253;
                        coords[40][0] = 0.5565708210253;
                        coords[40][1] = 0.2217145894873; 
                        coords[41][0] = 0.2217145894873;
                        coords[41][1] = 0.2217145894873;
                        coords[42][0] = 0.3860471669296; 
                        coords[42][1] = 0.2279056661408; 
                        coords[43][0] = 0.2279056661408; 
                        coords[43][1] = 0.3860471669296;
                        coords[44][0] = 0.3860471669296;
                        coords[44][1] = 0.3860471669296;                      
      } break;
          case CFPolyOrder::ORDER9:
      {
                        coords[0][0] = 0.9999999999999;
                        coords[0][1] = 0.0000000000001; 
                        coords[1][0] = 0.0000000000001;
                        coords[1][1] = 0.9999999999999; 
                        coords[2][0] = 0.0000000000001;
                        coords[2][1] = 0.0000000000001;
                        coords[3][0] = 0.9398863583577;
                        coords[3][1] = 0.0049848744634;
                        coords[4][0] =0.0543806683058;
                        coords[4][1] = 0.9386405618617; 
                        coords[5][0] = 0.0093940049164;
                        coords[5][1] = 0.0526424462697;
                        coords[6][0] = 0.0164345086362; 
                        coords[6][1] = 0.9469035517351; 
                        coords[7][0] = 0.9469487269862;
                        coords[7][1] = 0.0363373677167;
                        coords[8][0] = 0.0426604005768;
                        coords[8][1] = 0.0151224541799;
                        coords[9][0] = 0.0122269495439;
                        coords[9][1] = 0.8693773510664;
                        coords[10][0] = 0.8673696521047; 
                        coords[10][1] = 0.1204917285774;
                        coords[11][0] = 0.8456744021389;
                        coords[11][1] = 0.0157763967870;
                        coords[12][0] = 0.1395759632103;
                        coords[12][1] = 0.8448120870375;
                        coords[13][0] = 0.1317821743231;
                        coords[13][1] = 0.0135009605584;
                        coords[14][0] = 0.0157955126300;
                        coords[14][1] = 0.1455274938536;
                        coords[15][0] = 0.7365462884436;
                        coords[15][1] = 0.0155697540908;
                        coords[16][0] = 0.0139688430330;
                        coords[16][1] = 0.7379836894450;
                        coords[17][0] = 0.2547895186039; 
                        coords[17][1] = 0.7297615689771;
                        coords[18][0] = 0.7316386522555;
                        coords[18][1] = 0.2543076683315;
                        coords[19][0] = 0.0157253728951;
                        coords[19][1] = 0.2696239795791;
                        coords[20][0] = 0.2662302843647; 
                        coords[20][1] = 0.0144783956308;
                        coords[21][0] = 0.8673504065214; 
                        coords[21][1] = 0.0591679410400;
                        coords[22][0] = 0.0741493666957; 
                        coords[22][1] = 0.8634782575061;
                        coords[23][0] = 0.0159285948360;
                        coords[23][1] = 0.4191238955238;
                        coords[24][0] = 0.0156061028068;
                        coords[24][1] = 0.5809222921146;
                        coords[25][0] = 0.5910094817484;
                        coords[25][1] = 0.0159251452651;
                        coords[26][0] = 0.4034771496889; 
                        coords[26][1] = 0.5806700368104;
                        coords[27][0] = 0.5694745628526; 
                        coords[27][1] = 0.4149495146302;
                        coords[28][0] = 0.0678493700650;
                        coords[28][1] = 0.0761218678591;
                        coords[29][0] = 0.4265968590272; 
                        coords[29][1] = 0.0157509692312;
                        coords[30][0] = 0.0670982507890;
                        coords[30][1] = 0.7741898312421;
                        coords[31][0] = 0.7528310231480;
                        coords[31][1] = 0.0819119495639;
                        coords[32][0] = 0.7753727783557;
                        coords[32][1] = 0.1577128457292;
                        coords[33][0] = 0.1689073157787; 
                        coords[33][1] = 0.7503943099742;
                        coords[34][0] = 0.1687335832919; 
                        coords[34][1] = 0.0708311507268;
                        coords[35][0] = 0.0821244708436;
                        coords[35][1] = 0.1762996626771;
                        coords[36][0] = 0.6288705363345; 
                        coords[36][1] = 0.0807744953317; 
                        coords[37][0] = 0.0811413015266; 
                        coords[37][1] = 0.3054373589776; 
                        coords[38][0] = 0.2969112065080;
                        coords[38][1] = 0.6227485988871;
                        coords[39][0] = 0.0767542314171;
                        coords[39][1] = 0.6247247149546;
                        coords[40][0] = 0.6223022333845;
                        coords[40][1] = 0.3011485821166; 
                        coords[41][0] = 0.3103786288051;
                        coords[41][1] = 0.0779098365079;
                        coords[42][0] = 0.0819218215187; 
                        coords[42][1] = 0.4603633038351; 
                        coords[43][0] = 0.4717022665013; 
                        coords[43][1] = 0.0821554006797;
                        coords[44][0] = 0.4546603415250;
                        coords[44][1] = 0.4637565033890;
                        coords[45][0] = 0.1701091339237;
                        coords[45][1] = 0.6422277808188; 
                        coords[46][0] = 0.6406004329487;
                        coords[46][1] = 0.1898293537256; 
                        coords[47][0] = 0.1912267583717;
                        coords[47][1] = 0.1739955685343;
                        coords[48][0] = 0.1885315767070;                        
                        coords[48][1] = 0.4798914070406;
                        coords[49][0] = 0.4772929957691;
                        coords[49][1] = 0.3348356598119; 
                        coords[50][0] = 0.3126974621760;
                        coords[50][1] = 0.4957972197259;
                        coords[51][0] = 0.4961225945946; 
                        coords[51][1] = 0.1927553668904; 
                        coords[52][0] = 0.1928805312867;
                        coords[52][1] = 0.3161015807261;
                        coords[53][0] = 0.3360041453816;
                        coords[53][1] = 0.1894892801290;
                        coords[54][0] = 0.3337280550848;
                        coords[54][1] = 0.3343571021811;
                      
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Quadrature Points not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
    }
  

  coordsFluxPointsInit.resize(polyOrder+1);
switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
	  coordsFluxPointsInit[0] = 0.0;
      } break;
      case CFPolyOrder::ORDER1:
      {
	  coordsFluxPointsInit[0] = -1./sqrt(3.);
	  coordsFluxPointsInit[1] = +1./sqrt(3.);
      } break;
      case CFPolyOrder::ORDER2:
      {
	  coordsFluxPointsInit[0] = -sqrt(3./5.);
	  coordsFluxPointsInit[1] = 0.0;
	  coordsFluxPointsInit[2] = +sqrt(3./5.);
      } break;
      case CFPolyOrder::ORDER3:
      {
	coordsFluxPointsInit[0] = -sqrt(3./7.+2./7.*sqrt(6./5.));
	coordsFluxPointsInit[1] = -sqrt(3./7.-2./7.*sqrt(6./5.));
	coordsFluxPointsInit[2] = +sqrt(3./7.-2./7.*sqrt(6./5.));
	coordsFluxPointsInit[3] = +sqrt(3./7.+2./7.*sqrt(6./5.));
      } break;
      case CFPolyOrder::ORDER4:
      {
	coordsFluxPointsInit[0] = -sqrt(5.+2.*sqrt(10./7.))/3.;
	coordsFluxPointsInit[1] = -sqrt(5.-2.*sqrt(10./7.))/3.;
	coordsFluxPointsInit[2] = 0.0;
	coordsFluxPointsInit[3] = +sqrt(5.-2.*sqrt(10./7.))/3.;
	coordsFluxPointsInit[4] = +sqrt(5.+2.*sqrt(10./7.))/3.;
      } break;
      case CFPolyOrder::ORDER5:
      {
        coordsFluxPointsInit[0] = -0.9324695142031521;
	coordsFluxPointsInit[1] = -0.6612093864662645;
	coordsFluxPointsInit[2] = -0.2386191860831969;
	coordsFluxPointsInit[3] = 0.2386191860831969;
	coordsFluxPointsInit[4] = 0.6612093864662645;
	coordsFluxPointsInit[5] = 0.9324695142031521;
      } break;
      case CFPolyOrder::ORDER6:
      {
        coordsFluxPointsInit[0] = -0.949107912342759;
	coordsFluxPointsInit[1] = -0.741531185599394;
	coordsFluxPointsInit[2] = -0.405845151377397;
	coordsFluxPointsInit[3] = 0.0;
	coordsFluxPointsInit[4] = 0.405845151377397;
	coordsFluxPointsInit[5] = 0.741531185599394;
	coordsFluxPointsInit[6] = 0.949107912342759;
      } break;
      case CFPolyOrder::ORDER7:
      {
        coordsFluxPointsInit[0] = -0.960289856497536;
	coordsFluxPointsInit[1] = -0.796666477413627;
	coordsFluxPointsInit[2] = -0.525532409916329;
	coordsFluxPointsInit[3] = -0.183434642495650;
	coordsFluxPointsInit[4] = 0.183434642495650;
	coordsFluxPointsInit[5] = 0.525532409916329;
	coordsFluxPointsInit[6] = 0.796666477413627;
	coordsFluxPointsInit[7] = 0.960289856497536;
      } break;
      case CFPolyOrder::ORDER8:
      {
        coordsFluxPointsInit[0] = -0.968160239507626;
	coordsFluxPointsInit[1] = -0.836031107326636;
	coordsFluxPointsInit[2] = -0.613371432700590;
	coordsFluxPointsInit[3] = -0.324253423403809;
	coordsFluxPointsInit[4] = 0.0;
	coordsFluxPointsInit[5] = 0.324253423403809;
	coordsFluxPointsInit[6] = 0.613371432700590;
	coordsFluxPointsInit[7] = 0.836031107326636;
	coordsFluxPointsInit[8] = 0.968160239507626;
      } break;
      case CFPolyOrder::ORDER9:
      {
        coordsFluxPointsInit[0] = -0.973906528517172;
	coordsFluxPointsInit[1] = -0.865063366688985;
	coordsFluxPointsInit[2] = -0.679409568299024;
	coordsFluxPointsInit[3] = -0.433395394129247;
	coordsFluxPointsInit[4] = -0.148874338981631;
	coordsFluxPointsInit[5] = 0.148874338981631;
	coordsFluxPointsInit[6] = 0.433395394129247;
	coordsFluxPointsInit[7] = 0.679409568299024;
	coordsFluxPointsInit[8] = 0.865063366688985;
	coordsFluxPointsInit[9] = 0.973906528517172;
      } break;
      case CFPolyOrder::ORDER10:
      {
        coordsFluxPointsInit[0] = -0.978228658146057;
	coordsFluxPointsInit[1] = -0.887062599768095;
	coordsFluxPointsInit[2] = -0.730152005574049;
	coordsFluxPointsInit[3] = -0.519096129110681;
	coordsFluxPointsInit[4] = -0.269543155952345;
	coordsFluxPointsInit[5] = 0.0;
	coordsFluxPointsInit[6] = 0.269543155952345;
	coordsFluxPointsInit[7] = 0.519096129110681;
	coordsFluxPointsInit[8] = 0.730152005574049;
	coordsFluxPointsInit[9] = 0.887062599768095;
	coordsFluxPointsInit[10] = 0.978228658146057;
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Gauss Legendre not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
}
  
  solPntsLocalCoord2D = coords;


  m_flxPntsLocalCoord1D = coordsFluxPointsInit;

  flxPntsLocalCoord2D = getLocalCoords2D(polyOrder);

  resetFluxReconstructionElementData();
}

//////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
TriagFluxReconstructionElementData::TriagFluxReconstructionElementData(CFPolyOrder::Type polyOrder, 
								       Common::SafePtr< BasePointDistribution > solPntDist, 
								       Common::SafePtr< BasePointDistribution > flxPntDist)
{
  m_shape = CFGeoShape::TRIAG;
  m_dimensionality = DIM_2D;
  m_polyOrder = polyOrder;
  CFuint rows = (0.5*(polyOrder+1)*(polyOrder+2)) ;
  std::vector< std::vector< CFreal > > coords( 2 , std::vector<CFreal> (rows)); 
  coords.resize(rows, std::vector< CFreal > (2));
  solPntsLocalCoord2D.resize(rows, std::vector< CFreal > (2));
  CFuint nbfaces = 3  ;
  flxPntsLocalCoord2D.resize(nbfaces, std::vector< std::vector<CFreal> >(polyOrder+1));
  for (CFuint iFace=0; iFace <nbfaces; ++iFace){
    for (CFuint iFlx=0; iFlx<polyOrder+1; ++iFlx){
		  (flxPntsLocalCoord2D[iFace][iFlx]).resize(0);
	  }
  }
	for (CFuint iFace=0; iFace <nbfaces; ++iFace){
    for (CFuint iFlx=0; iFlx<polyOrder+1; ++iFlx){
		  flxPntsLocalCoord2D[iFace][iFlx].push_back(0.);
		  flxPntsLocalCoord2D[iFace][iFlx].push_back(0.);
	  }
  }


switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
                        coords[0][0] = 0.333333333333333; 
                        coords[0][1] = 0.333333333333333;  
      } break;
      case CFPolyOrder::ORDER1:
      {
                        coords[0][0] = 0.1666666666667; 
                        coords[0][1] = 0.1666666666667;
                        coords[1][0] = 0.6666666666667; 
                        coords[1][1] = 0.1666666666667;  
                        coords[2][0] = 0.1666666666667; 
                        coords[2][1] = 0.6666666666667;

    } break;

      case CFPolyOrder::ORDER2:              
      {
                        coords[0][0] =  0.091576213509780; 
                        coords[0][1] =  0.091576213509780; 
                        
                        coords[1][0] =  0.816847572980440; 
                        coords[1][1] =  0.091576213509780;  
                        
                        coords[2][0] =  0.091576213509780; 
                        coords[2][1] =  0.816847572980440;
                        
                        coords[3][0] =  0.445948490915964; 
                        coords[3][1] =  0.108103018168071; 
                        
                        coords[4][0] =  0.108103018168071; 
                        coords[4][1] =  0.445948490915964;
                        
                        coords[5][0] =  0.445948490915964; 
                        coords[5][1] =  0.445948490915964;


      } break;

      case CFPolyOrder::ORDER3:          
      {
                        coords[0][0] = 0.055564052669793;
                        coords[0][1] = 0.055564052669793;
                        
                        coords[1][0] = 0.888871894660413; 
                        coords[1][1] = 0.055564052669793;	
                        
                        coords[2][0] = 0.055564052669793; 
                        coords[2][1] = 0.888871894660413;  	

                        
                        coords[3][0] = 0.634210747745723; 
                        coords[3][1] = 0.070255540518384;  	

                        
                        coords[4][0] = 0.070255540518384; 
                        coords[4][1] = 0.634210747745723;	
                        
                        coords[5][0] = 0.295533711735893; 
                        coords[5][1] = 0.634210747745723; 	

                        
                        coords[6][0] = 0.295533711735893; 
                        coords[6][1] = 0.070255540518384;  	
                        
                        coords[7][0] = 0.070255540518384; 
                        coords[7][1] = 0.295533711735893; 	
                        
                        coords[8][0] = 0.634210747745723; 
                        coords[8][1] = 0.295533711735893; 	
                        
                        coords[9][0] = 0.333333333333333; 
                        coords[9][1] = 0.333333333333333; 


      } break;
      case CFPolyOrder::ORDER4:   
      {
                        coords[0][0] = 0.035870877695734; 
                        coords[0][1] = 0.035870877695734;

                        coords[1][0] = 0.928258244608533; 
                        coords[1][1] = 0.035870877695734;

                        coords[2][0] = 0.035870877695734; 
                        coords[2][1] = 0.928258244608533; 

                        coords[3][0] = 0.241729395767967;
                        coords[3][1] = 0.241729395767967;

                        coords[4][0] = 0.516541208464066; 
                        coords[4][1] = 0.241729395767967; 

                        coords[5][0] = 0.241729395767967; 
                        coords[5][1] = 0.516541208464066; 

                        coords[6][0] = 0.474308787777079;
                        coords[6][1] = 0.051382424445843; 

                        coords[7][0] = 0.051382424445843;
                        coords[7][1] = 0.474308787777079;

                        coords[8][0] = 0.474308787777079; 
                        coords[8][1] = 0.474308787777079;

                        coords[9][0] = 0.751183631106484; 
                        coords[9][1] = 0.047312487011716;

                        coords[10][0] = 0.047312487011716; 
                        coords[10][1] = 0.751183631106484;

                        coords[11][0] =  0.201503881881800;
                        coords[11][1] =  0.751183631106484;

                        coords[12][0] = 0.201503881881800; 
                        coords[12][1] =  0.047312487011716;

                        coords[13][0] = 0.047312487011716;
                        coords[13][1] = 0.201503881881800;

                        coords[14][0] = 0.751183631106484; 
                        coords[14][1] = 0.201503881881800;

      } break;

      case CFPolyOrder::ORDER5:
      {
                        coords[0][0] = 0.028112952182664; 
                        coords[0][1] = 0.028112952182664; 

                        coords[1][0] = 0.943774095634672;
                        coords[1][1] = 0.028112952182664; 

                        coords[2][0] = 0.028112952182664;
                        coords[2][1] = 0.943774095634672;

                        coords[3][0] = 0.177139098469317;
                        coords[3][1] = 0.177139098469317;

                        coords[4][0] = 0.645721803061365; 
                        coords[4][1] = 0.177139098469317;

                        coords[5][0] = 0.177139098469317;
                        coords[5][1] = 0.645721803061365;

                        coords[6][0] = 0.405508595867433;
                        coords[6][1] = 0.188982808265134;

                        coords[7][0] = 0.188982808265134;
                        coords[7][1] = 0.405508595867433;

                        coords[8][0] = 0.405508595867433;
                        coords[8][1] = 0.405508595867433;

                        coords[9][0] = 0.817900980028499; 
                        coords[9][1] = 0.033533207700614;

                        coords[10][0] = 0.033533207700614; 
                        coords[10][1] = 0.817900980028499;

                        coords[11][0] = 0.148565812270887;
                        coords[11][1] = 0.817900980028499;

                        coords[12][0] = 0.148565812270887; 
                        coords[12][1] = 0.033533207700614; 

                        coords[13][0] = 0.033533207700614;
                        coords[13][1] = 0.148565812270887;

                        coords[14][0] = 0.817900980028499;
                        coords[14][1] = 0.148565812270887;

                        coords[15][0] = 0.604978911775132; 
                        coords[15][1] = 0.037824789609186;

                        coords[16][0] = 0.037824789609186; 
                        coords[16][1] = 0.604978911775132;

                        coords[17][0] = 0.357196298615681;
                        coords[17][1] = 0.604978911775132;

                        coords[18][0] = 0.357196298615681; 
                        coords[18][1] = 0.037824789609186;

                        coords[19][0] = 0.037824789609186; 
                        coords[19][1] = 0.357196298615681; 

                        coords[20][0] = 0.604978911775132; 
                        coords[20][1] = 0.357196298615681;
      } break;


   case CFPolyOrder::ORDER6:
      {  
                        coords[0][0] =  0.0000000000000;
                        coords[0][1] =  0.9451704450174; 
                        coords[1][0] =  0.9451704450173; 
                        coords[1][1] =  0.0000000000000; 
                        coords[2][0] =  0.9289002405719;
                        coords[2][1] =  0.0685505797224;
                        coords[3][0] =  0.0685505797224;
                        coords[3][1] =  0.9289002405717;
                        coords[4][0] =  0.0243268355615;
                        coords[4][1] =  0.0243268355616; 
                        coords[5][0] =  0.1279662835335;
                        coords[5][1] =  0.0277838749488;
                        coords[6][0] =  0.0277838749488; 
                        coords[6][1] =  0.1279662835337; 
                        coords[7][0] =  0.0287083428360; 
                        coords[7][1] =  0.7498347588657;
                        coords[8][0] =  0.7498347588656;
                        coords[8][1] =  0.0287083428360;
                        coords[9][0] =  0.7228007909707;
                        coords[9][1] =  0.2497602062385;
                        coords[10][0] = 0.2497602062386; 
                        coords[10][1] = 0.7228007909707;
                        coords[11][0] = 0.0865562992839;
                        coords[11][1] = 0.8325513856997;
                        coords[12][0] = 0.8325513856998;
                        coords[12][1] = 0.0865562992839;
                        coords[13][0] = 0.3061619157672;
                        coords[13][1] = 0.0303526617491;
                        coords[14][0] = 0.0303526617491;
                        coords[14][1] = 0.3061619157675;
                        coords[15][0] = 0.4868610595047;
                        coords[15][1] = 0.4868610595047;
                        coords[16][0] = 0.6657904293017;
                        coords[16][1] = 0.1765456154219;
                        coords[17][0] = 0.1765456154221; 
                        coords[17][1] = 0.6657904293016;
                        coords[18][0] = 0.0293121007360; 
                        coords[18][1] = 0.5295657488669;
                        coords[19][0] = 0.5295657488667;
                        coords[19][1] = 0.0293121007360;
                        coords[20][0] = 0.1444673824391; 
                        coords[20][1] = 0.1444673824391;
                        coords[21][0] = 0.3299740111411; 
                        coords[21][1] = 0.5361815729050;
                        coords[22][0] = 0.5361815729052; 
                        coords[22][1] = 0.3299740111409;
                        coords[23][0] = 0.5511507516862;
                        coords[23][1] = 0.1437790861923;
                        coords[24][0] = 0.1437790861923;
                        coords[24][1] = 0.5511507516862;
                        coords[25][0] = 0.3348066587327; 
                        coords[25][1] = 0.1529619437161;
                        coords[26][0] = 0.1529619437161; 
                        coords[26][1] = 0.3348066587327;
                        coords[27][0] = 0.3430183498147; 
                        coords[27][1] = 0.3430183498147;
                      
      } break;


      case CFPolyOrder::ORDER7:
      {
                        coords[0][0] = 0.0242935351590;
                        coords[0][1] = 0.9493059293846; 
                        coords[1][0] = 0.0265193427722;
                        coords[1][1] = 0.0242695130640; 
                        coords[2][0] = 0.9492126023551;
                        coords[2][1] = 0.0265067966437;
                        coords[3][0] = 0.0033775763749;
                        coords[3][1] = 0.4767316412363;
                        coords[4][0] = 0.4757672298101;
                        coords[4][1] = 0.5198921829102; 
                        coords[5][0] = 0.5190783193471;
                        coords[5][1] = 0.0055912706202;
                        coords[6][0] = 0.8616839745321; 
                        coords[6][1] = 0.0133996048618; 
                        coords[7][0] = 0.1249209759926;
                        coords[7][1] = 0.8613054321334;
                        coords[8][0] = 0.0138565453861;
                        coords[8][1] = 0.1247733717358;
                        coords[9][0] = 0.0211887064222;
                        coords[9][1] = 0.8438438351223;
                        coords[10][0]  =0.8432296787219; 
                        coords[10][1] = 0.1354563645830;
                        coords[11][0] = 0.1354231797865;
                        coords[11][1] = 0.0213482820656;
                        coords[12][0] = 0.3088853510679;
                        coords[12][1] = 0.0221919663014;
                        coords[13][0] = 0.6685057595169;
                        coords[13][1] = 0.3089012879389;
                        coords[14][0] = 0.0226545012557;
                        coords[14][1] = 0.6691709943321;
                        coords[15][0] = 0.2808515408772;
                        coords[15][1] = 0.6924718155106;
                        coords[16][0] = 0.6922446749051;
                        coords[16][1] = 0.0268723345026;
                        coords[17][0] = 0.0268617447119; 
                        coords[17][1] = 0.2810093973222;
                        coords[18][0] = 0.1141778485470; 
                        coords[18][1] = 0.7973581413586;
                        coords[19][0] = 0.7974807922061;
                        coords[19][1] = 0.0879806508791;
                        coords[20][0] = 0.0892807293894; 
                        coords[20][1] = 0.1145020561128;
                        coords[21][0] = 0.1052487892455; 
                        coords[21][1] = 0.6686904119922;
                        coords[22][0] = 0.6663022280740; 
                        coords[22][1] = 0.2275051631832;
                        coords[23][0] = 0.2307803737547;
                        coords[23][1] = 0.1054572561221;
                        coords[24][0] = 0.1705059157540;
                        coords[24][1] = 0.5174064398658;
                        coords[25][0] = 0.5086593973043; 
                        coords[25][1] = 0.3170523855209;
                        coords[26][0] = 0.3141823862281; 
                        coords[26][1] = 0.1810706361659;
                        coords[27][0] = 0.4617460817864; 
                        coords[27][1] = 0.4678594539804;
                        coords[28][0] = 0.0693087496081; 
                        coords[28][1] = 0.4622856042085;
                        coords[29][0] = 0.4622856042085; 
                        coords[29][1] = 0.0724357805669;
                        coords[30][0] = 0.2578625857893;
                        coords[30][1] = 0.6131395039177;
                        coords[31][0] = 0.6112627766779;
                        coords[31][1] = 0.1300360834609;
                        coords[32][0] = 0.1305182135934;
                        coords[32][1] = 0.2581713828884;
                        coords[33][0] = 0.4281437991828; 
                        coords[33][1] = 0.2362005969817;
                        coords[34][0] = 0.3356995783730; 
                        coords[34][1] = 0.4311026308588;
                        coords[35][0] = 0.2305424298836; 
                        coords[35][1] = 0.3456013949376;

                      
      } break;

        case CFPolyOrder::ORDER8:
      {
                        coords[0][0] = 0.0000000000001;
                        coords[0][1] = 0.9999999999999; 
                        coords[1][0] = 0.9999999999999;
                        coords[1][1] = 0.0000000000001; 
                        coords[2][0] = 0.0000000000001;
                        coords[2][1] = 0.0000000000001;
                        coords[3][0] = 0.0573330873026;
                        coords[3][1] = 0.0151382269814;
                        coords[4][0] = 0.0573330873026;
                        coords[4][1] = 0.9275286857160; 
                        coords[5][0] = 0.9275286857160;
                        coords[5][1] = 0.0573330873026;
                        coords[6][0] = 0.0151382269814; 
                        coords[6][1] = 0.0573330873026; 
                        coords[7][0] = 0.9275286857160;
                        coords[7][1] = 0.0151382269814;
                        coords[8][0] = 0.0151382269814;
                        coords[8][1] = 0.9275286857160;
                        coords[9][0] = 0.8159625040711;
                        coords[9][1] = 0.1659719969565;
                        coords[10][0]  =0.8159625040711; 
                        coords[10][1] = 0.0180654989724;
                        coords[11][0] = 0.1659719969565;
                        coords[11][1] = 0.8159625040711;
                        coords[12][0] = 0.0180654989724;
                        coords[12][1] = 0.8159625040711;
                        coords[13][0] = 0.1659719969565;
                        coords[13][1] = 0.0180654989724;
                        coords[14][0] = 0.0180654989724;
                        coords[14][1] = 0.1659719969565;
                        coords[15][0] = 0.3165475556378;
                        coords[15][1] = 0.0186886898773;
                        coords[16][0] = 0.6647637544849;
                        coords[16][1] = 0.0186886898773;
                        coords[17][0] = 0.0186886898773; 
                        coords[17][1] = 0.6647637544849;
                        coords[18][0] = 0.0186886898773; 
                        coords[18][1] = 0.3165475556378;
                        coords[19][0] = 0.3165475556378;
                        coords[19][1] = 0.6647637544849;
                        coords[20][0] = 0.6647637544849; 
                        coords[20][1] = 0.3165475556378;
                        coords[21][0] = 0.0192662192492; 
                        coords[21][1] = 0.4903668903754;
                        coords[22][0] = 0.4903668903754; 
                        coords[22][1] = 0.0192662192492;
                        coords[23][0] = 0.4903668903754;
                        coords[23][1] = 0.4903668903754;
                        coords[24][0] = 0.0875134669581;
                        coords[24][1] = 0.8249730660837;
                        coords[25][0] = 0.0875134669581; 
                        coords[25][1] = 0.0875134669581;
                        coords[26][0] = 0.8249730660837; 
                        coords[26][1] = 0.0875134669581;
                        coords[27][0] = 0.0935526036219; 
                        coords[27][1] = 0.2079865423167;
                        coords[28][0] = 0.0935526036219; 
                        coords[28][1] = 0.6984608540613;
                        coords[29][0] = 0.2079865423167; 
                        coords[29][1] = 0.0935526036219;
                        coords[30][0] = 0.6984608540613;
                        coords[30][1] = 0.0935526036219;
                        coords[31][0] = 0.6984608540613;
                        coords[31][1] = 0.2079865423167;
                        coords[32][0] = 0.2079865423167;
                        coords[32][1] = 0.6984608540613;
                        coords[33][0] = 0.0974892983467; 
                        coords[33][1] = 0.5380088595149;
                        coords[34][0] = 0.3645018421383; 
                        coords[34][1] = 0.0974892983467;
                        coords[35][0] = 0.5380088595149;
                        coords[35][1] = 0.0974892983467;
                        coords[36][0] = 0.5380088595149; 
                        coords[36][1] = 0.3645018421383; 
                        coords[37][0] = 0.3645018421383; 
                        coords[37][1] = 0.5380088595149; 
                        coords[38][0] = 0.0974892983467;
                        coords[38][1] = 0.3645018421383;
                        coords[39][0] = 0.2217145894873;
                        coords[39][1] = 0.5565708210253;
                        coords[40][0] = 0.5565708210253;
                        coords[40][1] = 0.2217145894873; 
                        coords[41][0] = 0.2217145894873;
                        coords[41][1] = 0.2217145894873;
                        coords[42][0] = 0.3860471669296; 
                        coords[42][1] = 0.2279056661408; 
                        coords[43][0] = 0.2279056661408; 
                        coords[43][1] = 0.3860471669296;
                        coords[44][0] = 0.3860471669296;
                        coords[44][1] = 0.3860471669296;                      
      } break;
          case CFPolyOrder::ORDER9:
      {
                        coords[0][0] = 0.9999999999999;
                        coords[0][1] = 0.0000000000001; 
                        coords[1][0] = 0.0000000000001;
                        coords[1][1] = 0.9999999999999; 
                        coords[2][0] = 0.0000000000001;
                        coords[2][1] = 0.0000000000001;
                        coords[3][0] = 0.9398863583577;
                        coords[3][1] = 0.0049848744634;
                        coords[4][0] =0.0543806683058;
                        coords[4][1] = 0.9386405618617; 
                        coords[5][0] = 0.0093940049164;
                        coords[5][1] = 0.0526424462697;
                        coords[6][0] = 0.0164345086362; 
                        coords[6][1] = 0.9469035517351; 
                        coords[7][0] = 0.9469487269862;
                        coords[7][1] = 0.0363373677167;
                        coords[8][0] = 0.0426604005768;
                        coords[8][1] = 0.0151224541799;
                        coords[9][0] = 0.0122269495439;
                        coords[9][1] = 0.8693773510664;
                        coords[10][0] = 0.8673696521047; 
                        coords[10][1] = 0.1204917285774;
                        coords[11][0] = 0.8456744021389;
                        coords[11][1] = 0.0157763967870;
                        coords[12][0] = 0.1395759632103;
                        coords[12][1] = 0.8448120870375;
                        coords[13][0] = 0.1317821743231;
                        coords[13][1] = 0.0135009605584;
                        coords[14][0] = 0.0157955126300;
                        coords[14][1] = 0.1455274938536;
                        coords[15][0] = 0.7365462884436;
                        coords[15][1] = 0.0155697540908;
                        coords[16][0] = 0.0139688430330;
                        coords[16][1] = 0.7379836894450;
                        coords[17][0] = 0.2547895186039; 
                        coords[17][1] = 0.7297615689771;
                        coords[18][0] = 0.7316386522555;
                        coords[18][1] = 0.2543076683315;
                        coords[19][0] = 0.0157253728951;
                        coords[19][1] = 0.2696239795791;
                        coords[20][0] = 0.2662302843647; 
                        coords[20][1] = 0.0144783956308;
                        coords[21][0] = 0.8673504065214; 
                        coords[21][1] = 0.0591679410400;
                        coords[22][0] = 0.0741493666957; 
                        coords[22][1] = 0.8634782575061;
                        coords[23][0] = 0.0159285948360;
                        coords[23][1] = 0.4191238955238;
                        coords[24][0] = 0.0156061028068;
                        coords[24][1] = 0.5809222921146;
                        coords[25][0] = 0.5910094817484;
                        coords[25][1] = 0.0159251452651;
                        coords[26][0] = 0.4034771496889; 
                        coords[26][1] = 0.5806700368104;
                        coords[27][0] = 0.5694745628526; 
                        coords[27][1] = 0.4149495146302;
                        coords[28][0] = 0.0678493700650;
                        coords[28][1] = 0.0761218678591;
                        coords[29][0] = 0.4265968590272; 
                        coords[29][1] = 0.0157509692312;
                        coords[30][0] = 0.0670982507890;
                        coords[30][1] = 0.7741898312421;
                        coords[31][0] = 0.7528310231480;
                        coords[31][1] = 0.0819119495639;
                        coords[32][0] = 0.7753727783557;
                        coords[32][1] = 0.1577128457292;
                        coords[33][0] = 0.1689073157787; 
                        coords[33][1] = 0.7503943099742;
                        coords[34][0] = 0.1687335832919; 
                        coords[34][1] = 0.0708311507268;
                        coords[35][0] = 0.0821244708436;
                        coords[35][1] = 0.1762996626771;
                        coords[36][0] = 0.6288705363345; 
                        coords[36][1] = 0.0807744953317; 
                        coords[37][0] = 0.0811413015266; 
                        coords[37][1] = 0.3054373589776; 
                        coords[38][0] = 0.2969112065080;
                        coords[38][1] = 0.6227485988871;
                        coords[39][0] = 0.0767542314171;
                        coords[39][1] = 0.6247247149546;
                        coords[40][0] = 0.6223022333845;
                        coords[40][1] = 0.3011485821166; 
                        coords[41][0] = 0.3103786288051;
                        coords[41][1] = 0.0779098365079;
                        coords[42][0] = 0.0819218215187; 
                        coords[42][1] = 0.4603633038351; 
                        coords[43][0] = 0.4717022665013; 
                        coords[43][1] = 0.0821554006797;
                        coords[44][0] = 0.4546603415250;
                        coords[44][1] = 0.4637565033890;
                        coords[45][0] = 0.1701091339237;
                        coords[45][1] = 0.6422277808188; 
                        coords[46][0] = 0.6406004329487;
                        coords[46][1] = 0.1898293537256; 
                        coords[47][0] = 0.1912267583717;
                        coords[47][1] = 0.1739955685343;
                        coords[48][0] = 0.1885315767070;                        
                        coords[48][1] = 0.4798914070406;
                        coords[49][0] = 0.4772929957691;
                        coords[49][1] = 0.3348356598119; 
                        coords[50][0] = 0.3126974621760;
                        coords[50][1] = 0.4957972197259;
                        coords[51][0] = 0.4961225945946; 
                        coords[51][1] = 0.1927553668904; 
                        coords[52][0] = 0.1928805312867;
                        coords[52][1] = 0.3161015807261;
                        coords[53][0] = 0.3360041453816;
                        coords[53][1] = 0.1894892801290;
                        coords[54][0] = 0.3337280550848;
                        coords[54][1] = 0.3343571021811;
                      
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Quadrature Points not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
}


coordsFluxPointsInit.resize(polyOrder+1);
switch(polyOrder)
    {
      case CFPolyOrder::ORDER0:
      {
	  coordsFluxPointsInit[0] = 0.0;
      } break;
      case CFPolyOrder::ORDER1:
      {
	  coordsFluxPointsInit[0] = -1./sqrt(3.);
	  coordsFluxPointsInit[1] = +1./sqrt(3.);
      } break;
      case CFPolyOrder::ORDER2:
      {
	  coordsFluxPointsInit[0] = -sqrt(3./5.);
	  coordsFluxPointsInit[1] = 0.0;
	  coordsFluxPointsInit[2] = +sqrt(3./5.);
      } break;
      case CFPolyOrder::ORDER3:
      {
	coordsFluxPointsInit[0] = -sqrt(3./7.+2./7.*sqrt(6./5.));
	coordsFluxPointsInit[1] = -sqrt(3./7.-2./7.*sqrt(6./5.));
	coordsFluxPointsInit[2] = +sqrt(3./7.-2./7.*sqrt(6./5.));
	coordsFluxPointsInit[3] = +sqrt(3./7.+2./7.*sqrt(6./5.));
      } break;
      case CFPolyOrder::ORDER4:
      {
	coordsFluxPointsInit[0] = -sqrt(5.+2.*sqrt(10./7.))/3.;
	coordsFluxPointsInit[1] = -sqrt(5.-2.*sqrt(10./7.))/3.;
	coordsFluxPointsInit[2] = 0.0;
	coordsFluxPointsInit[3] = +sqrt(5.-2.*sqrt(10./7.))/3.;
	coordsFluxPointsInit[4] = +sqrt(5.+2.*sqrt(10./7.))/3.;
      } break;
      case CFPolyOrder::ORDER5:
      {
        coordsFluxPointsInit[0] = -0.9324695142031521;
	coordsFluxPointsInit[1] = -0.6612093864662645;
	coordsFluxPointsInit[2] = -0.2386191860831969;
	coordsFluxPointsInit[3] = 0.2386191860831969;
	coordsFluxPointsInit[4] = 0.6612093864662645;
	coordsFluxPointsInit[5] = 0.9324695142031521;
      } break;
      case CFPolyOrder::ORDER6:
      {
        coordsFluxPointsInit[0] = -0.949107912342759;
	coordsFluxPointsInit[1] = -0.741531185599394;
	coordsFluxPointsInit[2] = -0.405845151377397;
	coordsFluxPointsInit[3] = 0.0;
	coordsFluxPointsInit[4] = 0.405845151377397;
	coordsFluxPointsInit[5] = 0.741531185599394;
	coordsFluxPointsInit[6] = 0.949107912342759;
      } break;
      case CFPolyOrder::ORDER7:
      {
        coordsFluxPointsInit[0] = -0.960289856497536;
	coordsFluxPointsInit[1] = -0.796666477413627;
	coordsFluxPointsInit[2] = -0.525532409916329;
	coordsFluxPointsInit[3] = -0.183434642495650;
	coordsFluxPointsInit[4] = 0.183434642495650;
	coordsFluxPointsInit[5] = 0.525532409916329;
	coordsFluxPointsInit[6] = 0.796666477413627;
	coordsFluxPointsInit[7] = 0.960289856497536;
      } break;
      case CFPolyOrder::ORDER8:
      {
        coordsFluxPointsInit[0] = -0.968160239507626;
	coordsFluxPointsInit[1] = -0.836031107326636;
	coordsFluxPointsInit[2] = -0.613371432700590;
	coordsFluxPointsInit[3] = -0.324253423403809;
	coordsFluxPointsInit[4] = 0.0;
	coordsFluxPointsInit[5] = 0.324253423403809;
	coordsFluxPointsInit[6] = 0.613371432700590;
	coordsFluxPointsInit[7] = 0.836031107326636;
	coordsFluxPointsInit[8] = 0.968160239507626;
      } break;
      case CFPolyOrder::ORDER9:
      {
        coordsFluxPointsInit[0] = -0.973906528517172;
	coordsFluxPointsInit[1] = -0.865063366688985;
	coordsFluxPointsInit[2] = -0.679409568299024;
	coordsFluxPointsInit[3] = -0.433395394129247;
	coordsFluxPointsInit[4] = -0.148874338981631;
	coordsFluxPointsInit[5] = 0.148874338981631;
	coordsFluxPointsInit[6] = 0.433395394129247;
	coordsFluxPointsInit[7] = 0.679409568299024;
	coordsFluxPointsInit[8] = 0.865063366688985;
	coordsFluxPointsInit[9] = 0.973906528517172;
      } break;
      case CFPolyOrder::ORDER10:
      {
        coordsFluxPointsInit[0] = -0.978228658146057;
	coordsFluxPointsInit[1] = -0.887062599768095;
	coordsFluxPointsInit[2] = -0.730152005574049;
	coordsFluxPointsInit[3] = -0.519096129110681;
	coordsFluxPointsInit[4] = -0.269543155952345;
	coordsFluxPointsInit[5] = 0.0;
	coordsFluxPointsInit[6] = 0.269543155952345;
	coordsFluxPointsInit[7] = 0.519096129110681;
	coordsFluxPointsInit[8] = 0.730152005574049;
	coordsFluxPointsInit[9] = 0.887062599768095;
	coordsFluxPointsInit[10] = 0.978228658146057;
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Gauss Legendre not implemented for order "
                                      + StringOps::to_str(polyOrder) + ".");
      }
}

solPntsLocalCoord2D = coords;
m_flxPntsLocalCoord1D.resize(polyOrder+1);
m_flxPntsLocalCoord1D = coordsFluxPointsInit;

  flxPntsLocalCoord2D = getLocalCoords2D(polyOrder);


  //m_flxPntDistribution = flxPntDist;
  resetFluxReconstructionElementData();
}
//////////////////////////////////////////////////////////////////////



TriagFluxReconstructionElementData::~TriagFluxReconstructionElementData()
{
}

//////////////////////////////////////////////////////////////////////////////
CFreal TriagFluxReconstructionElementData::ComputeJacobi(CFuint N,CFuint alpha, CFuint beta,CFreal x)
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
CFreal TriagFluxReconstructionElementData::factorial(CFreal n)
{
  return (n==1. || n==0.) ? 1. : factorial(n-1.)*n;
}

//////////////////////////////////////////////////////////////////////

std::vector<CFreal> TriagFluxReconstructionElementData::getPercentage(CFPolyOrder::Type solOrder){
    CFAUTOTRACE;
    std::vector<CFreal> alpha;
    alpha.resize(solOrder+1);
    m_flxPntsLocalCoord1D.resize(solOrder+1);

    //std::vector< CFreal > m_flxPntsLocalCoord1D = getLocalCoords1D(solOrder);
    m_flxPntsLocalCoord1D = coordsFluxPointsInit;
    for (CFuint iFlx=0; iFlx< (solOrder+1) ; ++iFlx){
        alpha[iFlx] = (m_flxPntsLocalCoord1D[iFlx]+1)/2;
    }
		return alpha;
}
//////////////////////////////////////////////////////////////////////
std::vector<std::vector < std::vector < CFreal> > > TriagFluxReconstructionElementData::getLocalCoords2D(CFPolyOrder::Type solOrder)
{
  CFAUTOTRACE;
  std::vector<CFreal> alpha = getPercentage(solOrder);
  CFuint nbfaces = 3  ;
  std::vector< std::vector< std::vector < CFreal> > > coordsF( 2 , std::vector< std::vector<CFreal > > (solOrder+1)); 
  coordsF.resize(nbfaces, std::vector< std::vector<CFreal> >(solOrder+1));
  
  RealVector uFace1(2); uFace1[0] = 1.  ; uFace1[1] = 0. ; 
  RealVector uFace2(2); uFace2[0] = -1. ; uFace2[1] = 1. ; 
  RealVector uFace3(2); uFace3[0] = 0. ; uFace3[1] = -1. ; 
  RealVector start_Node1(2); start_Node1[0] = 0.; start_Node1[1] = 0.; 
  RealVector start_Node2(2); start_Node2[0] =  1.; start_Node2[1] = 0.; 
  RealVector start_Node3(2); start_Node3[0] =  0.; start_Node3[1] =  1.; 
  
  
// Structure of Flux points : Number of the Face, then the number of the Flx, then the KSI, ETA
	// Face 1
  for (int iFace=0; iFace <nbfaces; ++iFace){
    for (int iFlx=0; iFlx<solOrder+1; ++iFlx){
		  (coordsF[iFace][iFlx]).resize(0);
	  }
  }
	for (CFuint iFlx=0; iFlx<solOrder+1; ++iFlx){
		CFreal ksi_value = start_Node1[0]+alpha[iFlx]*uFace1[0];
		coordsF[0][iFlx].push_back(ksi_value);
		CFreal eta_value = start_Node1[1]+alpha[iFlx]*uFace1[1];
		coordsF[0][iFlx].push_back(eta_value);
	}
	
	// Face 2
	for (CFuint iFlx=0; iFlx<solOrder+1; ++iFlx){
		CFreal ksi_value = start_Node2[0]+alpha[iFlx]*uFace2[0];
		coordsF[1][iFlx].push_back(ksi_value);
		CFreal eta_value = start_Node2[1]+alpha[iFlx]*uFace2[1];
		coordsF[1][iFlx].push_back(eta_value);
	}
	// Face 3
	for (CFuint iFlx=0; iFlx<solOrder+1; ++iFlx){
		CFreal ksi_value = start_Node3[0]+alpha[iFlx]*uFace3[0];
		coordsF[2][iFlx].push_back(ksi_value);
		CFreal eta_value = start_Node3[1]+alpha[iFlx]*uFace3[1];
		coordsF[2][iFlx].push_back(eta_value);
	}

  return coordsF;
}
//////////////////////////////////////////////////////////////////////
void TriagFluxReconstructionElementData::createFaceNormals()
{
  CFAUTOTRACE;
  
  // number of faces
  const CFuint nbrFaces = 3;
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    // number of face nodes
    const CFuint nbrFaceNodes = m_faceNodeConn[iFace].size();

    m_faceNodeCoords[iFace].resize(nbrFaceNodes);
    for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
    {
      m_faceNodeCoords[iFace][iNode].resize(m_dimensionality);

      // node ID
      const CFuint nodeID = m_faceNodeConn[iFace][iNode];

      // set coordinate
      m_faceNodeCoords[iFace][iNode] = m_cellNodeCoords[nodeID];
    }
  }
  // compute the normals
  m_faceNormals.resize(nbrFaces);
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    // resize the variable
    m_faceNormals[iFace].resize(m_dimensionality);
    // compute normal
    m_faceNormals[iFace][KSI] =  (m_faceNodeCoords[iFace][1][ETA] - m_faceNodeCoords[iFace][1][ETA])*0.5;
    m_faceNormals[iFace][ETA] = -(m_faceNodeCoords[iFace][0][KSI] - m_faceNodeCoords[iFace][0][KSI])*0.5;
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceNodeConnectivityPerOrient()
{
  CFAUTOTRACE;

  // number of faces
  const CFuint nbrFaces = m_faceNodeConn.size();

  // number of possible orientations
  const CFuint nbrOrient = 6;
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
  cf_assert(iOrient == nbrOrient);
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createSolPolyExponents()
{
  CFAUTOTRACE;

  // define exponents
  m_solPolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < static_cast< CFuint >(m_polyOrder)+1; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < m_polyOrder+1-iKsi; ++iEta)
    {
      vector< CFint > solPolyExps(2);
      solPolyExps[KSI] = iKsi;
      solPolyExps[ETA] = iEta;
      m_solPolyExponents.push_back(solPolyExps);
    }
  }
//////////////
   /* // helper variable
  const CFuint polyOrderP1 = m_polyOrder + 1;

  // number of polynomial terms
  const CFuint nbrPolyTerms = (m_polyOrder + 1)*(m_polyOrder + 2)/2;

  // resize the variable
  m_solPolyExponents.resize(nbrPolyTerms);
  for (CFuint iTerm = 0; iTerm < nbrPolyTerms; ++iTerm)
  {
    m_solPolyExponents[iTerm].resize(2);
  }

  // define exponents
  CFuint iTerm = 0;
  for (CFuint iP = 0; iP < polyOrderP1; ++iP)
  {
    for (CFuint iY = 0; iY < iP+1; ++iY, ++iTerm)
    {
      m_solPolyExponents[iTerm][KSI] = iP-iY;
      m_solPolyExponents[iTerm][ETA] = iY;
    }
  }*/

  cf_assert(m_solPolyExponents.size() == static_cast<CFuint>((m_polyOrder+1)*(m_polyOrder+2)/2));
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createNodePolyExponents()
{
  CFAUTOTRACE;

  // define exponents
  m_nodePolyExponents.resize(0);
  for (CFuint iKsi = 0; iKsi < 2; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < 2-iKsi; ++iEta)
    {
      vector< CFint > nodePolyExps(2);
      nodePolyExps[KSI] = iKsi;
      nodePolyExps[ETA] = iEta;
      m_nodePolyExponents.push_back(nodePolyExps);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceFluxPntsConn() //Modif here
{
 CFAUTOTRACE;

  // number of flux points in 2D per Face
  const CFuint nbrFlxPnts2D = m_polyOrder+1;//flxPntsLocalCoord2D[0].size();

  // resize m_faceFlxPntConn
  m_faceFlxPntConn.resize(3);
  CFuint iFlx1 = 0; 
  // variable holding the face index
for(CFuint faceIdx= 0; faceIdx<3; ++faceIdx){
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx, ++iFlx1)
  {
    m_faceFlxPntConn[faceIdx].push_back(iFlx1); 
  }
}
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFluxPolyExponents()
{
  CFAUTOTRACE;
  // number of orientations
  const CFuint nbrOrients = 6;

  // number of flux points in 2D per Face
  const CFuint nbrFlxPnts2D = m_polyOrder+1;

  // create data structure
  m_faceFlxPntConnPerOrient.resize(nbrOrients);
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < 3; ++iFaceL)
  {
    for (CFuint iFaceR = iFaceL; iFaceR < 3; ++iFaceR, ++iOrient)
    {
      m_faceFlxPntConnPerOrient[iOrient].resize(2);
      for (CFuint iSol = 0; iSol < nbrFlxPnts2D; ++iSol)
      {
        m_faceFlxPntConnPerOrient[iOrient][LEFT ]
            .push_back(m_faceFlxPntConn[iFaceL][iSol  ]);
        m_faceFlxPntConnPerOrient[iOrient][RIGHT]
            .push_back(m_faceFlxPntConn[iFaceR][nbrFlxPnts2D-1-iSol]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::setInterpolationNodeSet(const CFPolyOrder::Type order,
                                                         vector< RealVector >& nodalSet)
{
  
  if(solPntsLocalCoord2D.size()!=0){
    const CFuint nbrSolPoint2D = (order+1)*(order+2)/2 ;

    // set solution point local coordinates
    nodalSet.resize(0);
    for (CFuint iSol = 0; iSol < nbrSolPoint2D; ++iSol)
    {
      RealVector node(2);
      node[KSI] = solPntsLocalCoord2D[iSol][0];
      node[ETA] = solPntsLocalCoord2D[iSol][1];
      nodalSet.push_back(node);
    }
    //cout<<" setInterpolationNodeSet --  end  "<<  endl;
  }
  else{
    //cout<<" Not good !!!!!!!!"<< endl;
  }

}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::setCFLConvDiffRatio()
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
      throw Common::NotImplementedException (FromHere(),"Higher-order triag FR cell not defined!");
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceOutputPntCellMappedCoords()
{
  const CFuint nbrFacePnts = m_polyOrder == 0 ? 1 :m_polyOrder;

  // face mapped coordinates of uniform distribution of points
  const CFreal dKsi = m_polyOrder == 0 ? 1.0 : 1.0/m_polyOrder;
  CFreal ksi = 0.0;
  m_faceOutputPntFaceMappedCoords.resize(0);
  for (CFuint iPnt = 0; iPnt < nbrFacePnts; ++iPnt, ksi += dKsi)
  {
    RealVector mapCoord(1);
    mapCoord[KSI] = ksi;
    m_faceOutputPntFaceMappedCoords.push_back(mapCoord);
  }

  // compute cell mapped coordinates for distribution on each face
  const CFuint nbrCellFaces = 3;
  m_faceOutputPntCellMappedCoords.resize(nbrCellFaces);
  for (CFuint iFace = 0; iFace < nbrCellFaces; ++iFace)
  {
    // current face node coordinates
    const vector<RealVector>& faceNodeCoords = m_faceNodeCoords[iFace];
    m_faceOutputPntCellMappedCoords[iFace].resize(0);
    for (CFuint iPnt = 0; iPnt < nbrFacePnts; ++iPnt)
    {
      const CFreal fun0 = 0.5*(1.0-m_faceOutputPntFaceMappedCoords[iPnt][KSI]);
      const CFreal fun1 = 0.5*(1.0+m_faceOutputPntFaceMappedCoords[iPnt][KSI]);
      m_faceOutputPntCellMappedCoords[iFace].push_back(fun0*faceNodeCoords[0]+fun1*faceNodeCoords[1]);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceOutputPntConn()
{
  CFAUTOTRACE;
  m_faceOutputPntConn.resize(0);

  for (CFuint iCell = 0; iCell < (CFuint) m_polyOrder; ++iCell)
  {
    vector<CFuint> cellNode(2);
    cellNode[0] = iCell;
    cellNode[1] = iCell+1;

    m_faceOutputPntConn.push_back(cellNode);
  }
}

//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFlxPntsLocalCoords()
{
  // number of flux points in 1D

  const CFuint nbrFlxPnts2D = m_polyOrder+1;

  // set flux point local coordinates
  m_flxPntsLocalCoords.resize(0);
  // loop over 1D Flx pnts for each face
  for (CFuint iFace=0; iFace <3; ++iFace){

    for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx)
    {
      RealVector flxCoords(2);
      flxCoords[KSI] = flxPntsLocalCoord2D[iFace][iFlx][0];
      flxCoords[ETA] = flxPntsLocalCoord2D[iFace][iFlx][1];
      m_flxPntsLocalCoords.push_back(flxCoords);
    }
  }
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createSolPntsLocalCoords()
{   
  // number of solution points in 1D
  const CFuint nbrSolPnts2D = (m_polyOrder+1)*(m_polyOrder+2)/2 ;

  // set solution point local coordinates
  m_solPntsLocalCoords.resize(0);
  for (CFuint iSol = 0; iSol< nbrSolPnts2D; ++iSol)
  {
    RealVector solCoords(2);
    solCoords[KSI] = solPntsLocalCoord2D[iSol][0]; 
    solCoords[ETA] = solPntsLocalCoord2D[iSol][1];
    m_solPntsLocalCoords.push_back(solCoords);
  }
  cf_assert(m_solPntsLocalCoords.size() == nbrSolPnts2D);
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createFaceFlxPntsFaceLocalCoords() 
{
  CFAUTOTRACE;

  // number of flux points in 1D
  const CFuint nbrFlxPnts1D = (m_polyOrder+1);
  CFuint nbfaces = 3;
  // set face flux point face local coordinates
  m_faceFlxPntsFaceLocalCoords.resize(0);
  for (CFuint iFlx=0; iFlx<nbrFlxPnts1D; ++iFlx)
  {  
    RealVector flxCoord(1);
    flxCoord[KSI] = ((m_flxPntsLocalCoord1D[iFlx]));
    m_faceFlxPntsFaceLocalCoords.push_back(flxCoord);
  }
}


//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createFaceFluxPntsConnPerOrient()
{
  CFAUTOTRACE;
  // number of orientations
  const CFuint nbrOrients = 6;

  // number of flux points in 2D per Face
  const CFuint nbrFlxPnts2D = m_polyOrder+1; //flxPntsLocalCoord2D.size();

  // create data structure
  m_faceFlxPntConnPerOrient.resize(nbrOrients);
  CFuint iOrient = 0;
  for (CFuint iFaceL = 0; iFaceL < 3; ++iFaceL)
  {
    for (CFuint iFaceR = iFaceL; iFaceR < 3; ++iFaceR, ++iOrient)
    {
      m_faceFlxPntConnPerOrient[iOrient].resize(2);
      for (CFuint iSol = 0; iSol < nbrFlxPnts2D; ++iSol)
      {
        m_faceFlxPntConnPerOrient[iOrient][LEFT ]
            .push_back(m_faceFlxPntConn[iFaceL][iSol]);
        m_faceFlxPntConnPerOrient[iOrient][RIGHT]
            .push_back(m_faceFlxPntConn[iFaceR][nbrFlxPnts2D-1-iSol]);
      }
    }
  }

}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createCellNodeCoords()
{
  CFAUTOTRACE;

  m_cellNodeCoords.resize(3);

  // first node
  m_cellNodeCoords[0].resize(2);
  m_cellNodeCoords[0][KSI] = 0.;
  m_cellNodeCoords[0][ETA] = 0.;

  // second node
  m_cellNodeCoords[1].resize(2);
  m_cellNodeCoords[1][KSI] = 1.; 
  m_cellNodeCoords[1][ETA] = 0.;

  // third node
  m_cellNodeCoords[2].resize(2);
  m_cellNodeCoords[2][KSI] = 0.;
  m_cellNodeCoords[2][ETA] = 1.;
}

//////////////////////////////////////////////////////////////////////
  
void TriagFluxReconstructionElementData::createFaceNodeConnectivity()
{
  CFAUTOTRACE;

  m_faceNodeConn.resize(3);

  m_faceNodeConn[0].resize(2);
  m_faceNodeConn[0][0] = 0;
  m_faceNodeConn[0][1] = 1;

  m_faceNodeConn[1].resize(2);
  m_faceNodeConn[1][0] = 1;
  m_faceNodeConn[1][1] = 2;

  m_faceNodeConn[2].resize(2);
  m_faceNodeConn[2][0] = 2;
  m_faceNodeConn[2][1] = 0;

}

//////////////////////////////////////////////////////////////////////
void TriagFluxReconstructionElementData::createFaceMappedCoordDir()
{
  CFAUTOTRACE;

  m_faceMappedCoordDir.resize(3);

  m_faceMappedCoordDir[0] = 2.;
  m_faceMappedCoordDir[1] = 2.;
  m_faceMappedCoordDir[2] = 2.; 

}
//////////////////////////////////////////////////////////////////////
void TriagFluxReconstructionElementData::createFlxSolDependencies()
{
  CFAUTOTRACE;
  
  const CFuint nbrSolPnts = solPntsLocalCoord2D.size();
  const CFuint nbrFlxPnts2D = 3*(m_polyOrder+1); //flxPntsLocalCoord2D[0].size();

  m_solFlxDep.resize(nbrSolPnts);
  m_solSolDep.resize(nbrSolPnts);  // also push back the current 
  m_flxSolDep.resize(nbrFlxPnts2D);

  CFuint iSol = 0;
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts2D; ++iFlx){
      for(CFuint iSol= 0; iSol< nbrSolPnts; ++iSol){
        m_flxSolDep[iFlx].push_back(iSol);
      }
    }
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol){
      for(CFuint iFlx= 0; iFlx< nbrFlxPnts2D; ++iFlx){
        m_solFlxDep[iSol].push_back(iFlx);
      }
    }
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol){
      for(CFuint iSol1= 0; iSol1< nbrSolPnts; ++iSol1){
        m_solSolDep[iSol].push_back(iSol1);
      }
    }
}
//////////////////////////////////////////////////////////////////////
void TriagFluxReconstructionElementData::createVandermondeMatrix()
{
  CFAUTOTRACE;
  
  const CFuint nbrSolPnts = m_solPntsLocalCoords.size();
  // number of solution points in 1D
  const CFuint nbrSolPnts1D = m_solPntsLocalCoord1D.size();
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
            double a = ((2.*m_solPntsLocalCoords[iSol][KSI])/(1.-m_solPntsLocalCoords[iSol][ETA]))-1.;
            double b = 2.*m_solPntsLocalCoords[iSol][ETA]-1.;
            double h1 = ComputeJacobi(iOrderKsi, 0., 0., a);
            double h2 = ComputeJacobi(iOrderEta, 2.*iOrderKsi+1., 0., b);
            m_vandermonde(iSol,modalDof)=sqrt(2.0)*h1*h2*pow((1.-b),iOrderKsi);
            //m_vandermonde(iSol,modalDof) = evaluateLegendre(-1.+2.*m_solPntsLocalCoords[iSol][KSI],iOrderKsi)*evaluateLegendre(-1.+2.*m_solPntsLocalCoords[iSol][ETA],iOrderEta);
            modalDof+=1;
          } 
        }
    }
    InvertMatrix(m_vandermonde,m_vandermondeInv);
  }
}

  
//////////////////////////////////////////////////////////////////////

void TriagFluxReconstructionElementData::createFaceIntegrationCoefs()
{
  CFAUTOTRACE;

  // number of flux points on a face
  const CFuint nbrFlxPnts = m_flxPntsLocalCoord1D.size();

  // resize m_faceIntegrationCoefs
  m_faceIntegrationCoefs.resize(nbrFlxPnts);

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
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_faceIntegrationCoefs[iFlx] = 0.0;

    const CFreal ksiFlx = m_flxPntsLocalCoord1D[iFlx];
    for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
    {
      // quadrature point local coordinate on the face
      const CFreal ksiQPnt = quadPntCoords[iQPnt][KSI];

      // evaluate polynomial value in quadrature point
      CFreal quadPntPolyVal = 1.;
      for (CFuint iFac = 0; iFac < nbrFlxPnts; ++iFac)
      {
        if (iFac != iFlx)
        {
          const CFreal ksiFac = m_flxPntsLocalCoord1D[iFac];
          quadPntPolyVal *= (ksiQPnt-ksiFac)/(ksiFlx-ksiFac);
        }
      }

      // add contribution of quadrature point to integration coefficient
      m_faceIntegrationCoefs[iFlx] += quadPntWheights[iQPnt]*quadPntPolyVal;
    }
  }
}

//////////////////////////////////////////////////////////////////////
void TriagFluxReconstructionElementData::createFluxPntFluxDim()
{
  CFAUTOTRACE;

  m_flxPntFlxDim.resize(3*(m_polyOrder+1));
  
    for (CFuint iFlx = 0; iFlx < (m_polyOrder+1); ++iFlx)
    { 
        m_flxPntFlxDim[m_faceFlxPntConn[0][iFlx]] = 0;//1; //along y
        m_flxPntFlxDim[m_faceFlxPntConn[1][iFlx]] = 1;//2; //oblique
        m_flxPntFlxDim[m_faceFlxPntConn[2][iFlx]] = 2;//0; //along x
    }
}
//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createCellAvgSolCoefs()
  {  // to be completed

    m_cellAvgSolCoefs.resize((m_polyOrder+1)*(m_polyOrder+2)/2);

    switch(m_polyOrder)
      {
        case CFPolyOrder::ORDER0:
        {
          m_cellAvgSolCoefs[0] = 1.0;
        } break;
        case CFPolyOrder::ORDER1:
        {
          m_cellAvgSolCoefs[0] = 0.333333333333333;
          m_cellAvgSolCoefs[1] = 0.333333333333333;
          m_cellAvgSolCoefs[2] = 0.333333333333333;

        } break;
        case CFPolyOrder::ORDER2:
        {
          m_cellAvgSolCoefs[0] = 0.109951743655333;
          m_cellAvgSolCoefs[1] = 0.109951743655333;
          m_cellAvgSolCoefs[2] = 0.109951743655333;
          m_cellAvgSolCoefs[3] = 0.223381589678000;
          m_cellAvgSolCoefs[4] = 0.223381589678000;
          m_cellAvgSolCoefs[5] = 0.223381589678000;

        } break;
        case CFPolyOrder::ORDER3:
        {
          m_cellAvgSolCoefs[0] = 0.041955512996649;
          m_cellAvgSolCoefs[1] = 0.041955512996649;
          m_cellAvgSolCoefs[2] = 0.041955512996649;
          m_cellAvgSolCoefs[3] = 0.112098412070887;
          m_cellAvgSolCoefs[4] = 0.112098412070887;
          m_cellAvgSolCoefs[5] = 0.112098412070887;
          m_cellAvgSolCoefs[6] = 0.112098412070887;
          m_cellAvgSolCoefs[7] = 0.112098412070887;
          m_cellAvgSolCoefs[8] = 0.112098412070887;
          m_cellAvgSolCoefs[9] = 0.201542988584730;
          
        } break;
        case CFPolyOrder::ORDER4:
        {
          m_cellAvgSolCoefs[0] = 0.017915455012303;
          m_cellAvgSolCoefs[1] = 0.017915455012303;
          m_cellAvgSolCoefs[2] = 0.017915455012303;
          m_cellAvgSolCoefs[3] = 0.127712195881265;
          m_cellAvgSolCoefs[4] = 0.127712195881265;
          m_cellAvgSolCoefs[5] = 0.127712195881265;
          m_cellAvgSolCoefs[6] = 0.076206062385535;
          m_cellAvgSolCoefs[7] = 0.076206062385535;
          m_cellAvgSolCoefs[8] = 0.076206062385535;
          m_cellAvgSolCoefs[9] = 0.055749810027115;
          m_cellAvgSolCoefs[10] = 0.055749810027115;
          m_cellAvgSolCoefs[11] = 0.055749810027115;
          m_cellAvgSolCoefs[12] = 0.055749810027115;
          m_cellAvgSolCoefs[13] = 0.055749810027115;
          m_cellAvgSolCoefs[14] = 0.055749810027115;
          
        } break;
        case CFPolyOrder::ORDER5:
        {
          m_cellAvgSolCoefs[0] = 0.010359374696538;
          m_cellAvgSolCoefs[1] = 0.010359374696538;
          m_cellAvgSolCoefs[2] = 0.010359374696538;
          m_cellAvgSolCoefs[3] = 0.075394884326738;
          m_cellAvgSolCoefs[4] = 0.075394884326738;
          m_cellAvgSolCoefs[5] = 0.075394884326738;
          m_cellAvgSolCoefs[6] = 0.097547802373242;
          m_cellAvgSolCoefs[7] = 0.097547802373242;
          m_cellAvgSolCoefs[8] = 0.097547802373242;
          m_cellAvgSolCoefs[9] = 0.028969269372473;
          m_cellAvgSolCoefs[10] = 0.028969269372473;
          m_cellAvgSolCoefs[11] = 0.028969269372473;
          m_cellAvgSolCoefs[12] = 0.028969269372473;
          m_cellAvgSolCoefs[13] = 0.028969269372473;
          m_cellAvgSolCoefs[14] = 0.028969269372473;
          m_cellAvgSolCoefs[15] = 0.046046366595935;
          m_cellAvgSolCoefs[16] = 0.046046366595935;
          m_cellAvgSolCoefs[17] = 0.046046366595935;
          m_cellAvgSolCoefs[18] = 0.046046366595935;
          m_cellAvgSolCoefs[19] = 0.046046366595935;
          m_cellAvgSolCoefs[20] = 0.046046366595935;

        } break;
        case CFPolyOrder::ORDER6:
        {
          m_cellAvgSolCoefs[0] = 0.083608212215637;

          m_cellAvgSolCoefs[1] = 0.005272170280495;
          m_cellAvgSolCoefs[2] = 0.005272170280495;
          m_cellAvgSolCoefs[3] = 0.005272170280495;

          m_cellAvgSolCoefs[4] = 0.044552936679504;
          m_cellAvgSolCoefs[5] = 0.044552936679504;
          m_cellAvgSolCoefs[6] = 0.044552936679504;

          m_cellAvgSolCoefs[7] = 0.033815712804198;
          m_cellAvgSolCoefs[8] = 0.033815712804198;
          m_cellAvgSolCoefs[9] = 0.033815712804198;

          m_cellAvgSolCoefs[10] = 0.015710461340183;
          m_cellAvgSolCoefs[11] = 0.015710461340183;
          m_cellAvgSolCoefs[12] = 0.015710461340183;
          m_cellAvgSolCoefs[13] = 0.015710461340183;
          m_cellAvgSolCoefs[14] = 0.015710461340183;
          m_cellAvgSolCoefs[15] = 0.015710461340183;

          m_cellAvgSolCoefs[16] = 0.028205136280616;
          m_cellAvgSolCoefs[17] = 0.028205136280616;
          m_cellAvgSolCoefs[18] = 0.028205136280616;
          m_cellAvgSolCoefs[19] = 0.028205136280616;
          m_cellAvgSolCoefs[20] = 0.028205136280616;
          m_cellAvgSolCoefs[21] = 0.028205136280616;


          m_cellAvgSolCoefs[22] = 0.066995957127830;
          m_cellAvgSolCoefs[23] = 0.066995957127830;
          m_cellAvgSolCoefs[24] = 0.066995957127830;
          m_cellAvgSolCoefs[25] = 0.066995957127830;
          m_cellAvgSolCoefs[26] = 0.066995957127830;
          m_cellAvgSolCoefs[27] = 0.066995957127830;

        } break;
        case CFPolyOrder::ORDER7:
        {

        } break;
        case CFPolyOrder::ORDER8:
        {

        } break;
        case CFPolyOrder::ORDER9:
        {

        } break;
        case CFPolyOrder::ORDER10:
        {

        } break;
      }
    };

//////////////////////////////////////////////////////////////////////
  
  void TriagFluxReconstructionElementData::createCellCenterDerivCoefs(){
      CFAUTOTRACE;

  // center coordinate
  vector< RealVector > centerCoord(1,RealVector(2));
  centerCoord[0][KSI] = 0.33;
  centerCoord[0][ETA] = 0.33;

  vector< vector< vector< CFreal > > > polyDerivs =
      getSolPolyDerivsAtNode(centerCoord);

  // number of solution points
  const CFuint nbrSolPnts = getNbrOfSolPnts();

  // set polynomial derivatives
  m_cellCenterDerivCoefs.resize(2);
  m_cellCenterDerivCoefs[KSI].resize(nbrSolPnts);
  m_cellCenterDerivCoefs[ETA].resize(nbrSolPnts);
  for (CFuint iPoly = 0; iPoly < nbrSolPnts; ++iPoly)
  {
    m_cellCenterDerivCoefs[KSI][iPoly] = polyDerivs[0][KSI][iPoly];
    m_cellCenterDerivCoefs[ETA][iPoly] = polyDerivs[0][ETA][iPoly];
  }
};

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

