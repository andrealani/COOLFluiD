#include "MathTools/MathChecks.hh"

#include "Framework/MethodCommandProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/PeriodicNonMatching.hh"
#include <math.h>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PeriodicNonMatching, CellCenterFVMData, FiniteVolumeModule> periodicNonMatchingFVMCCProvider("PeriodicNonMatchingFVMCC");

//////////////////////////////////////////////////////////////////////////////

void PeriodicNonMatching::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
}
      
//////////////////////////////////////////////////////////////////////////////

PeriodicNonMatching::PeriodicNonMatching(const std::string& name) :
  FVMCC_BC(name)
{    
  addConfigOptionsTo(this);
  
  _threshold = 1e-4;
  setParameter("Threshold",&_threshold);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicNonMatching::~PeriodicNonMatching()
{
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicNonMatching::setup()
{
  
/// creation of the vector<PeriodicFace> of the boundaries (bottom and top)
  FVMCC_BC::setup();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
	      geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  
  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;
  geoData.trs = trs;

  const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();

  // mapping the global face IDs to local IDs in the TRS and (it is not necessary!)

  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
     CFLogDebugMed( "iFace = " << iFace << "\n");
     // build the GeometricEntity
     geoData.idx = iFace;
     GeometricEntity *const face = geoBuilder->buildGE();
     const CFuint faceGlobalID = face->getID();
     _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);

     // release the GeometricEntity
     geoBuilder->releaseGE();
  }

  _globalToLocalTRSFaceID.sortKeys();
  //const CFuint nbTrFaces = nbTrsFaces/2;

  
  vector<int> vectorIndex;
  for(int r=0; r<nbTrsFaces; r++){
    vectorIndex.push_back(r);
  } //end for
 
  int k = 0;
  
  for(int s=0;s<2;s++){

    geoData.idx = k;
    GeometricEntity *const face = geoBuilder->buildGE();
    const vector<Node*>& nodes = *face->getNodes();  // std::vector<Node*>* getNodes()  see Framework/GeometricEntity
    PeriodicFace PF;
    PF.setupFirstNode(*(nodes[0]));
    PF.setupSecondNode(*(nodes[1]));
    std::vector<PeriodicFace> PeriodicFaceA;
    PeriodicFaceA.push_back(PF); 
    geoBuilder->releaseGE(); 
    Node nodeFirst(PF.getFirstNode());
    Node nodeLast(PF.getSecondNode());
  
    for(int i=0;i<nbTrsFaces;i++){
   
      int j=vectorIndex[i];
      geoData.idx = j;
      GeometricEntity *const face = geoBuilder->buildGE();
      const vector<Node*>& nodes = *face->getNodes(); 
      PeriodicFace PF_S;
      PF_S.setupFirstNode(*(nodes[0]));
      PF_S.setupSecondNode(*(nodes[1]));
      geoBuilder->releaseGE();
  
      if(MathTools::MathChecks::isEqualWithError(nodeFirst[XX],PF_S.getFirstNode()[XX], _threshold) &&     MathTools::MathChecks::isEqualWithError(nodeFirst[YY],PF_S.getFirstNode()[YY], _threshold)){
   
        Node temp(PF_S.getFirstNode());
        PF_S.setupFirstNode(PF_S.getSecondNode());
        PF_S.setupSecondNode(temp);
        nodeFirst=PF_S.getFirstNode();
        vectorIndex.erase(i);
        PeriodicFaceA.insert(PeriodicFaceBottom.begin(),PF_S);
      } //end if

      if(MathTools::MathChecks::isEqualWithError(nodeLast[XX],PF_S.getFirstNode()[XX], _threshold) &&     MathTools::MathChecks::isEqualWithError(nodeLast[YY],PF_S.getFirstNode()[YY], _threshold)){
  
        nodeLast=PF_S.getSecondNode();
        vectorIndex.erase(i);
        PeriodicFaceA.push_back(PF_S);
      } //end if

      if(MathTools::MathChecks::isEqualWithError(nodeFirst[XX],PF_S.getSecondNode()[XX], _threshold) &&     MathTools::MathChecks::isEqualWithError(nodeFirst[YY],PF_S.getSecondNode()[YY], _threshold)){
        
        nodeFirst=PF_S.getFirstNode();
        vectorIndex.erase(i);
        PeriodicFaceA.insert(PeriodicFaceBottom.begin(),PF_S);
      } //end if

      if(MathTools::MathChecks::isEqualWithError(nodeLast[XX],PF_S.getSecondNode()[XX], _threshold) &&     MathTools::MathChecks::isEqualWithError(nodeLast[YY],PF_S.getSecondNode()[YY], _threshold)){
        
        Node temp(PF_S.getFirstNode());
        PF_S.setupFirstNode(PF_S.getSecondNode());
        PF_S.setupSecondNode(temp);
        nodeLast=PF_S.getSecondNode();
        vectorIndex.erase(i);
        PeriodicFaceA.push_back(PF_S);
      } //end if

    } //end for
 
    if(s==0){
  
      PeriodicFaceBottom(PeriodicFaceA);
    }
    else{
    
    PeriodicFaceTop(PeriodicFaceA);
    }
    k=vectorIndex[0];
  
  } //end for

 /// translation and rotation of the bottom boundary towards the top boundary
  CFreal X_bottom_first = PeriodicFaceBottom[0].getFirstNode()[XX];
  CFreal Y_bottom_first = PeriodicFaceBottom[0].getFirstNode()[YY];
  CFreal X_bottom_last = PeriodicFaceBottom[PeriodicFaceBottom.end()].getSecondNode()[XX];
  CFreal Y_bottom_last = PeriodicFaceBottom[PeriodicFaceBottom.end()].getSecondNode()[YY];
  CFreal X_top_first = PeriodicFaceTop[0].getFirstNode()[XX];
  CFreal Y_top_first = PeriodicFaceTop[0].getFirstNode()[YY];
  CFreal X_top_last = PeriodicFaceTop[PeriodicFaceTop.end()].getSecondNode()[XX];
  CFreal Y_top_last = PeriodicFaceTop[PeriodicFaceTop.end()].getSecondNode()[YY];
  CFreal gap_X = X_bottom_first - X_top_first;
  CFreal gap_y = Y_bottom_first - Y_top_first;
  CFreal alfa_top = atan((Y_top_last-Y_top_first)/(X_top_last-X_top_first));
  CFreal alfa_bottom = atan((Y_bottom_last-Y_bottom_first)/(X_bottom_last-X_bottom_first));
  CFreal alfa = alfa_bottom-alfa_top;

  CFuint nFacesTop = PeriodicFaceTop.size();
  CFreal point_X_first = PeriodicFaceBottom[0].getFirstNode()[XX] ;
  CFreal point_Y_first = PeriodicFaceBottom[0].getFirstNode()[YY];
  Node& point(point_X,point_Y);
  PeriodicFaceTop[0].setFirstNode(point);
  for(int l=1;l<nFacesBottom;l++){
    CFreal point_X = (PeriodicFaceTop[0].getFirstNode()[XX] + gap_X) - point_X_first;
    CFreal point_Y = (PeriodicFaceTop[0].getFirstNode()[YY] + gap_Y) - point_Y_first;
    CFreal module = pow(pow(point_X,2) + pow(point_Y,2),0.5);
    CFreal inclination_old = atan(point_Y/point_X);
    CFreal inclination_new = inclination_old - alfa;
    if(inclination_new>=0 && inclination_new<pi/2){
      CFreal sign_X = 1; 
      CFreal sign_Y = 1;
    }
    if(inclination_new>=pi/2 && inclination_new<pi){
      CFreal sign_X = -1; 
      CFreal sign_Y = 1;
    }
    if(inclination_new<=0 && inclination_new<-pi/2){
      CFreal sign_X = 1; 
      CFreal sign_Y = -1;
    }
    if(inclination_new<=-pi/2 && inclination_new<-pi){
      CFreal sign_X = -1; 
      CFreal sign_Y = -1;
    }
    CFreal point_X_new = sign_X*(module/(pow((1+pow(tan(inclination_new),2)),0.5));
    CFreal point_Y_new = sign_Y*(module/(pow((1+pow(tan(inclination_new),2)),0.5))*tan(inclination_new);
    Node& point(point_X_new,point_Y_new);
    PeriodicFaceTop[l].setFirstNode(point);
  } //end for

/// construction of the array PeriodicFaceNew
  nFaceBottom = PeriodicFaceBottom.size();
  nFaceTop = PeriodicFaceTop.size();

  int indexBottom = 0;
  int indexTop = 0;
  int indexNew =0;

  while(indexBottom != nFaceBottom && indexTop != nFaceTop)
  {

    PeriodicFace faceBottom = PeriodicFaceBottom[indexBottom];
    PeriodicFace faceTop = PeriodicFaceTop[indexTop];
    PeriodicFace faceNew = PeriodicFaceNew[k];

    // definition of the point asked for the calculations
    CFreal X_B_2 = faceBottom.getSecondNode()[XX];
    CFreal X_T_2 = faceTop.getSecondNode()[XX];
    CFreal X_N_1 = faceNew.getFirstNode()[XX];
    CFreal Y_B_2 = faceBottom.getSecondNode()[YY];
    CFreal Y_T_2 = faceTop.getSecondNode()[YY];
    CFreal Y_N_1 = faceNew.getFirstNode()[YY];

    // calculate point BT
    CFreal delta = 0.5*pow(pow(Y_T_2-Y_B_2,2)+pow(X_T_2-X_T_2,2),0.5);
    CFreal alfa = atan((Y_T_2-Y_B_2)/(X_T_2-X_B_2));
    CFreal X_BT = X_T_2 - delta*cos(alfa);
    CFreal Y_BT = Y_T_2 + delta*sin(alfa);

    // calculate poit BP
    CFreal deltaBP = ((X_B_2-X_N_1)*(X_BT-X_N_1)+(Y_B_2-Y_N_1)*(Y_BT-Y_N_1))/(pow(pow(Y_BT-Y_N_1,2)+pow(X_BT-X_N_1,2),0.5));
    CFreal alfaBP = atan((Y_BT-Y_N_1)/(X_BT-X_N_1));
    CFreal X_BP = X_N_1 - deltaBP*cos(alfaBP);
    CFreal Y_BP = Y_N_1 + deltaBP*sin(alfaBP);

    // calculate poit TP
    CFreal deltaTP = ((X_T_2-X_N_1)*(X_BT-X_N_1)+(Y_T_2-Y_N_1)*(Y_BT-Y_N_1))/(pow(pow(Y_BT-Y_N_1,2)+pow(X_BT-X_N_1,2),0.5));
    //CFreal alfaTP = atan((Y_BT-Y_N_1)/(X_BT-X_N_1));  alfaTP=alfaBP
    CFreal X_TP = X_N_1 - deltaTP*cos(alfaBP);
    CFreal Y_TP = Y_N_1 + deltaTP*sin(alfaBP);

    if(MathTools::MathChecks::isGreaterWithError(deltaBP,deltaTP, _threshold))  
    {
    Node nodeNew;
    nodeNew[XX]=X_TP;
    nodeNew[YY]=Y_TP;    
    faceNew.setSecondNode(nodeNew);
    PeriodicFaceNew.push_back(faceNew);
    // if you want you can erase PeriodicFaceTop[indexTop]
    indexTop++;
    indexNew++;
    }
    
    if(MathTools::MathChecks::isSmallerWithError(deltaBP,deltaTP, _threshold)) 
    {
    Node nodeNew;
    nodeNew[XX]=X_BP;
    nodeNew[YY]=Y_BP;
    faceNew.setSecondNode(nodeNew);
    PeriodicFaceNew.push_back(faceNew);
    // if you want you can erase PeriodicFaceBottom[indexBottom]
    indexBottom++;
    indexNew++;
    }

    if(MathTools::MathChecks::isEqualWithError(deltaBP,deltaTP, _threshold))
    {
    Node nodeNew;
    nodeNew[XX]=X_TP;
    nodeNew[YY]=Y_TP;
    faceNew.setSecondNode(nodeNew);
    PeriodicFaceNew.push_back(faceNew);
    // if you want you can erase PeriodicFaceTop[indexTop]
    indexTop++;
    indexBottom++;
    indexNew++;
    }

} //end while

/*
*
* to do: is the "setGhostState" function to do?
*
*/
void PeriodicNonMatching::setGhostState(GeometricEntity *const face)
{
   
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD





