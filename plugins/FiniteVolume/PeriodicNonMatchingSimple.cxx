#include "MathTools/MathChecks.hh"

#include "Framework/MethodCommandProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/PeriodicNonMatchingSimple.hh"
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

MethodCommandProvider<PeriodicNonMatchingSimple, CellCenterFVMData, FiniteVolumeModule> periodicNonMatchingSimpleFVMCCProvider("PeriodicNonMatchingSimpleFVMCC");

//////////////////////////////////////////////////////////////////////////////

void PeriodicNonMatchingSimple::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
}
      
//////////////////////////////////////////////////////////////////////////////

PeriodicNonMatchingSimple::PeriodicNonMatchingSimple(const std::string& name) :
  FVMCC_BC(name)
{    
  addConfigOptionsTo(this);
  
  _threshold = 1e-4;
  setParameter("Threshold",&_threshold);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicNonMatchingSimple::~PeriodicNonMatchingSimple()
{
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicNonMatchingSimple::setup()
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
  
// to do: create and sort of the vector PeriodicFaceBottom PeriodicFaceTop 

/// construction of the array PeriodicFaceNew
  nFaceBottom = PeriodicFaceBottom.size();
  nFaceTop = PeriodicFaceTop.size();

  int indexBottom = 0;
  int indexTop = 0;
  
  while(indexBottom != nFaceBottom && indexTop != nFaceTop)
  {

    PeriodicFace faceBottom = PeriodicFaceBottom[indexBottom];
    PeriodicFace faceTop = PeriodicFaceTop[indexTop];
    
    // definition of the point asked for the calculations
    CFreal X_B_2 = faceBottom.getSecondNode()[XX]; // +deltaX
    CFreal X_T_2 = faceTop.getSecondNode()[XX];
    CFreal Y_B_2 = faceBottom.getSecondNode()[YY];
    CFreal Y_T_2 = faceTop.getSecondNode()[YY];
    CFreal X_B_1 = faceBottom.getFirstNode()[XX]; // +deltaX
    CFreal X_T_1 = faceTop.getFirstNode()[XX];
    CFreal Y_B_1 = faceBottom.getFirstNode()[YY];
    CFreal Y_T_1 = faceTop.getFirstNode()[YY];
   
    // calculate length of bottom to top proiection
    CFreal TopLength =(pow(pow(Y_T_2-Y_T_1,2)+pow(X_T_2-X_T_1,2),0.5);
    CFreal deltaBtoT = ((X_T_2-X_T_1)*(X_B_2-X_T_1)+(Y_T_2-Y_T_1)*(Y_B_2-Y_T_1))/TopLength;

    // calculate length of top to bottom proiection
    CFreal BottomLength = (pow(pow(Y_B_2-Y_B_1,2)+pow(X_B_2-X_B_1,2),0.5))
    CFreal deltaTtoB = ((X_B_2-X_B_1)*(X_T_2-X_B_1)+(Y_B_2-Y_B_1)*(Y_T_2-Y_B_1))/BottomLength;

    if(MathTools::MathChecks::isSmallerWithError(deltaBtoT,TopLength, _threshold))  
    {
    alfa = atan((Y_T_2-Y_T_1)/(X_T_2-X_T_1));
    Node nodeNew;
    nodeNew[XX]=X_T_1 + deltaBtoT*cos(alfa);
    nodeNew[YY]=Y_T_1 + deltaBtoT*sin(alfa); 
    // do other things   
    indexBottom++;
    }
    
    if(MathTools::MathChecks::isSmallerWithError(deltaTtoB,BottomLength, _threshold))  
    {
    alfa = atan((Y_B_2-Y_B_1)/(X_B_2-X_B_1));
    Node nodeNew;
    nodeNew[XX]=X_B_1 + deltaTtoB*cos(alfa);
    nodeNew[YY]=Y_B_1 + deltaTtoB*sin(alfa); 
    // do other things   
    indexTop++;
    }

    if(MathTools::MathChecks::isEqualWithError(deltaBtoT,TopLength, _threshold))
    {
    // do other things
    indexTop++;
    indexBottom++;
    }
  } //end while

    
/*
*
* to do: is the "setGhostState" function to do?
*
*/
void PeriodicNonMatchingSimple::setGhostState(GeometricEntity *const face)
{
   
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD



