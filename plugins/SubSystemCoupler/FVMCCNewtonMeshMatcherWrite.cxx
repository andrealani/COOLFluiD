#include <boost/progress.hpp>

#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Environment/DirPaths.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "SubSystemCoupler/FVMCCNewtonMeshMatcherWrite.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCCNewtonMeshMatcherWrite, SubSysCouplerData, SubSystemCouplerModule> FVMCCNewtonMeshMatcherWriteProvider("FVMCCNewtonMeshMatcherWrite");


//////////////////////////////////////////////////////////////////////////////

FVMCCNewtonMeshMatcherWrite::FVMCCNewtonMeshMatcherWrite(const std::string& name) :
  NewtonMeshMatcherWrite(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCCNewtonMeshMatcherWrite::~FVMCCNewtonMeshMatcherWrite()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCNewtonMeshMatcherWrite::executeWrite(const CFuint iProc)
{

  CFAUTOTRACE;

  // Names of the interfaces, subsystems
  const std::string interfaceName = getCommandGroupName();
  vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);
  SubSysCouplerData::GeoEntityIdx matchingFace;

  // If the geometry is non matching, then get the threshold for acceptance
  // else set the threshold to a very high value
  const bool isNonMatchingGeometry = getMethodData().isInterfaceGeometryNonMatching(interfaceName);
  CFreal nonMatchingGeometryThreshold = MathTools::MathConsts::CFrealMax();
  if(isNonMatchingGeometry) nonMatchingGeometryThreshold = getMethodData().getNonMatchingGeometryThreshold(interfaceName);

  // Get the geometric entity builder
  Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  for (CFuint iTRS=0; iTRS < otherTrsNames.size(); iTRS++)
  {
    CFLogInfo("Matching for TRS [" << otherTrsNames[iTRS] << "]\n");
    const vector<std::string> socketCoordNames =
      getMethodData().getOtherCoupledCoordName(interfaceName,otherTrsNames[iTRS],iProc);
    const vector<std::string> socketAcceptNames =
      getMethodData().getOtherCoupledAcceptedName(interfaceName,otherTrsNames[iTRS],iProc);
    const vector<std::string> socketDataNames =
      getMethodData().getOtherCoupledDataName(interfaceName,otherTrsNames[iTRS],iProc);

    for(CFuint iType=0; iType < socketCoordNames.size(); iType++)
    {
      CFLogInfo("Matching [" << socketCoordNames[iType] << "]\n");

      DataHandle< RealVector> interfaceCoords =
        _sockets.getSocketSink<RealVector>(socketCoordNames[iType])->getDataHandle();
      DataHandle< CFreal> isAccepted =
        _sockets.getSocketSink<CFreal>(socketAcceptNames[iType])->getDataHandle();
      DataHandle< RealVector> interfaceData =
        _sockets.getSocketSink<RealVector>(socketDataNames[iType])->getDataHandle();


      // Get the coupledGeoEntities and resize
      const CFuint otherNbStates = interfaceCoords.size();
      SubSysCouplerData::CoupledGeoEntities* coupledGeoEntities =
        getMethodData().getCoupledInterfaces(interfaceName,iTRS, iType, iProc);

      (*coupledGeoEntities).resize(otherNbStates);

      // Counter for the number of rejected states
      CFuint rejectedStates = 0;

      // For each state, find the interpolated value
      // on the other boundary
      boost::progress_display progress (otherNbStates);
      for (CFuint i = 0; i < otherNbStates; ++i, ++progress)
      {
        CFuint idx = i-rejectedStates;
        bool accepted = false;
        _minimumDistanceOnFace = MathTools::MathConsts::CFrealMax();

        _coord.resize(interfaceCoords[i].size());
        _coord = interfaceCoords[i];

        // PreSelect a given number of faces whose centroid is close to the point to project
        facePreSelection();

        // Pair each fluid point with the closest wet structural element
        nodeToElementPairing(matchingFace);

        // Acceptation of the projection depending on
        // the distance between point and its projection
        if(nonMatchingGeometryThreshold > _minimumDistanceOnFace) accepted = true;

        if(accepted)
        {
          CFLogDebugMin("Accepted: " <<_coord << "\n");
          cf_assert(matchingFace.first.isNotNull());

          isAccepted[i] = _minimumDistanceOnFace;
          geoData.trs = matchingFace.first;
          geoData.idx = matchingFace.second;
          geoData.isBFace = true;
          GeometricEntity& matchingFaceGeo = *geoBuilder->buildGE();

          (*coupledGeoEntities)[idx].first = matchingFace.first;
          (*coupledGeoEntities)[idx].second = matchingFace.second;

          // Compute the shape function at the projected point
          CFuint nbNodesInFace = matchingFaceGeo.getNodes()->size();
          (*coupledGeoEntities)[idx].third.resize(nbNodesInFace);
          (*coupledGeoEntities)[idx].fourth.resize(PhysicalModelStack::getActive()->getDim());

          (*coupledGeoEntities)[idx].third = matchingFaceGeo.computeGeoShapeFunctionAtMappedCoord(_mappedCoordMin);
          _tempCoord = matchingFaceGeo.computeCoordFromMappedCoord(_mappedCoordMin);
          (*coupledGeoEntities)[idx].fourth = _tempCoord;

          CFLogDebugMin("...and associated to: " << _tempCoord << "\n");
          CFLogDebugMin("...idx: " << idx << "\n");
          CFLogDebugMin("...Distance: " << _minimumDistanceOnFace << "\n");
          CFLogDebugMin("...ShapeFunction: " << (*coupledGeoEntities)[idx].third << "\n");
          std::vector<State*>* states = matchingFaceGeo.getStates();
          CFLogDebugMin("...nbStates: " << states->size() << "\n");
          CFLogDebugMin("...State[0]: " << (*states)[0]->getCoordinates() << "\n");
          CFLogDebugMin("...State[1]: " << (*states)[1]->getCoordinates() << "\n");
          if(states->size() > 2) { CFLogDebugMin("...State[2]: " << (*states)[2]->getCoordinates() << "\n"); }
          if(states->size() > 3) { CFLogDebugMin("...State[3]: " << (*states)[3]->getCoordinates() << "\n"); }

          //release the geometric entity
          geoBuilder->releaseGE();
        }
        else
        {
//           CFout << " =================== Rejected: " <<_coord << "\n";
          isAccepted[i] = -1.;
          rejectedStates++;
        }

      } //end loop over otherStates

      ///We know the number of accepted states -> can now resize the datahandle
      interfaceData.resize(otherNbStates - rejectedStates);

      ///Write the file with the info isAccepted
      writeIsAcceptedFile(socketAcceptNames[iType]);
    } //end of loop over data transfer coord type
  } // end loop over the OtherTRS
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCNewtonMeshMatcherWrite::facePreSelection()
{
  CFAUTOTRACE;

  RealVector centroid(_coord.size());
  CFreal maxDistanceInSelection = MathTools::MathConsts::CFrealMax();
  CFreal distance;
  CFuint iMax = 0;

  ///Initialize _faceSelection
  _faceSelection.resize(_nbSelectedFaces);
  for(CFuint iFace=0;iFace < _nbSelectedFaces;iFace++)
  {
    _faceSelection[iFace].first = CFNULL;
  }

  /// Loop over all the faces of all the TRS's of this command
  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  vector< SafePtr<TopologicalRegionSet> >::iterator iTRS;

  Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  for (iTRS = trs.begin(); iTRS != trs.end(); ++iTRS)
  {
    const CFuint nbGeos = (*iTRS)->getLocalNbGeoEnts();
    const std::string currentTRSName = (*iTRS)->getName();
    geoData.trs = (*iTRS);
    geoData.isBFace = true;

    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
    {
      // build the GeometricEntity
      geoData.idx = iGeoEnt;
      _currentFace = geoBuilder->buildGE();

      // Compute distance from point to centroid of the face
      centroid = _currentFace->computeCentroid();

      distance = 0.;
      for(CFuint iDim = 0;iDim < _coord.size();iDim++)
      {
        distance += (centroid[iDim]-_coord[iDim])*(centroid[iDim]-_coord[iDim]);
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();

      // If distance is smaller than the maximum of the preselected distance
      // then unselect the previous maximum and select this face
      // Then, recompute wich face is the further away...
      if (distance < maxDistanceInSelection)
      {
        //Set the new value at the position of the previous face
        _faceSelection[iMax].first = (*iTRS);
        _faceSelection[iMax].second = iGeoEnt;

        //Recompute the furthest face from the point
        maxDistanceInSelection = 0.;
        bool hasNull = false;
        for(CFuint iFace=0;iFace < _nbSelectedFaces;iFace++)
        {
          if(_faceSelection[iFace].first.isNotNull())
          {
            geoData.trs = _faceSelection[iFace].first;
            geoData.idx = _faceSelection[iFace].second;
            GeometricEntity& iFaceGeo = *geoBuilder->buildGE();

            centroid = iFaceGeo.computeCentroid();
            distance = 0.;
            for(CFuint iDim = 0;iDim < _coord.size();iDim++)
            {
              distance += (centroid[iDim]-_coord[iDim])*(centroid[iDim]-_coord[iDim]);
            }

            if((distance>maxDistanceInSelection) && (!hasNull))
            {
              iMax = iFace;
              maxDistanceInSelection = distance;
            }

            // reset the trs to the current trs
            geoData.trs = (*iTRS);
            // release the geometric entity
            geoBuilder->releaseGE();
          }
          else
          {
            hasNull = true;
            iMax = iFace;
            maxDistanceInSelection = MathTools::MathConsts::CFrealMax();
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCNewtonMeshMatcherWrite::nodeToElementPairing(SubSysCouplerData::GeoEntityIdx& matchingFace)
{

  Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  /// Loop over all the preSelected faces
  for (CFuint iFace = 0; iFace< _faceSelection.size(); ++iFace)
  {
    if(_faceSelection[iFace].first.isNotNull())
    {

      geoData.trs = _faceSelection[iFace].first;
      geoData.idx = _faceSelection[iFace].second;
      geoData.isBFace = true;
      _currentFace = &(*geoBuilder->buildGE());
  /*std::cout << "=====================================" << std::endl;
  std::cout << "Node 0: " << *(_currentFace->getNode(0)) << std::endl;
  std::cout << "Node 1: " << *(_currentFace->getNode(1)) << std::endl;
  std::cout << "Node 2: " << *(_currentFace->getNode(2)) << std::endl;
  std::cout << "=====================================" << std::endl;*/
      /// Check if node defined by coord has his projection in face...
      computeProjection();

      if (_bestDistance < _minimumDistanceOnFace)
      {
        matchingFace.first = geoData.trs;
        matchingFace.second = geoData.idx;

        _mappedCoordMin.resize(_bestCoord.size());
        _mappedCoordMin = _bestCoord;
        _minimumDistanceOnFace = _bestDistance;
      }
      // release the geometric entity
      geoBuilder->releaseGE();
    }
  } // end of loop over faces
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshMatcher

  } // namespace Numerics

} // namespace COOLFluiD
