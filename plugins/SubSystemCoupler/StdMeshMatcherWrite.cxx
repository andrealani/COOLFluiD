#include <boost/progress.hpp>

#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MeshData.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "SubSystemCoupler/StdMeshMatcherWrite.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdMeshMatcherWrite, SubSysCouplerData, SubSystemCouplerModule> StdMeshMatcherWriteProvider("StdMeshMatcherWrite");

//////////////////////////////////////////////////////////////////////////////

StdMeshMatcherWrite::StdMeshMatcherWrite(const std::string& name) :
  CouplerCom(name),
  _sockets(),
  _matchingFace(static_cast<Framework::TopologicalRegionSet*>(CFNULL),CFNULL),
  _shapeFunctionAtCoord()
{
}

//////////////////////////////////////////////////////////////////////////////

StdMeshMatcherWrite::~StdMeshMatcherWrite()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdMeshMatcherWrite::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite::configure ( Config::ConfigArgs& args )
{
  CouplerCom::configure(args);

  typedef std::vector<Common::SafePtr<CommandGroup> > InterfaceList;
  InterfaceList interfaces = getMethodData().getInterfaces();

  try {

  InterfaceList::iterator itr = interfaces.begin();
  for(; itr != interfaces.end(); ++itr) {

    const std::vector<std::string>& comNames = (*itr)->getComNames();

    // check if this command applies to this Interface
    if(count(comNames.begin(),comNames.end(),getName()) != 0) {

      const std::string interfaceName = (*itr)->getName();

      for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(); ++iProc)
      {
        const vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);
        for (CFuint iTRS=0; iTRS< otherTrsNames.size(); ++iTRS)
        {
          // Create the datahandle for the otherSubSystem Coordinates
          const std::string otherTrsName = otherTrsNames[iTRS];

          // Get the datahandle for the coordinates of the OTHER SubSystem
          vector<std::string> socketCoordNames = getMethodData().getOtherCoupledCoordName(interfaceName,otherTrsName,iProc);
          vector<std::string> socketAcceptNames = getMethodData().getOtherCoupledAcceptedName(interfaceName,otherTrsName,iProc);
          vector<std::string> socketDataNames = getMethodData().getOtherCoupledDataName(interfaceName,otherTrsName,iProc);
          for(CFuint iType=0;iType < socketCoordNames.size();iType++)
          {
            _sockets.createSocketSink<RealVector>(socketCoordNames[iType]);
            _sockets.createSocketSink<CFreal>(socketAcceptNames[iType]);
            _sockets.createSocketSink<RealVector>(socketDataNames[iType]);
          }
        } // other trs
      } //loop over processors
    } // if check
  } // interfaces

  } catch (Exception& e)
  {
    CFout << e.what() << "\n" << CFendl;
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite::execute()
{
  CFAUTOTRACE;

  for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(); ++iProc) {
    executeWrite(iProc);
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite::executeWrite(const CFuint iProc)
{
  CFAUTOTRACE;

  // Get the names of the interfaces, subsystems
  const std::string interfaceName = getCommandGroupName();
  vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);

  // If the geometry is non matching, then get the threshold for acceptance
  // else set the threshold to a very high value
  const bool isNonMatchingGeometry = getMethodData().isInterfaceGeometryNonMatching(interfaceName);
  CFreal nonMatchingGeometryThreshold = MathTools::MathConsts::CFrealMax();
  if(isNonMatchingGeometry) nonMatchingGeometryThreshold = getMethodData().getNonMatchingGeometryThreshold(interfaceName);

  // Loop over the TRSs of the Interface
  for (CFuint iTRS=0; iTRS < otherTrsNames.size(); iTRS++)
  {
    CFLogInfo ("Matching for TRS [" << otherTrsNames[iTRS] << "]\n");

    const vector<std::string> socketCoordNames =
      getMethodData().getOtherCoupledCoordName(interfaceName,otherTrsNames[iTRS],iProc);
    const vector<std::string> socketAcceptNames =
      getMethodData().getOtherCoupledAcceptedName(interfaceName,otherTrsNames[iTRS],iProc);
    const vector<std::string> socketDataNames =
      getMethodData().getOtherCoupledDataName(interfaceName,otherTrsNames[iTRS],iProc);

    for(CFuint iType=0; iType < socketCoordNames.size();iType++)
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
      for (CFuint iState = 0; iState < otherNbStates; ++iState, ++progress) {

        CFint nodeID = -1;
        CFuint idx = iState-rejectedStates;
        bool accepted = false;
        _minimumDistanceOnFace = MathTools::MathConsts::CFrealMax();
        _minimumDistanceOffFace = MathTools::MathConsts::CFrealMax();

        const RealVector coord = interfaceCoords[iState];
        RealVector coordProj(coord.size());

        // Pair each fluid point with the closest wet structural element
        nodeToElementPairing(coord,nodeID,coordProj);

        // if the projection of the point is outside all faces, take the closest node
        if (nodeID != -1){
          if(nonMatchingGeometryThreshold > _minimumDistanceOffFace) accepted = true;

          if(accepted)
          {
            isAccepted[iState] = _minimumDistanceOffFace;
            // Resize
            (*coupledGeoEntities)[idx].third.resize(_shapeFunctionAtCoord.size());
            (*coupledGeoEntities)[idx].fourth.resize(PhysicalModelStack::getActive()->getDim());
            // set values
            (*coupledGeoEntities)[idx].first = _matchingFace.first;
            (*coupledGeoEntities)[idx].second = _matchingFace.second;
            (*coupledGeoEntities)[idx].third = 0.;
            (*coupledGeoEntities)[idx].third[nodeID] = 1.;
            (*coupledGeoEntities)[idx].fourth = coordProj;
            CFLogDebugMin("Node: "<< coord << " does not project on a face and was projected on coord: " << coordProj << "\n");
          }
          else
          {
            isAccepted[iState] = -1.;
            rejectedStates++;
            CFout << "Node doesn't project on a face and is rejected because distance is: "<< _minimumDistanceOffFace << "\n";
          }
        }
        // else continue
        else {
          if(nonMatchingGeometryThreshold > _minimumDistanceOnFace) accepted = true;

          if(accepted)
          {
            isAccepted[iState] = _minimumDistanceOnFace;
            // Determine the natural coordinates of the projection of the fluid point
            // and compute the shape function at this location
            // Resize
            (*coupledGeoEntities)[idx].third.resize(_shapeFunctionAtCoord.size());
            (*coupledGeoEntities)[idx].fourth.resize(PhysicalModelStack::getActive()->getDim());
            // set values
            (*coupledGeoEntities)[idx].first = _matchingFace.first;
            (*coupledGeoEntities)[idx].second = _matchingFace.second;
            (*coupledGeoEntities)[idx].third = _shapeFunctionAtCoord;
            (*coupledGeoEntities)[idx].fourth = coordProj;
            CFLogDebugMin ("Node: "<< coord << " projects on a face at coord: " << coordProj << "\n");
          }
          else
          {
            isAccepted[iState] = -1.;
            rejectedStates++;
            CFout << "Node projects on a face but is rejected because distance is: " << _minimumDistanceOnFace << "\n";
          }
        } //end else
      } //end loop over otherStates

      ///We know the number of accepted states -> can now resize the datahandle
      interfaceData.resize(otherNbStates - rejectedStates);

      ///Write the file with the info isAccepted
      writeIsAcceptedFile(socketAcceptNames[iType]);
    } //end of loop over data transfer coord type
  } // end loop over the OtherTRS
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite::nodeToElementPairing(const RealVector& coord, CFint& nodeID, RealVector& coord_Proj)
{
  CFAUTOTRACE;

  const CFuint dim = coord.size();
  RealVector normal(dim);
  RealVector x0(dim);
  RealVector x1(dim);
  RealVector v(dim);
  RealVector v1(dim);
  RealVector w(dim);
  RealVector project(dim);
  RealVector project2(dim);
  RealVector projectCoord(dim);
  CFreal tempV, tempW;
  bool isOnFace(false);

  /// Loop over all the faces of all the TRS's of this command
  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  vector< SafePtr<TopologicalRegionSet> >::iterator iTRS;

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  for (iTRS = trs.begin(); iTRS != trs.end(); ++iTRS) {

      const CFuint nbGeos = (*iTRS)->getLocalNbGeoEnts();
      const std::string currentTRSName = (*iTRS)->getName();
      geoData.trs = (*iTRS);

      for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {

        // build the GeometricEntity
        geoData.idx = iGeoEnt;
        GeometricEntity& currFace = *geoBuilder->buildGE();

        /// Check if node defined by coord has his projection in face...
        // 1 - compute the normal to the face
        // compute the face normal in 2D
        cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_2D);
        ///@todo modify this for 3D

        x0 = *(currFace.getNode(0));
        x1 = *(currFace.getNode(1));
        v1 = x1-x0;

        normal = currFace.computeAvgCellNormal();
//         normal[0] = -v1[1];
//         normal[1] = v1[0];

        // 2 - project vector xP-x0 (=v) on the normal from x0 -> w = (n*v)*v
        v = coord - x0;
        w = (normal * v);

        // 3 - Obtain the vector xProj-x0 (= v - w)
        project = v-w ;
        projectCoord = project + x0 ;
        project2 = projectCoord - x1;
        // 4 - Now we have the point projected on the plane of the face
        //     We have to check if the point is inside the face
        //  We check that (xProj -x0) < (x1-x0)

        if ((project.norm2() < v1.norm2()) && (project2.norm2() < v1.norm2()))
          {
          isOnFace = true;
          // 5 - if it falls inside the face, then OK...
          //     but with concave surfaces, there might be more than one...
          //     so continue and select the face for which the distance
          //     to ??? is minimum!!
          //     Distance to be minimized: averaged distance from the points
          //                               minimum distance between xP and any node
          //                               ...

          // Compute the distance of projected point to the nodes and take minimum
          v = projectCoord - x0;
          w = projectCoord - x1;
          tempV = sqrt(v[0]*v[0] + v[1]*v[1]);
          tempW = sqrt(w[0]*w[0] + w[1]*w[1]);
          if ((tempV < _minimumDistanceOnFace) || (tempW < _minimumDistanceOnFace)) {
            coord_Proj = projectCoord;
            _matchingFace.first = (*iTRS);
            _matchingFace.second = iGeoEnt;
            _shapeFunctionAtCoord.resize(currFace.nbNodes());
            _shapeFunctionAtCoord = currFace.computeShapeFunctionAtCoord(coord_Proj);
/*CFout << "Face iGeoEnt: " << iGeoEnt << "\n";
CFout << "Node 0: " << x0 << "\n";
CFout << "Node 1: " << x1 << "\n";
CFout << "Coord Proj: " << coord_Proj << "\n";
CFout << "_shapeFunctionAtCoord: " << _shapeFunctionAtCoord << "\n";*/
            if(tempV < tempW)
            {
              _minimumDistanceOnFace = tempV;
            }
            else
            {
              _minimumDistanceOnFace = tempW;
            }
            //once is has been projected on a face, no need of this
            _minimumDistanceOffFace = -1.;

            }
          }
        else
          {
            if(!isOnFace)
            {
              // 6 - the point projection might fall outside all faces
              // Compute the coordinates of projected point
              coord_Proj = project + x0 ;

              v = coord_Proj - x0;
              w = coord_Proj - x1;
              tempV = sqrt(v[0]*v[0] + v[1]*v[1]);
              tempW = sqrt(w[0]*w[0] + w[1]*w[1]);
              _shapeFunctionAtCoord.resize(currFace.nbNodes());
              _shapeFunctionAtCoord = 0.;
              if (tempV < _minimumDistanceOffFace)
              {
                nodeID = 0;
                _matchingFace.first = (*iTRS);
                _matchingFace.second = iGeoEnt;
                _shapeFunctionAtCoord[nodeID] = 1.;
                _minimumDistanceOffFace = tempV;
              }
              if (tempW < _minimumDistanceOffFace)
              {
                nodeID = 1;
                _matchingFace.first = (*iTRS);
                _matchingFace.second = iGeoEnt;
                _shapeFunctionAtCoord[nodeID] = 1.;
                _minimumDistanceOffFace = tempW;
              }
            }
          }

    //release the GeometricEntity
    geoBuilder->releaseGE();

    } // end of loop over faces
  } // end of loop over TRS's

  if(isOnFace == true){
    nodeID = -1;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite::writeIsAcceptedFile(const std::string dataHandleName)
{
  CFAUTOTRACE;

// CFout << "Writing the isAccepted file: " << dataHandleName << "\n";

  DataHandle< CFreal > isAccepted = _sockets.getSocketSink<CFreal>(dataHandleName)->getDataHandle();

  boost::filesystem::path dpath (dataHandleName);
  boost::filesystem::path nameOutputFile =
    Environment::DirPaths::getInstance().getResultsDir() / dpath;

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(nameOutputFile);

  const CFuint nbStates = isAccepted.size();

  fout << nbStates << "\n";
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    fout << isAccepted[iState] << "\n";
  }

  fhandle->close();
// CFout << "Done Writing the isAccepted file: " << dataHandleName << "\n";
}


//////////////////////////////////////////////////////////////////////////////


    } // namespace MeshMatcher

  } // namespace Numerics

} // namespace COOLFluiD
