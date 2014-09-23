// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "MeshTools/MeshToolsFVM.hh"
#include "MeshTools/ComputeWallDistanceVector2CC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeWallDistanceVector2CC, DataProcessingData, MeshToolsFVMModule> ComputeWallDistanceVector2CCProvider("ComputeWallDistanceVector2CC");

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVector2CC::ComputeWallDistanceVector2CC(const std::string& name) :
  ComputeWallDistance(name),
  socket_gstates("gstates"),
  socket_normals("normals")
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVector2CC::~ComputeWallDistanceVector2CC()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVector2CC::setup()
{
  CFAUTOTRACE;

  ComputeWallDistance::setup();

   const CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
   m_tempNode.resize(dim);
   m_midNode.resize(dim);
   m_tempGhostNode.resize(dim);
   m_faceNormal.resize(dim);


}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeWallDistanceVector2CC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = ComputeWallDistance::needsSockets();

  result.push_back(&socket_gstates);
  result.push_back(&socket_normals);

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVector2CC::execute()
{
  CFAUTOTRACE;

  CFLog(INFO, "Computing wall distances using Vector's method...\n");
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  const CFuint nbStates = states.size();
  CFreal minimumDistance;
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    minimumDistance = MathTools::MathConsts::CFrealMax();
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size() ; ++iTRS)
    {
      Common::SafePtr<Framework::TopologicalRegionSet> const faces =
         MeshDataStack::getActive()->getTrs(_boundaryTRS[iTRS]);
      Common::SafePtr<std::vector<CFuint> > nodesInTrs = faces->getNodesInTrs();

      Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
          geoBuilder = getMethodData().getFaceTrsGeoBuilder();

      SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
      geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

      FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
      geoData.isBFace = true;
      geoData.trs = faces;

      const CFuint nbFaces = faces->getLocalNbGeoEnts();

      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
        CFLogDebugMed( "Computing iFace = " << iFace << "\n");

        // build the GeometricEntity
        geoData.idx = iFace;

        GeometricEntity& currFace = *geoBuilder->buildGE();

        const CFuint dim = PhysicalModelStack::getActive()->getDim();
        const CFuint faceID = currFace.getID();
        const CFuint startID = faceID*dim;

        // set the current normal
        for (CFuint i = 0; i < dim; ++i) {
          m_faceNormal[i] = normals[startID + i];
        }

        if(dim == 2){
          // compute the original position of the ghost state @see ComputeDummyState
          const Node& firstNode = *currFace.getNode(0);
          const Node& secondNode = *currFace.getNode(1);

          RealVector faceVector = secondNode - firstNode;
          RealVector nodeStateVector = states[iState]->getCoordinates() - firstNode;
          CFreal project = MathFunctions::innerProd(faceVector, nodeStateVector)/faceVector.norm2();

          if(project < faceVector.norm2())
          {
            faceVector.normalize();
            RealVector projectCoord = firstNode + project * faceVector;
            CFreal stateFaceDistance = MathFunctions::getDistance(projectCoord, states[iState]->getCoordinates());
            if(stateFaceDistance < minimumDistance) minimumDistance = stateFaceDistance;
          }
        }

        if(dim == 3){
          // compute the original position of the ghost state @see ComputeDummyState
          const Node& firstNode = *currFace.getNode(0);
	  
	      // a condition has to be added to check that the projection is internal to the face
          RealVector nodeStateVector = states[iState]->getCoordinates() - firstNode;
          CFreal stateFaceDistance = MathFunctions::innerProd(m_faceNormal, nodeStateVector)/m_faceNormal.norm2();
          if(stateFaceDistance < minimumDistance) minimumDistance = stateFaceDistance;
        }
        geoBuilder->releaseGE();
      }
    }

    wallDistance[iState] = minimumDistance;
  }

  printToFile();
  CFLog(INFO, "Wall distances computation finished...\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
