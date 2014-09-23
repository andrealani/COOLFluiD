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
#include "MeshTools/ComputeWallDistanceVectorCC.hh"

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

MethodCommandProvider<ComputeWallDistanceVectorCC, DataProcessingData, MeshToolsFVMModule> ComputeWallDistanceVectorCCProvider("ComputeWallDistanceVectorCC");

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVectorCC::ComputeWallDistanceVectorCC(const std::string& name) :
  ComputeWallDistance(name),
  socket_gstates("gstates"),
  socket_normals("normals")
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVectorCC::~ComputeWallDistanceVectorCC()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVectorCC::setup()
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
ComputeWallDistanceVectorCC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = ComputeWallDistance::needsSockets();

  result.push_back(&socket_gstates);
  result.push_back(&socket_normals);

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVectorCC::execute()
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

        // compute the original position of the ghost state @see ComputeDummyState
        const Node& firstNode = *currFace.getNode(0);
        const CFreal k = - MathFunctions::innerProd(m_faceNormal, firstNode);
        const CFreal n2 = MathFunctions::innerProd(m_faceNormal, m_faceNormal);
        cf_assert(std::abs(n2) > 0.0);
        State *const innerState = currFace.getState(0);
        m_innerNode = &innerState->getCoordinates();
        m_t = (MathFunctions::innerProd(m_faceNormal,*m_innerNode) + k)/n2;
        m_tempGhostNode = (*m_innerNode) - 2.*m_t*m_faceNormal;

        // this middle node is by construction on the boundary face
        m_midNode = 0.5*(*m_innerNode + m_tempGhostNode);

        //Compute the distance between the state and the middle node
        CFreal stateFaceDistance = MathFunctions::getDistance(m_midNode, states[iState]->getCoordinates());

        if(stateFaceDistance < minimumDistance) minimumDistance = stateFaceDistance;

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
