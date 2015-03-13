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

  CFLog(VERBOSE, "ComputeWallDistanceVector2CC::execute() START\n");

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector nodeStateVector(dim);
  RealVector faceCentroid(dim);
  
  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
	geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;
  
  const CFuint nbStates = states.size();
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    CFreal minimumDistance = MathTools::MathConsts::CFrealMax();
    CFreal minStateFaceDistance =  MathTools::MathConsts::CFrealMax();
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size() ; ++iTRS)
    {
      geoData.trs = MeshDataStack::getActive()->getTrs(_boundaryTRS[iTRS]);
      
      const CFuint nbFaces = geoData.trs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
        CFLogDebugMed( "Computing iFace = " << iFace << "\n");
	
        // build the GeometricEntity
        geoData.idx = iFace;
        GeometricEntity& currFace = *geoBuilder->buildGE();
	
	faceCentroid = 0.;
	const CFuint nbNodesInFace = currFace.nbNodes();
	for (CFuint n = 0; n < nbNodesInFace; ++n) {
	  faceCentroid += *currFace.getNode(n);
	}
	const CFreal ovNbNodesInFace = 1./(CFreal)nbNodesInFace;
	faceCentroid *= ovNbNodesInFace;
	
	const Node& firstNode = *currFace.getNode(0);
	nodeStateVector = states[iState]->getCoordinates() - firstNode;
	const CFreal stateFaceDistance = MathFunctions::getDistance(states[iState]->getCoordinates(), faceCentroid);
	if (stateFaceDistance < minStateFaceDistance) {
	  for (CFuint i = 0; i < dim; ++i) {
            const CFuint faceID = currFace.getID();
            const CFuint startID = faceID*dim;
	    m_faceNormal[i] = normals[startID + i];
	  }
	  // normal to boundary face is always pointing outward with respect to the computational domain 
	  // a "-" sign needs to be considered
	  minimumDistance = -MathFunctions::innerProd(m_faceNormal, nodeStateVector)/m_faceNormal.norm2();
	  minStateFaceDistance = stateFaceDistance;
	}
	
	geoBuilder->releaseGE();
      }
    }
    
    wallDistance[iState] = std::abs(minimumDistance);
  }

  printToFile();
  
  CFLog(VERBOSE, "ComputeWallDistanceVector2CC::execute() END\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
