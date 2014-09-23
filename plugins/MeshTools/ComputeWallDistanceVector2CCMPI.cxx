// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"

#ifdef CF_HAVE_MPI
#include "Common/MPI/MPIStructDef.hh"
#endif

#include "Common/Stopwatch.hh"

#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/TRSDistributeData.hh"

#include "MeshTools/MeshToolsFVM.hh"
#include "MeshTools/ComputeWallDistanceVector2CCMPI.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeWallDistanceVector2CCMPI, DataProcessingData, MeshToolsFVMModule>
ComputeWallDistanceVector2CCMPIProvider("ComputeWallDistanceVector2CCMPI");

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVector2CCMPI::ComputeWallDistanceVector2CCMPI(const std::string& name) :
  ComputeWallDistance(name),
#ifdef CF_HAVE_MPI
  m_comm(MPI_COMM_WORLD),
#endif
  m_myRank(0),
  m_nbProc(1)
{
}
    
//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVector2CCMPI::~ComputeWallDistanceVector2CCMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVector2CCMPI::setup()
{
  CFAUTOTRACE;

  ComputeWallDistance::setup();
  
#ifdef CF_HAVE_MPI
  m_comm   = PE::GetPE().GetCommunicator();
  m_myRank = PE::GetPE().GetRank();
  m_nbProc = PE::GetPE().GetProcessorCount();
#endif
  
  // initialize the wall distance
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
  wallDistance = MathTools::MathConsts::CFrealMax();
}

//////////////////////////////////////////////////////////////////////////////

vector<Common::SafePtr<BaseDataSocketSink> >
ComputeWallDistanceVector2CCMPI::needsSockets()
{
  return ComputeWallDistance::needsSockets();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVector2CCMPI::execute()
{
  CFAUTOTRACE;
  
  // this has only to be run at setup for now
  CFLog(VERBOSE, "ComputeWallDistanceVector2CCMPI::execute() START\n");
  
  Stopwatch<WallTime> stp;
  stp.start();
  
  CFLog(INFO,"ComputeWallDistanceVector2CCMPI::execute() computing distance to the wall ...\n");
  
  //  (PhysicalModelStack::getActive()->getDim() == DIM_3D) ? execute3D() : execute2D();
  execute3D();
  
  CFLog(INFO,"ComputeWallDistanceVector2CCMPI::execute() took " << stp.read() << "s\n");
  CFLog(VERBOSE, "ComputeWallDistanceVector2CCMPI::execute() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVector2CCMPI::execute3D()
{
  CFAUTOTRACE;
    
  DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_always_assert(_boundaryTRS.size() > 0);
  
  RealVector faceCentroid(dim);
  vector<CFreal> trsFaceData;
  
  // loop over all processors
  for (CFuint root = 0; root < m_nbProc; ++root) {
    
    // loop over all wall boundary TRSs
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size(); ++iTRS) {  
      CFLog(VERBOSE, "ComputeWallDistanceVector2CCMPI::execute() Processing TRS named " << _boundaryTRS[iTRS] << "\n");
      
      SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs(_boundaryTRS[iTRS]);
      
      // broadcast number of faces in the current TRS to allow for buffer memory preallocation
      CFuint nbLocalTrsFaces = faces->getLocalNbGeoEnts();
      
#ifdef CF_HAVE_MPI
      MPI_Bcast(&nbLocalTrsFaces, 1, MPI_UNSIGNED, root, m_comm);
#endif
      
      // something is gonna be packed to be sent only if the number of faces 
      // inside the current TRS in the root processor is positive
      if (nbLocalTrsFaces > 0) {
	// trsFaceData stores the coordinates of the wall face centroid
	const CFuint sizeData = nbLocalTrsFaces*dim;
	trsFaceData.resize(sizeData); 
	
	if (m_myRank == root) {
	  // loop over faces in the current TRS and store connectivity info into TRSDistributeData
	  CFuint count = 0;
	  for (CFuint iFace = 0; iFace < nbLocalTrsFaces; ++iFace) {
	    // compute the face centroid
	    faceCentroid = 0.;
	    const CFuint nbNodesInFace = faces->getNbNodesInGeo(iFace);
	    for (CFuint n = 0; n < nbNodesInFace; ++n) {
	      const CFuint nodeID = faces->getNodeID(iFace, n);
	      faceCentroid += *nodes[nodeID];
	    }
	    const CFreal ovNbNodesInFace = 1./(CFreal)nbNodesInFace;
	    faceCentroid *= ovNbNodesInFace;
	    
	    for (CFuint iDim = 0; iDim < dim; ++iDim, ++count) {
	      trsFaceData[count] = faceCentroid[iDim];
	    }
	  }
	  cf_assert(count == trsFaceData.size());
	}
	
#ifdef CF_HAVE_MPI
	MPI_Bcast(&trsFaceData[0], sizeData, MPI_DOUBLE, root, m_comm);
#endif
	
	// compute the distance to the wall starting from broadcast data 
	computeWallDistance3D(trsFaceData);
      }
    }
  }
  
  // printToFile();
}
    
//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVector2CCMPI::computeWallDistance3D(std::vector<CFreal>& data)
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbStates = states.size();
  RealVector faceCentroid(dim, (CFreal*)CFNULL);
  const CFuint nbFaces = data.size()/dim;
  
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    cf_assert(iState == states[iState]->getLocalID());
    const RealVector& stateCoord = states[iState]->getCoordinates(); 
    
    // initialize the minimum distance to a huge value
    CFreal minimumDistance = MathTools::MathConsts::CFrealMax();
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFreal start = iFace*dim;
      cf_assert(start < data.size());
      faceCentroid.wrap(dim, &data[start]);
      
      // compute distance between face centroid and cell centroid
      const CFreal stateFaceDistance = MathFunctions::getDistance(faceCentroid, stateCoord);
      minimumDistance = std::min(stateFaceDistance, minimumDistance);
    }
    
    wallDistance[iState] = std::min(wallDistance[iState], minimumDistance);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVector2CCMPI::execute2D()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_always_assert(_boundaryTRS.size() > 0);
  cf_always_assert(dim == DIM_2D);
  
  // loop over all processors
  for (CFuint root = 0; root < m_nbProc; ++root) {
    
    // loop over all wall boundary TRSs
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size(); ++iTRS) {
      SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs(_boundaryTRS[iTRS]);
      
      // broadcast number of faces in the current TRS to allow for buffer memory preallocation
      CFuint nbLocalTrsFaces = faces->getLocalNbGeoEnts();
      
#ifdef CF_HAVE_MPI
      MPI_Bcast(&nbLocalTrsFaces, 1, MPI_UNSIGNED, root, m_comm);
#endif
      
      CFuint nodeSize = 0;
      CFuint connSize = 0;
      // something is gonna be packed to be sent only if the number of faces 
      // inside the current TRS in the root processor is positive
      if (nbLocalTrsFaces > 0) {
	int sizes[3];
	TRSDistributeData trsData; // data to be distributed
	
	if (m_myRank == root) {
	  SafePtr<vector<CFuint> > nodesInTrs = faces->getNodesInTrs();
	  const CFuint nbNodesInTrs = nodesInTrs->size();
	  nodeSize = nbNodesInTrs*dim;
	  // connectivity size overestimated to allow for memory preallocation in 3D
	  connSize = (dim == DIM_2D) ? nbLocalTrsFaces*2 : nbLocalTrsFaces*4;
	  
	  // unidimensional node coordinates storage (each node is unique)
	  trsData.trsNodes.reserve(nodeSize);
	  trsData.trsNodeConn.reserve(connSize);
	  // 3D mesh could have both triangles and quads faces, so we keep track of
	  // the number of nodes for each individual face
	  trsData.trsNbNodesInFace.reserve(nbLocalTrsFaces);
	  
	  // map local node IDs (within the mesh partition) with the local ID 
	  // in the node storage for the current TRS (in the range [0,nbNodesInTrs])
	  CFMap<CFuint, CFuint> localNodeIDs2TrsNodeIDs(nbNodesInTrs);
	  for (CFuint n = 0; n < nbNodesInTrs; ++n){
	    const CFuint nodeID = (*nodesInTrs)[n];
	    localNodeIDs2TrsNodeIDs.insert(nodeID,n);
	    // store coordinates of the node in "n" position
	    for (CFuint iDim =0; iDim < dim; ++iDim) {
	      trsData.trsNodes.push_back((*nodes[nodeID])[iDim]);
	    }
	  }
	  localNodeIDs2TrsNodeIDs.sortKeys();
	  cf_assert(trsData.trsNodes.size() == trsData.trsNodes.capacity());
	  
	  // loop over faces in the current TRS and store connectivity info into TRSDistributeData
	  for (CFuint iFace = 0; iFace < nbLocalTrsFaces; ++iFace) {
	    const CFuint nbNodesInGeo = faces->getNbNodesInGeo(iFace);
	    trsData.trsNbNodesInFace.push_back(nbNodesInGeo);
	    for (CFuint iNode = 0; iNode < nbNodesInGeo; ++iNode) {
	      const CFuint nodeID = faces->getNodeID(iFace, iNode);
	      // here you push back the ID of the node within the node storage
	      // corresponding to the current TRS
	      bool found = false;
	      trsData.trsNodeConn.push_back(localNodeIDs2TrsNodeIDs.find(nodeID, found));
	      cf_assert(found);
	    }
	  }
	  
	  sizes[0] = nodeSize;
	  sizes[1] = connSize;
	  sizes[2] = nbLocalTrsFaces;
	}
	
#ifdef CF_HAVE_MPI
	// broadcast the sizes of the storages to make possible to preallocate the distributed arrays
	MPI_Bcast(&sizes[0], 3, MPI_INT, root, m_comm);
#endif
	
	if (m_myRank != root) {
	  cf_assert(sizes[0] > 0);
	  cf_assert(sizes[1] > 0);
	  cf_assert(sizes[2] > 0);
	  
	  trsData.trsNodes.resize(sizes[0]);
	  trsData.trsNodeConn.resize(sizes[1]);
	  trsData.trsNbNodesInFace.resize(sizes[2]);
	}
	
#ifdef CF_HAVE_MPI
	MPIStruct ms;
	MPIStructDef::buildMPIStruct<double,CFuint,CFuint>(&trsData.trsNodes[0],
							   &trsData.trsNodeConn[0],
							   &trsData.trsNbNodesInFace[0], sizes, ms);
	MPI_Bcast(ms.start, 1, ms.type, root, m_comm);
	
#endif
	// compute the distance to the wall starting from broadcast data 
	computeWallDistance2D(trsData);
      }
    }
  }
  
  // printToFile();
}
    
//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVector2CCMPI::computeWallDistance2D(Framework::TRSDistributeData& data)
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  const CFuint nbStates = states.size();
  CFreal minimumDistance = MathTools::MathConsts::CFrealMax();
  RealVector faceVector10(dim);
  //  RealVector faceVector20(dim);
  RealVector node0(dim, (CFreal*)CFNULL);
  RealVector node1(dim, (CFreal*)CFNULL);
  //   RealVector node2(dim, (CFreal*)CFNULL);
  RealVector nodeStateVector(dim);
  RealVector projectCoord(dim);
  
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    CFuint nodeCount0 = 0;
    CFuint nodeCount1 = nodeCount0 + 1;
    // CFuint nodeCount2 = nodeCount1 + 1;
    
    // initialize the minimum distance to a huge value
    minimumDistance = MathTools::MathConsts::CFrealMax();
    
    const CFuint nbFaces = data.trsNbNodesInFace.size();
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      cf_assert(nodeCount0 < data.trsNodeConn.size());
      cf_assert(nodeCount1 < data.trsNodeConn.size());
      const CFuint nodeID0 = data.trsNodeConn[nodeCount0];
      const CFuint nodeID1 = data.trsNodeConn[nodeCount1];
      const CFuint coordID0 = nodeID0*dim;
      const CFuint coordID1 = nodeID1*dim;
      cf_assert(coordID0 < data.trsNodes.size());
      cf_assert(coordID1 < data.trsNodes.size());
      
      node0.wrap(dim, &data.trsNodes[coordID0]);
      node1.wrap(dim, &data.trsNodes[coordID1]);
      const RealVector& stateCoord = states[iState]->getCoordinates();
      
      faceVector10 = node1 - node0;
      nodeStateVector = stateCoord - node0;
      const CFuint faceVector10Norm = faceVector10.norm2();
      const CFreal project = MathFunctions::innerProd(faceVector10, nodeStateVector)/faceVector10Norm;
      if (project < faceVector10Norm) {
	faceVector10.normalize();
	projectCoord = node0 + project * faceVector10;
	const CFreal stateFaceDistance = MathFunctions::getDistance(projectCoord, stateCoord);
	if(stateFaceDistance < minimumDistance) minimumDistance = stateFaceDistance;
      }
      
      const CFuint nbNodesNFace = data.trsNbNodesInFace[iFace];
      nodeCount0 += nbNodesNFace;
      nodeCount1 += nbNodesNFace;
      // nodeCount2 += 4; // this is safe because it is unused in 2D
      // this is consistent with definition of "connSize"
    }
    
    wallDistance[iState] = std::min(wallDistance[iState], minimumDistance);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
