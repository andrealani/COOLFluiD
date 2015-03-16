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
#include "Common/CFPrintContainer.hh"

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

void ComputeWallDistanceVector2CCMPI::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >
     ("CentroidBased", "Flag to select algorithm based on wall face centroid (limited usability!).");
}
    
//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVector2CCMPI::ComputeWallDistanceVector2CCMPI(const std::string& name) :
  ComputeWallDistance(name),
#ifdef CF_HAVE_MPI
  m_comm(MPI_COMM_WORLD),
#endif
  m_myRank(0),
  m_nbProc(1)
{
  addConfigOptionsTo(this);
  
  _centroidBased = false;
  setParameter("CentroidBased",&_centroidBased);
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
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // AL: gory fix to use centroid-based algorithm 
  if (_centroidBased && dim == DIM_3D) {execute3D(); return;}
  RealVector faceCentroidArray(dim);
  m_minStateFaceDistance.resize(socket_states.getDataHandle().size());
  m_minStateFaceDistance = MathTools::MathConsts::CFrealMax();
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle <CFreal> normals = socket_normals.getDataHandle();
  cf_always_assert(_boundaryTRS.size() > 0);
  cf_always_assert(dim == DIM_2D || dim == DIM_3D);
  
  // loop over all processors
  for (CFuint root = 0; root < m_nbProc; ++root) {
    
    // loop over all wall boundary TRSs
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size(); ++iTRS) {
      SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs(_boundaryTRS[iTRS]);
      
      // broadcast number of faces in the current TRS to allow for buffer memory preallocation
      CFuint nbLocalTrsFaces = faces->getLocalNbGeoEnts();
      
#ifdef CF_HAVE_MPI
      MPI_Bcast(&nbLocalTrsFaces, 1, MPIStructDef::getMPIType(&nbLocalTrsFaces), root, m_comm);
#endif
      
      CFuint nodeSize = 0;
      CFuint connSize = 0;
      // something is gonna be packed to be sent only if the number of faces 
      // inside the current TRS in the root processor is positive
      if (nbLocalTrsFaces > 0) {
	int sizes[5];
	TRSFaceDistributeData trsData; // data to be distributed
	
	if (m_myRank == root) {
	  // list of the local node IDs which are referenced by this TRS
	  SafePtr<vector<CFuint> > nodesInTrs = faces->getNodesInTrs();
	  const CFuint nbNodesInTrs = nodesInTrs->size();
	  cf_assert(nbNodesInTrs > 0);
	  nodeSize = nbNodesInTrs*dim;
	  // connectivity size overestimated to allow for memory preallocation in 3D
	  connSize = (dim == DIM_2D) ? nbLocalTrsFaces*2 : nbLocalTrsFaces*4;
	  
	  // unidimensional node coordinates storage (each node is unique)
	  // coordinate x,y(,z) of the TRS nodes
	  trsData.trsNodes.reserve(nodeSize);
	  // face-node connectivity in the TRS 
	  trsData.trsNodeConn.reserve(connSize);
	  // 3D mesh could have both triangle and quad faces, so we keep track of
	  // the number of nodes for each individual face
	  trsData.trsNbNodesInFace.reserve(nbLocalTrsFaces);
	  
	  // map local node IDs (within the mesh partition) with the local ID 
	  // in the node storage for the current TRS (in the range [0,nbNodesInTrs-1])
	  CFMap<CFuint, CFuint> localNodeIDs2TrsNodeIDs(nbNodesInTrs);
	  for (CFuint n = 0; n < nbNodesInTrs; ++n){
	    const CFuint nodeID = (*nodesInTrs)[n];
	    localNodeIDs2TrsNodeIDs.insert(nodeID,n);
	    // store coordinates of the node (with ID = nodeID) in "n" position
	    for (CFuint iDim =0; iDim < dim; ++iDim) {
	      trsData.trsNodes.push_back((*nodes[nodeID])[iDim]);
	    }
	  }
	  localNodeIDs2TrsNodeIDs.sortKeys();
	  cf_assert(trsData.trsNodes.size() == trsData.trsNodes.capacity());
	  
	  // loop over faces in the current TRS and store related data into TRSFaceDistributeData
	  trsData.faceCenters.reserve(nbLocalTrsFaces*dim);
	  for (CFuint iFace = 0; iFace < nbLocalTrsFaces; ++iFace) {
	    const CFuint nbNodesInGeo = faces->getNbNodesInGeo(iFace);
	    cf_assert((nbNodesInGeo == 2 && dim == DIM_2D) || 
		      ((nbNodesInGeo == 3 || nbNodesInGeo == 4) && dim == DIM_3D)); 
	    trsData.trsNbNodesInFace.push_back(nbNodesInGeo);
	    
	    // compute the face centroid
	    faceCentroidArray = 0.;
	    for (CFuint iNode = 0; iNode < nbNodesInGeo; ++iNode) {
	      const CFuint nodeID = faces->getNodeID(iFace, iNode);
	      faceCentroidArray += *nodes[nodeID];
	      // here you push back the ID of the node within the node storage
	      // corresponding to the current TRS
	      bool found = false;
	      trsData.trsNodeConn.push_back(localNodeIDs2TrsNodeIDs.find(nodeID, found));
	      cf_assert(found);
	    }
	    // fix for 3D to account for possibly oversized storage: 
	    // add "-1" as last ID for the current face
	    if (dim == DIM_3D && nbNodesInGeo < 4) {trsData.trsNodeConn.push_back(-1);}
	    
	    const CFreal ovNbNodesInGeo = 1./(CFreal)nbNodesInGeo;
	    faceCentroidArray *= ovNbNodesInGeo;
	    
	    for (CFuint iDim = 0; iDim < dim; ++iDim) {
	      trsData.faceCenters.push_back(faceCentroidArray[iDim]);
	    }
	  }
	  cf_assert(trsData.faceCenters.size() == trsData.faceCenters.capacity());
	  
	  sizes[0] = nodeSize; 
	  sizes[1] = connSize;
	  sizes[2] = nbLocalTrsFaces;
	  sizes[3] = sizes[4] = nbLocalTrsFaces*dim;
	  
	  // face normal is needed for computing projection
	  trsData.faceNormals.reserve(nbLocalTrsFaces*dim);
	  for (CFuint iFace = 0; iFace < nbLocalTrsFaces; ++iFace) {
	    const CFuint start = faces->getLocalGeoID(iFace)*dim;
	    for (CFuint iDim = 0; iDim < dim; ++iDim) {
	      trsData.faceNormals.push_back(normals[start+iDim]);
	    }
	  }
	  cf_assert(trsData.faceNormals.size() == trsData.faceNormals.capacity());
	}
	
#ifdef CF_HAVE_MPI
	// broadcast the sizes of the storages to make possible to preallocate the distributed arrays
	MPI_Bcast(&sizes[0], 5, MPIStructDef::getMPIType(&sizes[0]), root, m_comm);
#endif
	
	if (m_myRank != root) {
	  cf_assert(sizes[0] > 0);
	  cf_assert(sizes[1] > 0);
	  cf_assert(sizes[2] > 0);
	  cf_assert(sizes[3] > 0);
	  cf_assert(sizes[4] > 0);
	  
	  trsData.trsNodes.resize(sizes[0]);
	  trsData.trsNodeConn.resize(sizes[1]);
	  trsData.trsNbNodesInFace.resize(sizes[2]);
	  trsData.faceCenters.resize(sizes[3]);
	  trsData.faceNormals.resize(sizes[4]);
	}
	
#ifdef CF_HAVE_MPI
	MPIStruct ms;
	MPIStructDef::buildMPIStruct<CFreal,CFint,CFuint,CFreal,CFreal>
	  (&trsData.trsNodes[0], &trsData.trsNodeConn[0], &trsData.trsNbNodesInFace[0], 
	   &trsData.faceCenters[0], &trsData.faceNormals[0], sizes, ms);
	
	MPI_Bcast(ms.start, 1, ms.type, root, m_comm);
#endif
	// compute the distance to the wall starting from broadcast data 
	computeWallDistance(trsData);
      }
    }
  }
  
  if (m_nbProc == 1) {
    printToFile();
  }  
  
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
      MPI_Bcast(&nbLocalTrsFaces, 1, MPIStructDef::getMPIType(&nbLocalTrsFaces), root, m_comm);
#endif
      
      // something is gonna be packed to be sent only if the number of faces 
      // inside the current TRS in the root processor is positive
      if (nbLocalTrsFaces > 0) {
	// trsFaceData stores the coordinates of the wall face centroid
	const CFuint sizeData = nbLocalTrsFaces*dim;
	trsFaceData.resize(sizeData); 
		
	if (m_myRank == root) {
	  // loop over faces in the current TRS and store connectivity info into TRSFaceDistributeData
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
	// CFLog(INFO, CFPrintContainer<vector<CFreal> >("trsFaceData    = ", &trsFaceData));
#endif
	
	// compute the distance to the wall starting from broadcast data 
	computeWallDistance3D(trsFaceData);
      }
    }
  }
  
  if (m_nbProc == 1) {
    printToFile();
  }
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

void ComputeWallDistanceVector2CCMPI::computeWallDistance(TRSFaceDistributeData& data)
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  const CFuint nbStates = states.size();
  RealVector node0(dim, (CFreal*)CFNULL);
  RealVector fcenter(dim, (CFreal*)CFNULL);
  RealVector fnormal(dim, (CFreal*)CFNULL);
  RealVector nodeStateVector(dim);
  
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    CFuint nodeCount0 = 0;
    bool updateDistance = false;
    
    // initialize the minimum distance to a huge value
    CFreal minimumDistance = MathTools::MathConsts::CFrealMax();
    cf_assert(minimumDistance > 0.);
    
    const CFuint nbFaces = data.trsNbNodesInFace.size();
    cf_assert(nbFaces > 0);
    cf_assert(nbFaces == data.faceCenters.size()/dim);
    cf_assert(nbFaces == data.faceNormals.size()/dim);
    
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      cf_assert(nodeCount0 < data.trsNodeConn.size());
      const CFuint nodeID0 = data.trsNodeConn[nodeCount0];
      const CFuint coordID0 = nodeID0*dim;
      cf_assert(coordID0 < data.trsNodes.size());
      const CFuint faceID0 = iFace*dim;
      cf_assert(faceID0 < data.faceNormals.size());
      cf_assert(faceID0 < data.faceCenters.size());
      
      node0.wrap(dim, &data.trsNodes[coordID0]);
      fcenter.wrap(dim, &data.faceCenters[faceID0]);
      RealVector& centerCoord = states[iState]->getCoordinates();
      // if (m_myRank == 0) cout << iState <<  " => " << centerCoord << endl;
      nodeStateVector = centerCoord - node0;
      // if distance between the current cell center and the boundary face center is minimal
      // compute projection of vector joining cell center & first face node onto the face normal direction
      const CFreal stateFaceDistance = MathFunctions::getDistance(centerCoord, fcenter);
      cf_assert(m_minStateFaceDistance[iState] > 0.);
      if (stateFaceDistance < m_minStateFaceDistance[iState]) {
	fnormal.wrap(dim, &data.faceNormals[faceID0]);
	minimumDistance = std::abs(-MathFunctions::innerProd(fnormal, nodeStateVector)/fnormal.norm2());
	m_minStateFaceDistance[iState] = stateFaceDistance;
	updateDistance = true;
      }
      
      nodeCount0 += (dim == DIM_3D) ? 4 : 2; // this is consistent with definition of "connSize"
    }
    
    if (updateDistance) {wallDistance[iState] = minimumDistance;} // std::min(wallDistance[iState], minimumDistance);
    CFLog(DEBUG_MAX, "ComputeWallDistanceVector2CCMPI::computeWallDistance3D() => " <<
	  "wallDistance[ " <<  iState << " ] = " << wallDistance[iState] << "\n");
  }
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
