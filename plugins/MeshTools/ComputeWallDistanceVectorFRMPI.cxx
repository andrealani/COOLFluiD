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


#include <vector>
#include <cmath>
#include "MeshTools/MeshToolsFR.hh"
#include "MeshTools/ComputeWallDistanceVectorFRMPI.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolver.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::FluxReconstructionMethod;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeWallDistanceVectorFRMPI, DataProcessingData, MeshToolsFRModule>
ComputeWallDistanceVectorFRMPIProvider("ComputeWallDistanceVectorFRMPI");

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVectorFRMPI::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("AcceptableDistance","Distance");
}
    
//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVectorFRMPI::ComputeWallDistanceVectorFRMPI(const std::string& name) :
  ComputeWallDistance(name),
#ifdef CF_HAVE_MPI
  m_comm(MPI_COMM_WORLD),
#endif
  m_myRank(0),
  m_nbProc(1),
  socket_flxPntNormals("flxPntNormals")
{
  addConfigOptionsTo(this);
  
  m_acceptableDistance = 0.;
  setParameter("AcceptableDistance",&m_acceptableDistance);

}
    
//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceVectorFRMPI::~ComputeWallDistanceVectorFRMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVectorFRMPI::setup()
{
  CFAUTOTRACE;

  ComputeWallDistance::setup();
  
#ifdef CF_HAVE_MPI
  const std::string nsp = this->getMethodData().getNamespace();
  m_comm   = PE::GetPE().GetCommunicator(nsp);
  m_myRank = PE::GetPE().GetRank(nsp);
  m_nbProc = PE::GetPE().GetProcessorCount(nsp);
#endif


  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  
  // initialize the wall distance
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
  wallDistance = MathTools::MathConsts::CFrealMax();
  DataHandle<bool> nodeisAD = socket_nodeisAD.getDataHandle();
  nodeisAD.resize(0);
//  nodeisAD.resize(socket_nodes.getDataHandle().size());
//  nodeisAD=false;

  DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle();
  nodeDistance.resize(0);
//  nodeDistance.resize(socket_nodes.getDataHandle().size());
//  nodeDistance=MathTools::MathConsts::CFrealMax();
  
  
  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<FluxReconstructionSolver> fr = spaceMethod.d_castTo<FluxReconstructionSolver>();
  cf_assert(fr.isNotNull());

  m_frData = fr->getData();

}


//////////////////////////////////////////////////////////////////////////////

vector<Common::SafePtr<BaseDataSocketSink> >
ComputeWallDistanceVectorFRMPI::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ComputeWallDistance::needsSockets();
  result.push_back(&socket_flxPntNormals);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVectorFRMPI::execute()
{
  CFAUTOTRACE;
  
  // this has only to be run at setup for now
  CFLog(VERBOSE, "ComputeWallDistanceVectorFRMPI::execute() START\n");
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = m_frData->getFRLocalData();
  
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  
  const CFuint nbFaceFlxPnts = flxLocalCoords->size();
  
  // AL: gory fix to use centroid-based algorithm 
  if (dim == DIM_3D) {execute3D(); return;}
  
  CFLog(INFO, "ComputeWallDistanceVectorFRMPI::execute() => Computing distance to the wall ...\n");
  
  Stopwatch<WallTime> stp;
  stp.start();
  
  RealVector faceCentroidArray(dim);
  m_minStateFaceDistance.resize(socket_states.getDataHandle().size());
  m_minStateFaceDistance = MathTools::MathConsts::CFrealMax();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle <CFreal> normals = socket_flxPntNormals.getDataHandle();
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
	    const CFuint start = faces->getLocalGeoID(iFace)*dim*nbFaceFlxPnts;
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
        //computeWallDistanceExtrapolate(trsData);
      }
    }
  }
  
  CFLog(INFO, "ComputeWallDistanceVectorFRMPI::execute() => took " << stp.read() << "s\n");
  
  if (m_nbProc == 1) {
    printToFile();
  }  
  
  CFLog(VERBOSE, "ComputeWallDistanceVectorFRMPI::execute() END\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVectorFRMPI::execute3D()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "ComputeWallDistanceVectorFRMPI::execute3D() => Computing distance to the wall ...\n");
  
  Stopwatch<WallTime> stp;
  stp.start();
  
  DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_always_assert(_boundaryTRS.size() > 0);
  
  RealVector faceCentroid(dim);
  vector<CFreal> trsFaceData;
  
  // loop over all processors
  for (CFuint root = 0; root < m_nbProc; ++root) {
    
    // loop over all wall boundary TRSs
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size(); ++iTRS) {  
      CFLog(VERBOSE, "ComputeWallDistanceVectorFRMPI::execute() => Processing TRS [" 
	    << _boundaryTRS[iTRS] << "] from P[" << root << "] => START\n");
      
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
	cf_assert(sizeData > 0);
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
	      cf_assert(nodeID < nodes.size());
	      faceCentroid += *nodes[nodeID];
	    }
	    const CFreal ovNbNodesInFace = 1./(CFreal)nbNodesInFace;
	    faceCentroid *= ovNbNodesInFace;
	    
	    for (CFuint iDim = 0; iDim < dim; ++iDim, ++count) {
	      cf_assert(count < trsFaceData.size());
	      trsFaceData[count] = faceCentroid[iDim];
	    }
	  }
	  cf_assert(count == trsFaceData.size());
	}
	
#ifdef CF_HAVE_MPI
	cf_assert(trsFaceData.size() > 0);
	MPI_Bcast(&trsFaceData[0], sizeData, MPI_DOUBLE, root, m_comm);
	// CFLog(INFO, CFPrintContainer<vector<CFreal> >("trsFaceData    = ", &trsFaceData));
#endif
	
	// compute the distance to the wall starting from broadcast data 
	computeWallDistance3D(trsFaceData);
      }
      
      CFLog(VERBOSE, "ComputeWallDistanceVectorFRMPI::execute() => Processing TRS [" 
	    << _boundaryTRS[iTRS] << "] from P[" << root << "] => END\n");
    }
  }
  
  CFLog(INFO, "ComputeWallDistanceVectorFRMPI::execute3D() => took " << stp.read() << "s\n");
    
  if (m_nbProc == 1) {
    printToFile();
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceVectorFRMPI::computeWallDistance3D(std::vector<CFreal>& data)
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
  DataHandle<bool> nodeisAD = socket_nodeisAD.getDataHandle(); 
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbStates = states.size();
  RealVector faceCentroid(dim, (CFreal*)CFNULL);
  const CFuint nbFaces = data.size()/dim;
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  
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
    
    cf_assert(iState < wallDistance.size());
    wallDistance[iState] = std::min(wallDistance[iState], minimumDistance);
    
    
    const CFuint nbNodesInCell = cells->getNbNodesInGeo(iState);
    for (CFuint in = 0; in < nbNodesInCell; ++in) {
      // local ID of the cell node
      const CFuint cellNodeID = cells->getNodeID(iState, in);
      cf_assert(cellNodeID < nodeisAD.size());
      nodeisAD[cellNodeID] = (std::min(wallDistance[iState], minimumDistance) < m_acceptableDistance  ? true : false);
    }
  }
}
    

//////////////////////////////////////////////////////////////////////////////

  void  ComputeWallDistanceVectorFRMPI::computeWallDistance(TRSFaceDistributeData& data)
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();

 const CFuint dim = PhysicalModelStack::getActive()->getDim();
  	//cout<<"yes2"<<endl;

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  
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
      nodeStateVector = centerCoord - node0;
      // if distance between the current cell center and the boundary face center is minimal
      // compute projection of vector joining cell center & first face node onto the face normal direction
      const CFreal stateFaceDistance = MathFunctions::getDistance(centerCoord, fcenter);

      cf_assert(m_minStateFaceDistance[iState] > 0.);
      if (stateFaceDistance < m_minStateFaceDistance[iState]) {
	fnormal.wrap(dim, &data.faceNormals[faceID0]);

	minimumDistance = fabs(-MathFunctions::innerProd(fnormal, nodeStateVector)/fnormal.norm2());
	m_minStateFaceDistance[iState] = stateFaceDistance;
	updateDistance = true;
      }

      nodeCount0 += (dim == DIM_3D) ? 4 : 2; // this is consistent with definition of "connSize"
    }

    if (updateDistance) 
    { 
      wallDistance[iState] = minimumDistance; 
    } 
    
    CFLog(VERBOSE, "ComputeWallDistanceVectorFRMPI::computeWallDistance2D() => " <<
	  "miwallDistance[ " <<  iState << " ] = " << wallDistance[iState] << "\n");
  }
} 
  
//////////////////////////////////////////////////////////////////////////////
  
void  ComputeWallDistanceVectorFRMPI::computeWallDistanceExtrapolate(TRSFaceDistributeData& data){
	      //cout<<"yes"<<endl;

  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();  
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");

  const CFuint nbCells = states.size();

  RealVector Extrapolated(nodes.size()); Extrapolated = 0.;
  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    CFuint nbNodesInSideRegion = 0;
    bool foundCell = false;
    const CFuint nbNodesInCell = cells->getNbNodesInGeo(iCell);
    for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
      CFuint nodeIDinC=cells->getNodeID(iCell, iNodeC);
      if ( (nodeisAD[nodeIDinC] == true)  && Extrapolated[nodeIDinC] == 0.){
	for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
	  CFuint nodeIDinC=cells->getNodeID(iCell, iNodeC);
	  foundCell = true;
	}
	if(foundCell){
	  for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
	    CFuint nodeIDinC2=cells->getNodeID(iCell, iNodeC);
	    if(nodeisAD[nodeIDinC2]){
	      nbNodesInSideRegion= nbNodesInSideRegion+1;
	    }
	  }
	  if(nbNodesInSideRegion>= 2){
	    for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
	      CFuint nodeIDinC3=cells->getNodeID(iCell, iNodeC);
	      if(nodeIDinC3 != nodeIDinC){
		nodeisAD[nodeIDinC3] = true;
		Extrapolated[nodeIDinC3] = 1.; 
	      }
	    }
	  }
	}
      }
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
