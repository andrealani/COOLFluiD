// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/DirPaths.hh"

#include "Framework/BadFormatException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "FluxReconstructionMethod/StdSetup.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BasePointDistribution.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"

#include "Common/ConnectivityTable.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSetup,FluxReconstructionSolverData,FluxReconstructionModule >
  stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_nstates("nstates"),
  socket_nstatesProxy("nstatesProxy"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_nodes("nodes"),
  socket_gradients("gradients"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  socket_normals("normals")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StdSetup::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  StdSetup::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_nstatesProxy);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gradients);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "StdSetup::execute() => START\n");
  
  // number of equations
  const CFuint nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get datahandle for nodes
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  // get datahandle for states in nodes and resize
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  nstates.resize(nodes.size());

  // create states in nodes
  for (CFuint i = 0; i < nstates.size(); ++i) {
    nstates[i].resize(nbrEqs);
  }

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  DataHandle<ProxyDofIterator< RealVector >* > nstatesProxy = socket_nstatesProxy.getDataHandle();

  // const CFuint nbStates = states.size();
  nstatesProxy.resize(1);
  
  const CFuint nbStates = states.size();
  
  // set node to state mapping
  // m_nodeIDToStateID.resize(nbStates);
  // for (CFuint stateID=0;stateID<nbStates;++stateID) {
  //   //  const CFuint nodeID = states[stateID]->getCoordinates().getLocalID();
  //   // cf_assert(nodeID < nbStates);
  //   m_nodeIDToStateID[nodeID] = stateID;
  // }
  // nstatesProxy[0] =
  //   new DofDataHandleIterator< RealVector,State, GLOBAL >(states,&m_nodeIDToStateID);
  
  // get number of gradients (assumed equal to the number of physical variables)
  const CFuint nbrGrads = nbrEqs;
  
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated.resize(nbStates);
  isUpdated = false;
  
  // dimensionality
  const CFuint dim = PhysicalModelStack::getActive()->getDim();  

  // get datahandle
  DataHandle< std::vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // resize gradients
  gradients.resize(nbStates);
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    gradients[iState].resize(nbrGrads);
    for (CFuint iGrad = 0; iGrad < nbrGrads; ++iGrad) {
      gradients[iState][iGrad].resize(dim);
    }
  }
  
  // CREATE ADDITIONAL MESH DATASTRUCTURE

  // get the start indexes of the range of faces with a certain orientation
  createFaceOrientationStartIndexes();

  // compute the unit normals and the size of the projection vectors in face flux points
  computeFaceJacobianVectorSizes();

  // set the boundary condition type of the boundary faces
  setBndFacesBCType();
    
  CFLog(VERBOSE, "StdSetup::execute() => END\n");
}

/////////////////////////////////////////////////////////////////////////////

void StdSetup::computeFaceJacobianVectorSizes()
{
  CFAUTOTRACE;

  // get TRS list
  vector<SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  // number of TRSs
  const CFuint nbTRSs = trsList.size();

  // count total number of faces
  CFuint totNbrFaces = 0;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    if (trsList[iTRS]->getName() != "InnerCells")
    {
      totNbrFaces += trsList[iTRS]->getLocalNbGeoEnts();
    }
  }

  // get the face flux point projection vector size data handle
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();

  // resize the datahandle for the face Jacobian vector sizes
  faceJacobVecSizeFaceFlxPnts.resize(totNbrFaces);

  // face flux points data
  std::vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  SafePtr< std::vector< RealVector > >                 faceFlxPntsFaceLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();

  const CFuint nbrFaceFlxPnts = faceFlxPntsFaceLocalCoords->size();

  // get Inner Cells TRS
  SafePtr<TopologicalRegionSet> cellTRS = MeshDataStack::getActive()->getTrs("InnerCells");

  // get face builder and set cells trs
  SafePtr< GeometricEntityPool<FaceToCellGEBuilder> > geoBuilder = getMethodData().getFaceBuilder();
  FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cellTRS;

  // set face jacobian vector sizes in face flux points
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    if (trsList[iTRS]->getName() != "InnerCells")
    {
      // get number of faces in this TRS
      const CFuint nbrFaces = trsList[iTRS]->getLocalNbGeoEnts();

      // set face trs
      geoData.facesTRS = trsList[iTRS];
      geoData.isBoundary = trsList[iTRS]->hasTag("boundary") || trsList[iTRS]->hasTag("partition");

      // loop over faces
      for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
      {
        // build the GeometricEntity
        geoData.idx = iFace;
        GeometricEntity& face = *geoBuilder->buildGE();

        // get face ID
        const CFuint faceID = face.getID();

        // get the face Jacobian vectors
        vector< RealVector > faceJacobVecs = face.computeFaceJacobDetVectorAtMappedCoords(*faceFlxPntsFaceLocalCoords);

        // loop over face flux points
        for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
        {
          faceJacobVecSizeFaceFlxPnts[faceID].push_back(faceJacobVecs[iFlx].norm2());
        }

        // release the GeometricEntity
        geoBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::createFaceOrientationStartIndexes()
{
  CFAUTOTRACE;

  // START INDEXES FOR INNER CFGeoEnt::CELLS
  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get connectivity where the face start indexes are stored
  SafePtr< ConnectivityTable< CFuint > > innerFacesStartIdxsConnTable
    = MeshDataStack::getActive()->getConnectivity("innerFacesStartIdxs");

  // get number of orientations
  const CFuint nbrFaceOrientsPlus1 = innerFacesStartIdxsConnTable->nbRows();

  // resize innerFacesStartIdxs
  innerFacesStartIdxs.resize(nbrFaceOrientsPlus1);

  // put indexes in innerFacesStartIdxs
  for (CFuint iIdx = 0; iIdx < nbrFaceOrientsPlus1; ++iIdx)
  {
    innerFacesStartIdxs[iIdx] = (*innerFacesStartIdxsConnTable)(iIdx,0);
  }

  // START INDEXES FOR BOUNDARY CFGeoEnt::CELLS
  // get bndFacesStartIdxs from FluxReconstructionSolverData
  map< std::string , vector< vector< CFuint > > >& bndFacesStartIdxs = getMethodData().getBndFacesStartIdxs();

  // get TRS list
  vector<Common::SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  // number of TRSs
  const CFuint nbTRSs = trsList.size();

  // loop over boundary TRSs
  CFuint iBndTRS = 0;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    const std::string trsName = trsList[iTRS]->getName();

    if (trsList[iTRS]->hasTag("boundary"))
    {

      // get connectivity where the face start indexes are stored
      SafePtr< ConnectivityTable< CFuint > > boundaryFacesStartIdxsConnTable
        = MeshDataStack::getActive()->getConnectivity(trsName+"boundaryFacesStartIdxs");

      // get number of TRs in this TRS
      const CFuint nbTRs = trsList[iTRS]->getNbTRs();
      // number of boundary face orientations + 1
      const CFuint nbBndFaceOrientP1 = boundaryFacesStartIdxsConnTable->nbRows()/nbTRs;

      // array over TRs and startIdxs
      vector< vector< CFuint > > trStartIdxs(nbTRs);

      // loop over TRs
      CFuint iStartIdx = 0;
      for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
      {
        // resize trStartIdxs[iTR]
        trStartIdxs[iTR].resize(nbBndFaceOrientP1);

        // loop over start indexes
        for (CFuint iIdx = 0; iIdx < nbBndFaceOrientP1; ++iIdx, ++iStartIdx)
        {
          trStartIdxs[iTR][iIdx] = (*boundaryFacesStartIdxsConnTable)(iStartIdx,0);
        }
      }
      cf_assert(iStartIdx == boundaryFacesStartIdxsConnTable->nbRows());

      // store trStartIdxs in boundary TRS map
      bndFacesStartIdxs[trsName] = trStartIdxs;

      // increment boundary TRS counter
      ++iBndTRS;
    } 
  }

  // START INDEXES FOR Partition CFGeoEnt::CELLS
  // get the face start indexes
  vector< CFuint >& partitionFacesStartIdxs = getMethodData().getPartitionFacesStartIdxs();

  // get connectivity where the face start indexes are stored
  SafePtr< ConnectivityTable< CFuint > > partitionFacesStartIdxsConnTable
    = MeshDataStack::getActive()->getConnectivity("partitionFacesStartIdxs");

  // get number of orientations
  const CFuint nbrPartFaceOrientsPlus1 = partitionFacesStartIdxsConnTable->nbRows();

  // resize innerFacesStartIdxs
  partitionFacesStartIdxs.resize(nbrPartFaceOrientsPlus1);

  // put indexes in innerFacesStartIdxs
  for (CFuint iIdx = 0; iIdx < nbrPartFaceOrientsPlus1; ++iIdx)
  {
    partitionFacesStartIdxs[iIdx] = (*partitionFacesStartIdxsConnTable)(iIdx,0);
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setBndFacesBCType()
{
  CFAUTOTRACE;

  // get TRS list
  vector<SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  // number of TRSs
  const CFuint nbTRSs = trsList.size();

  // count total number of faces
  CFuint totNbrFaces = 0;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    if (trsList[iTRS]->getName() != "InnerCells")
    {
      totNbrFaces += trsList[iTRS]->getLocalNbGeoEnts();
    }
  }

  // get faceToInFaceIdxOrientOrBCIdx connectivity
  SafePtr< ConnectivityTable<CFuint> > faceToInFaceIdxOrientOrBCIdx
      = MeshDataStack::getActive()->getConnectivity("faceToInFaceIdxOrientOrBCIdx");

  // get Inner Cells TRS
  SafePtr<TopologicalRegionSet> cellTRS = MeshDataStack::getActive()->getTrs("InnerCells");

  // get bc TRS names variable from the method data
  SafePtr< vector< vector< std::string > > > bcTRSNames = getMethodData().getBCTRSNameStr();

  // loop over TRSs to set BC types and booleans
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    if (trsList[iTRS]->getName() != "InnerCells" && trsList[iTRS]->getName() != "InnerFaces")
    {
      // get number of faces in this TRS
      const CFuint nbrFaces = trsList[iTRS]->getLocalNbGeoEnts();

      // prepares to loop over faces by getting the GeometricEntityPool
      SafePtr< GeometricEntityPool<FaceToCellGEBuilder> > geoBuilder = getMethodData().getFaceBuilder();

      FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
      geoData.facesTRS = trsList[iTRS];
      geoData.cellsTRS = cellTRS;
      geoData.isBoundary = true;

      // determine boundary condition index
      CFuint bcIdx = 0;
      bool foundBCIdx = false;
      for (CFuint iBC = 0; iBC < bcTRSNames->size() && !foundBCIdx; ++iBC)
      {
        for (CFuint iBCTRS = 0; iBCTRS < (*bcTRSNames)[iBC].size() && !foundBCIdx; ++iBCTRS)
        {
          if ((*bcTRSNames)[iBC][iBCTRS] == trsList[iTRS]->getName())
          {
            foundBCIdx = true;
            bcIdx = iBC;
          }
        }
      }

      // check if boundary TRS has been found
      if (!foundBCIdx && !trsList[iTRS]->hasTag("partition"))
      {
	//throw Common::ShouldNotBeHereException (FromHere(),"Boundary condition corresponding to boundary face not found...");
      }

      // loop over faces
      for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
      {
        // build the GeometricEntity
        geoData.idx = iFace;
        GeometricEntity& face = *geoBuilder->buildGE();

        // set BC type
        cf_assert(faceToInFaceIdxOrientOrBCIdx->nbCols(face.getID()) == 1);
        (*faceToInFaceIdxOrientOrBCIdx)(face.getID(),0) = bcIdx;

        // release the GeometricEntity
        geoBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setInvGeoJacobMatrixXJacobDet(vector< RealMatrix >& matr, GeometricEntity& cell,
                                             const vector< RealVector >& mappedCoord)
{
  // get the dimensionality
  const CFuint dim = cell.getDimensionality();

  // number of points
  const CFuint nbrPnts = mappedCoord.size();

  // get Jacobian matrices and determinants
  vector< RealMatrix > jacobMatr =
      cell.computeGeometricShapeFunctionJacobianMatrix     (mappedCoord);
  std::valarray< CFreal >   jacobDet  =
      cell.computeGeometricShapeFunctionJacobianDeterminant(mappedCoord);

  // loop over points to set the matrices
  matr.resize(nbrPnts);
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    // resize matr[iPnt]
    matr[iPnt].resize(dim,dim);

    // transpose the matrix
    RealMatrix tMatr(dim,dim);
    jacobMatr[iPnt].transpose(tMatr);

    // invert the matrix
    if (dim == 2)
    {
      inverter2.invert(tMatr,matr[iPnt]);
    }
    else if (dim == 3)
    {
      inverter3.invert(tMatr,matr[iPnt]);
    }
    else
    {
      throw Common::NotImplementedException (FromHere(),"StdSetup::setGeoJacobianMatrix only implemented for 2D and 3D");
    }

    // multiply with Jacobian determinant
    matr[iPnt] *= jacobDet[iPnt];
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

