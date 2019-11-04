// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/DirPaths.hh"

#include "Framework/BadFormatException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "FluxReconstructionMethod/SetupCUDA.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BasePointDistribution.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"

#include "Common/ConnectivityTable.hh"
#include "Common/NotImplementedException.hh"

#ifdef CF_HAVE_CUDA
#include "Framework/CudaDeviceManager.hh"
#include "Framework/CudaDeviceManager.hh"
#include "Common/CUDA/CFVec.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SetupCUDA,FluxReconstructionSolverData,FluxReconstructionModule >
  setupCUDAProvider("SetupCUDA");

//////////////////////////////////////////////////////////////////////////////

SetupCUDA::SetupCUDA(const std::string& name) :
  StdSetup(name),
  socket_solPntNormals("solPntNormals"),
  socket_flxPntNormals("flxPntNormals"),
  socket_faceDir("faceDir"),
  m_face(CFNULL),
  m_flxLocalCoords(),
  m_faceJacobVecs(),
  m_faceMappedCoordDir(),
  m_faceFlxPntConnPerOrient()
{
}

//////////////////////////////////////////////////////////////////////////////

SetupCUDA::~SetupCUDA()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  SetupCUDA::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = StdSetup::providesSockets();
  result.push_back(&socket_solPntNormals);
  result.push_back(&socket_flxPntNormals);
  result.push_back(&socket_faceDir);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SetupCUDA::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "SetupCUDA::execute() => START\n");
  
  StdSetup::execute();
  
  createNormalSockets();
    
  CFLog(VERBOSE, "SetupCUDA::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void SetupCUDA::createNormalSockets()
{
  CFAUTOTRACE;
  
  // get dimensions
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // get InnerCells TopologicalRegionSet
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get geometric entity builder
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get SpectralFDElementData
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // compute total number of states
  cf_assert(elemType->size() == frLocalData.size());
  CFuint totNbrStates = 0;
  CFuint totNbrCells = 0; 
  const CFuint nbrElemTypes = elemType->size();
  
  std::vector< std::vector< std::vector< CFuint > > > dimList(nbrElemTypes);
  
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();
    
    totNbrStates += (*elemType)[iElemType].getNbElems()*nbSolPnts;
    totNbrCells += (*elemType)[iElemType].getNbElems();
  
    // create a list of the dimensions in which the deriv will be calculated
    dimList[iElemType].resize(dim);
    for (CFuint iDim = 0; iDim < dim; ++iDim)
    {
      dimList[iElemType][iDim].resize(nbSolPnts);
      for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
      {
        dimList[iElemType][iDim][iSolPnt] = iDim;
      }
    }
  }
  
  // get the face flux point projection vector size data handle
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();

  // resize the datahandle for the face Jacobian vector sizes
  solPntNormals.resize(totNbrStates*dim*dim);

CFLog(INFO, "new size: " << solPntNormals.size() << ", datahandle size: " << socket_solPntNormals.getDataHandle().size() << "\n");

  // loop over element types
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // number of states in this element type
    const CFuint nbrStates = frLocalData[iElemType]->getNbrOfSolPnts();

    // solution points mapped coordinates
    SafePtr< vector< RealVector > > solPntsLocalCoords = frLocalData[iElemType]->getSolPntsLocalCoords();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      GeometricEntity *const cell = geoBuilder->buildGE();

      // get the states in this cell
      vector< State* >* states = cell->getStates();

      for (CFuint iDim = 0; iDim < dim; ++iDim)
      {
        vector<RealVector> temp = cell->computeMappedCoordPlaneNormalAtMappedCoords(dimList[iElemType][iDim],*solPntsLocalCoords);
          
        for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
        {
          const CFuint solID = (*states)[iSol]->getLocalID();
            
          for (CFuint jDim = 0; jDim < dim; ++jDim)
          {
            solPntNormals[solID*dim*dim+iDim*dim+jDim] = temp[iSol][jDim];
          }
        }
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }

  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");
  
  // get face builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();

  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  const CFuint nbFaceFlxPnts = flxLocalCoords->size();

  // get the face flux point projection vector size data handle
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  DataHandle< CFint > faceDir = socket_faceDir.getDataHandle();

  const CFuint totNbFlxPnts = frLocalData[0]->getNbrOfFlxPnts();

  // resize flx pnt normal socket
  flxPntNormals.resize(nbFaces*nbFaceFlxPnts*dim);//innerFacesStartIdxs[nbrFaceOrients]
  faceDir.resize(totNbrCells*totNbFlxPnts);

  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();

  // get flux point mapped coordinate directions per orient
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDirPerOrient();

  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConnPerOrient = frLocalData[0]->getFaceFlxPntConnPerOrient();
  
  m_faceJacobVecs.resize(nbFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < nbFaceFlxPnts; ++iFlx)
  {
    m_faceJacobVecs[iFlx].resize(dim);
  }

  // loop over different orientations
  for (CFuint orient = 0; orient < nbrFaceOrients; ++orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[orient+1];

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      const CFuint leftCellID = m_face->getNeighborGeo(LEFT)->getID();
      const CFuint rightCellID = m_face->getNeighborGeo(RIGHT)->getID();

      // compute face Jacobian vectors
      m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

      // Loop over flux points to set the normal vectors
      for (CFuint iFlxPnt = 0; iFlxPnt < nbFaceFlxPnts; ++iFlxPnt)
      {
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[orient][LEFT][iFlxPnt];
        const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[orient][RIGHT][iFlxPnt];

        for (CFuint iDim = 0; iDim < dim; ++iDim)
        {
          // set unit normal vector
          flxPntNormals[m_face->getID()*nbFaceFlxPnts*dim+iFlxPnt*dim+iDim] = m_faceJacobVecs[iFlxPnt][iDim];
        }
     
        faceDir[leftCellID*totNbFlxPnts+flxPntIdxL] = (*m_faceMappedCoordDir)[orient][LEFT];
        faceDir[rightCellID*totNbFlxPnts+flxPntIdxR] = (*m_faceMappedCoordDir)[orient][RIGHT];
      }
      
      m_faceBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

