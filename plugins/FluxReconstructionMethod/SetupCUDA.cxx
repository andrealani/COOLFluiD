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
  SetupExtra(name),
  socket_gradientsCUDA("gradientsCUDA"),
  socket_gradientsAVCUDA("gradientsAVCUDA"),        
  socket_faceDir("faceDir")
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
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = SetupExtra::providesSockets();
  result.push_back(&socket_gradientsCUDA);
  result.push_back(&socket_gradientsAVCUDA);
  result.push_back(&socket_faceDir);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SetupCUDA::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "SetupCUDA::execute() => START\n");
  
  SetupExtra::execute();
  
  createCUDASockets();
    
  CFLog(VERBOSE, "SetupCUDA::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void SetupCUDA::createCUDASockets()
{
  CFAUTOTRACE;
  
  // get dimensions
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // number of equations
  const CFuint nbrEqs = PhysicalModelStack::getActive()->getNbEq();

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
  
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();
    
    totNbrStates += (*elemType)[iElemType].getNbElems()*nbSolPnts;
    totNbrCells += (*elemType)[iElemType].getNbElems();
  }
  
  // get the face flux point projection vector size data handle
  DataHandle< CFreal > gradients = socket_gradientsCUDA.getDataHandle();
  DataHandle< CFreal > gradientsAV = socket_gradientsAVCUDA.getDataHandle();

  gradients.resize(totNbrStates*dim*nbrEqs);
  gradientsAV.resize(totNbrStates*dim*nbrEqs);

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

  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  const CFuint nbFaceFlxPnts = flxLocalCoords->size();

  // get the face flux point projection vector size data handle
  DataHandle< CFint > faceDir = socket_faceDir.getDataHandle();

  const CFuint totNbFlxPnts = frLocalData[0]->getNbrOfFlxPnts();

  // resize flx pnt normal socket
  faceDir.resize(totNbrCells*totNbFlxPnts);

  // get flux point mapped coordinate directions per orient
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDirPerOrient();

  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConnPerOrient = frLocalData[0]->getFaceFlxPntConnPerOrient();

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

      // Loop over flux points to set the normal vectors
      for (CFuint iFlxPnt = 0; iFlxPnt < nbFaceFlxPnts; ++iFlxPnt)
      {
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[orient][LEFT][iFlxPnt];
        const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[orient][RIGHT][iFlxPnt];
     
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

