// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/DirPaths.hh"

#include "Framework/BadFormatException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "FluxReconstructionMethod/SetupExtra.hh"
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

MethodCommandProvider< SetupExtra,FluxReconstructionSolverData,FluxReconstructionModule >
  setupExtraProvider("SetupExtra");

//////////////////////////////////////////////////////////////////////////////

SetupExtra::SetupExtra(const std::string& name) :
  StdSetup(name),
  socket_solPntNormals("solPntNormals"),
  socket_flxPntNormals("flxPntNormals"),
  socket_cellVolumes("cellVolumes"),
  m_face(CFNULL),
  m_flxLocalCoords(),
  m_faceJacobVecs(),
  m_faceMappedCoordDir(),
  m_faceFlxPntConnPerOrient()
{
}

//////////////////////////////////////////////////////////////////////////////

SetupExtra::~SetupExtra()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  SetupExtra::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = StdSetup::providesSockets();
  result.push_back(&socket_solPntNormals);
  result.push_back(&socket_flxPntNormals);
  result.push_back(&socket_cellVolumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SetupExtra::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "SetupExtra::execute() => START\n");
  
  StdSetup::execute();
  
  createNormalSockets();
    
  CFLog(VERBOSE, "SetupExtra::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void SetupExtra::createNormalSockets()
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
    
  //Setting ndimplus, needed for Triag (and tetra, prism)

  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();
  CFuint m_ndimplus;
  if (elemShape == CFGeoShape::TRIAG)
    {
      m_ndimplus=3;	
    }
  else if (elemShape == CFGeoShape::TETRA)
    {
      m_ndimplus=4;
    }  
  else
    {
      m_ndimplus=0;
    }

  std::vector< std::vector< std::vector< CFuint > > > dimList(nbrElemTypes);
  
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();
    
    totNbrStates += (*elemType)[iElemType].getNbElems()*nbSolPnts;
    totNbrCells += (*elemType)[iElemType].getNbElems();
  
    // create a list of the dimensions in which the deriv will be calculated
    dimList[iElemType].resize(dim+m_ndimplus);
    for (CFuint iDim = 0; iDim < dim+m_ndimplus; ++iDim)
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
  
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  cellVolumes.resize(totNbrCells);

  // resize the datahandle for the face Jacobian vector sizes
  solPntNormals.resize(totNbrStates*(dim+m_ndimplus)*dim);
  CFuint nbNegJacob = 0;
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
      
      cellVolumes[elemIdx] = cell->computeVolume();
      //cf_assert(cellVolumes[elemIdx] > 0.0);
      if (cellVolumes[elemIdx] < 0.0)
      {
          CFLog(INFO, "WARNING: NEGATIVE VOLUME FOUND!\n");
        cellVolumes[elemIdx] = -cellVolumes[elemIdx];
        nbNegJacob++;
      }

      // get the states in this cell
      vector< State* >* states = cell->getStates();

      for (CFuint iDim = 0; iDim < dim+m_ndimplus; ++iDim)
      {
        vector<RealVector> temp = cell->computeMappedCoordPlaneNormalAtMappedCoords(dimList[iElemType][iDim],*solPntsLocalCoords);
          
        for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
        {
          const CFuint solID = (*states)[iSol]->getLocalID();
            
          for (CFuint jDim = 0; jDim < dim; ++jDim)
          {
            solPntNormals[solID*(dim+m_ndimplus)*dim+iDim*dim+jDim] = temp[iSol][jDim];
          }
        }
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
  CFLog(INFO, "NB NEGATIVE VOLUMES: " << nbNegJacob << "\n");

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

  const CFuint totNbFlxPnts = frLocalData[0]->getNbrOfFlxPnts();

  // resize flx pnt normal socket
  flxPntNormals.resize(nbFaces*nbFaceFlxPnts*dim);//innerFacesStartIdxs[nbrFaceOrients]

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

      // compute face Jacobian vectors
      m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

      // Loop over flux points to set the normal vectors
      for (CFuint iFlxPnt = 0; iFlxPnt < nbFaceFlxPnts; ++iFlxPnt)
      {
        for (CFuint iDim = 0; iDim < dim; ++iDim)
        {
          // set unit normal vector
          flxPntNormals[m_face->getID()*nbFaceFlxPnts*dim+iFlxPnt*dim+iDim] = m_faceJacobVecs[iFlxPnt][iDim];
        }
      }
      
      m_faceBuilder->releaseGE();
    }
  }
  
  /// add the bnd face normals in socket flxPntNormals
  
  // get bndFacesStartIdxs from FluxReconstructionMethodData
  map< std::string , vector< vector< CFuint > > >&
      bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  
  // get bc TRS names variable from the method data
  SafePtr< std::vector< std::vector< std::string > > > bcTRSNames = getMethodData().getBCTRSNameStr();
  
  // loop over all wall boundary TRSs
  for(CFuint iBc = 0; iBc < bcTRSNames->size(); ++iBc) 
  {
  for(CFuint iTRS = 0; iTRS < (*bcTRSNames)[iBc].size(); ++iTRS) 
  {
    SafePtr<TopologicalRegionSet> faceTrs = MeshDataStack::getActive()->getTrs((*bcTRSNames)[iBc][iTRS]);
  
    vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

    // number of face orientations (should be the same for all TRs)
    cf_assert(bndFacesStartIdxs.size() != 0);
    const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

    // number of TRs
    const CFuint nbTRs = faceTrs->getNbTRs();
    cf_assert(bndFacesStartIdxs.size() == nbTRs);

    // get the geodata of the face builder and set the TRSs
    FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
    geoData.cellsTRS = cells;
    geoData.facesTRS = faceTrs;
    geoData.isBoundary = true;
  
    // loop over TRs
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
    {
      // loop over different orientations
      for (CFuint orient = 0; orient < nbOrients; ++orient)
      {   
        // start and stop index of the faces with this orientation
        const CFuint startFaceIdx = bndFacesStartIdxs[iTR][orient  ];
        const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][orient+1];

        // loop over faces with this orientation
        for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
        {
          // build the face GeometricEntity
          geoData.idx = faceID;
          m_face = m_faceBuilder->buildGE();
        
          // compute face Jacobian vectors
          m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
        
          // Loop over flux points to set the normal vectors
          for (CFuint iFlxPnt = 0; iFlxPnt < nbFaceFlxPnts; ++iFlxPnt)
          {
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              // set unit normal vector
              flxPntNormals[m_face->getID()*nbFaceFlxPnts*dim+iFlxPnt*dim+iDim] = m_faceJacobVecs[iFlxPnt][iDim];
            }
          }
        
          // release the face
          m_faceBuilder->releaseGE();
        }
      }
    }
  }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

