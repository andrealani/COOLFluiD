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
  socket_solPntNormals("solPntNormals")
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
  const CFuint nbrElemTypes = elemType->size();
  
  std::vector< std::vector< std::vector< CFuint > > > dimList(nbrElemTypes);
  
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();
    
    totNbrStates += (*elemType)[iElemType].getNbElems()*nbSolPnts;
  
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
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD


