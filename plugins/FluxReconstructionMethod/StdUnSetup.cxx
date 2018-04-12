// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionMethod/StdUnSetup.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodCommandProvider< StdUnSetup,FluxReconstructionSolverData,FluxReconstructionModule >
  stdUnSetupProvider("StdUnSetup");
  
//////////////////////////////////////////////////////////////////////////////
  
StdUnSetup::StdUnSetup(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  socket_gradients("gradients"),
  socket_gradientsAV("gradientsAV"),
  socket_posPrev("posPrev")
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  StdUnSetup::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StdUnSetup::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  result.push_back(&socket_gradients);
  result.push_back(&socket_gradientsAV);
  result.push_back(&socket_posPrev);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;
  
  // Force deallocate gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  gradients.resize(0);
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();
  gradientsAV.resize(0);
  DataHandle< CFreal > posPrev = socket_posPrev.getDataHandle();
  posPrev.resize(0);
  
  // Force deallocate socket_faceJacobVecSizeFaceFlxPnts
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  faceJacobVecSizeFaceFlxPnts.resize(0);

  // VARIABLES IN FLUXRECONSTRUCTIONSOLVERDATA

  // clear frLocalData
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    delete frLocalData[iElemType];
  }
  frLocalData.resize(0);
  
  // clear start index of inner faces with a certain orientation
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
  innerFacesStartIdxs.resize(0);
  
  // clear start index of partition faces with a certain orientation
  vector< CFuint >& partitionFacesStartIdxs = getMethodData().getPartitionFacesStartIdxs();
  partitionFacesStartIdxs.resize(0);

  // clear start index of boundary faces with a certain orientation
  map<std::string, vector< vector< CFuint > > >& bndFacesStartIdxs = getMethodData().getBndFacesStartIdxs();
  bndFacesStartIdxs.clear();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

