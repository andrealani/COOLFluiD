// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "MeshTools/MeshToolsFVM.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "MeshTools/ComputeWallDistanceFVMCC.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeWallDistanceFVMCC, DataProcessingData, MeshToolsFVMModule> ComputeWallDistanceFVMCCProvider("ComputeWallDistanceFVMCC");

/////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceFVMCC::ComputeWallDistanceFVMCC(const std::string& name) :
  ComputeWallDistance(name),
  socket_gstates("gstates"),
  _cellBuilder(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////

ComputeWallDistanceFVMCC::~ComputeWallDistanceFVMCC()
{
}

/////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeWallDistanceFVMCC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = ComputeWallDistance::needsSockets();

  result.push_back(&socket_gstates);

  return result;
}

//////////////////////////////////////////////////////////////////////

void ComputeWallDistanceFVMCC::setup()
{
  CFAUTOTRACE;

  ComputeWallDistance::setup();

  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
  cf_assert(fvmcc.isNotNull());

  _cellBuilder = fvmcc->getData()->getCellTrsGeoBuilder();
}

//////////////////////////////////////////////////////////////////////

void ComputeWallDistanceFVMCC::execute()
{
  CFAUTOTRACE;

  CFLog(INFO, "Computing wall distances using simple state-node distance...\n");
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();

  const CFuint nbNodes = nodes.size();
  RealVector nodalDistances(nbNodes);
  CFreal minimumDistance;
  CFreal distance;

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
  {
    minimumDistance = MathTools::MathConsts::CFrealMax();
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size() ; ++iTRS)
    {
      Common::SafePtr<std::vector<CFuint> > nodesInTrs = MeshDataStack::getActive()->getTrs(_boundaryTRS[iTRS])->getNodesInTrs();
      for(CFuint iBoundaryNode = 0; iBoundaryNode < nodesInTrs->size() ; ++iBoundaryNode)
      {
//       std::cout << "Comparing node: " << *(states[iState]->getCoordinates()) << " with node: " << *(nodes[(*nodesInTrs)[iBoundaryNode]]) << std::endl;
      _tmpVector = *(nodes[(*nodesInTrs)[iBoundaryNode]]);
      _tmpVector -= *(nodes[iNode]);
      distance = _tmpVector.norm2();
      minimumDistance = min(distance, minimumDistance);
      }
    }
    nodalDistances[iNode] = minimumDistance;
  }

  // set the builders of cells
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  Common::SafePtr<CellTrsGeoBuilder> geoBuilderPtr = _cellBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder->getDataGE();
  geoData.trs = cells;

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {

    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = _cellBuilder->buildGE();
    const CFuint nbNodesInCell = currCell->getNodes()->size();

    CFreal averageDistance = 0.;
    for(CFuint iNode=0; iNode < nbNodesInCell; iNode++)
    {
      const CFuint nodeID = currCell->getNode(iNode)->getLocalID();
      averageDistance += nodalDistances[nodeID];
    }
    averageDistance /= nbNodesInCell;
    wallDistance[iCell] = averageDistance;
    _cellBuilder->releaseGE();
  }

  printToFile();
  CFLog(INFO, "Wall distances computation finished...\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
