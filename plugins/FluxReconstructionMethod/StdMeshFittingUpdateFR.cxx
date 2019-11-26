#include "FluxReconstructionMethod/FluxReconstructionSolver.hh" 
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh" 
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/StdSetup.hh"


#include "FluxReconstructionMethod/StdMeshFittingUpdateFR.hh"

//#include "FiniteVolume/FVMCC_PolyRec.hh"

#include "Framework/ComputeDummyStates.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeCalculator.hh"
#include "Framework/Node.hh"
#include "Framework/SetElementStateCoord.hh"

#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdMeshFittingUpdateFR, FluxReconstructionSolverData, FluxReconstructionModule> 
StdMeshFittingUpdateFRProvider("StdMeshFittingUpdateFR");

//////////////////////////////////////////////////////////////////////////////

StdMeshFittingUpdateFR::StdMeshFittingUpdateFR(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_volumes("volumes"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  socket_normals("normals")
{
}

//////////////////////////////////////////////////////////////////////////////

StdMeshFittingUpdateFR::~StdMeshFittingUpdateFR()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StdMeshFittingUpdateFR::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_volumes);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  CFLog(INFO, "StdSetup NeedsSockets\n");
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void StdMeshFittingUpdateFR::execute()
{
  CFAUTOTRACE;
  CFLog(INFO, " StdMeshFittingUpdateFR::execute => start \n");

  // CREATE ADDITIONAL MESH DATASTRUCTURE
  // compute the cell volumes if necessary (actually, the Jacobian determinants in the states are put in the volumes socket)
  if (getMethodData().createVolumesSocket())
  {
    computeStatesVolumes();
  }

  // compute the unit normals and the size of the projection vectors in face flux points
  computeFaceJacobianVectorSizes();
  CFLog(INFO, " StdMeshFittingUpdateFR::execute => end \n");

}

/////////////////////////////////////////////////////////////////////////////

void StdMeshFittingUpdateFR::computeStatesVolumes()
{
  CFAUTOTRACE;

  // get the datahandle of cell volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

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
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    totNbrStates += (*elemType)[iElemType].getNbElems()*frLocalData[iElemType]->getNbrOfSolPnts();
  }

  // resize the vector for the cell volumes per state
  volumes.resize(totNbrStates);

  // loop over element types
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // number of states in this element type
    const CFuint nbrStates = frLocalData[iElemType]->getNbrOfSolPnts();

    // solution points mapped coordinates
    SafePtr< vector< RealVector > > solPntsLocalCoords
        = frLocalData[iElemType]->getSolPntsLocalCoords();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      GeometricEntity *const cell = geoBuilder->buildGE();

      // get the states in this cell
      vector< State* >* states = cell->getStates();

      std::valarray<CFreal> jacobDet = cell->computeGeometricShapeFunctionJacobianDeterminant(*solPntsLocalCoords);

      for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
      {
        const CFuint solID = (*states)[iSol]->getLocalID();
        if (jacobDet[iSol] < 0.0)
        {
          const std::string message = "Negative Jacobian determinant (" + StringOps::to_str(jacobDet[iSol]) +
                                   ") in cell with ID " + StringOps::to_str(elemIdx);
          throw BadFormatException (FromHere(),message);
        }
        volumes[solID] = jacobDet[iSol];
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

void StdMeshFittingUpdateFR::computeFaceJacobianVectorSizes()
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

    } // namespace FiniteVolume


} // namespace COOLFluiD
