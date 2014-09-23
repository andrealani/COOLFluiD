#include "Environment/DirPaths.hh"
#include "MathTools/MatrixInverterT.hh"

#include "Framework/BadFormatException.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SpectralFV/BaseBndFaceTermComputer.hh"
#include "SpectralFV/BaseFaceTermComputer.hh"
#include "SpectralFV/BaseVolTermComputer.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/FaceDiffusiveFlux.hh"
#include "SpectralFV/RiemannFlux.hh"
#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/StdSetup.hh"

#include "Common/ConnectivityTable.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSetup,SpectralFVMethodData,SpectralFVModule >
  stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_nstates("nstates"),
  socket_nstatesProxy("nstatesProxy"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_nodes("nodes"),
  socket_faceNormTransfMatrices("faceNormTransfMatrices"),
  socket_faceSurf("faceSurf"),
  socket_volumes("volumes"),
  socket_gradients("gradients")
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
  result.push_back(&socket_nstatesProxy           );
  result.push_back(&socket_nstates                );
  result.push_back(&socket_isUpdated              );
  result.push_back(&socket_faceNormTransfMatrices );
  result.push_back(&socket_faceSurf               );
  result.push_back(&socket_gradients              );
  result.push_back(&socket_volumes                );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  CFLog(INFO,"\n\nStdSetup::execute\n\n");

  // get datahandle for nodes
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  // get datahandle for states in nodes and resize
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  nstates.resize(nodes.size());

  // create states in nodes
  for (CFuint i = 0; i < nstates.size(); ++i) {
    nstates[i].resize(PhysicalModelStack::getActive()->getNbEq());
  }

  DataHandle<ProxyDofIterator< RealVector >* > nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  nstatesProxy.resize(1);
  nstatesProxy[0] = new DofDataHandleIterator<RealVector,RealVector>(nstates);

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  // number of states
  const CFuint nbrStates = states.size();

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated.resize(nbrStates);
  isUpdated = false;

  // create gradients variable
  if (getMethodData().hasDiffTerm())
  {
    // get number of gradients
    const CFuint nbrGrads = PhysicalModelStack::getActive()->getNbEq();

    // dimensionality
    const CFuint dim = PhysicalModelStack::getActive()->getDim();

    // get datahandle
    DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

    // resize gradients
    gradients.resize(nbrStates);
    for (CFuint iState = 0; iState < nbrStates; ++iState)
    {
      gradients[iState].resize(nbrGrads);
      for (CFuint iGrad = 0; iGrad < nbrGrads; ++iGrad)
      {
        gradients[iState][iGrad].resize(dim);
      }
    }
  }

  // CREATE ADDITIONAL MESH DATASTRUCTURE
  // compute the cell volumes if necessary
  if (getMethodData().createVolumesSocket())
  {
    computeStatesVolumes();
  }

  // compute the face surfaces
  computeFaceSurfaces();

  // get the start indexes of the range of faces with a certain orientation
  createFaceOrientationStartIndexes();

  // compute the cell face normal transformation matrices
  computeCellTransformationMatrices();

  // set the boundary condition type of the boundary faces
  setBndFacesBCType();

  // set maxNbStatesData  (must be setup after volume-, face- and boundarytermcomputers!!!)
 /// @todo broken after release 2009.3
//   PhysicalModelStack::getActive()->getImplementor()->setMaxNbStatesData(getMethodData().getMaxNbrStatesData());
}

/////////////////////////////////////////////////////////////////////////////

void StdSetup::computeStatesVolumes()
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

  // get SpectralFVElementData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // compute total number of states
  cf_assert(elemType->size() == svLocalData.size());
  CFuint totNbrStates = 0;
  const CFuint nbrElemTypes = elemType->size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    totNbrStates += (*elemType)[iElemType].getNbElems()*svLocalData[iElemType]->getNbrOfCVs();
  }

  // resize the vector for the cell volumes per state
  volumes.resize(totNbrStates);

  // loop over element types
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // number of CVs in this element type
    const CFuint nbrCVs = svLocalData[iElemType]->getNbrOfCVs();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      GeometricEntity *const cell = geoBuilder->buildGE();

      // compute the volume of the cell
      const CFreal cellVolume = cell->computeVolume();
      cf_assert(cellVolume > 0.0);

      // get the states in this cell
      vector< State* >* states = cell->getStates();

      for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
      {
        const CFuint cvID = (*states)[iCV]->getLocalID();
        volumes[cvID] = cellVolume;
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

void StdSetup::computeFaceSurfaces()
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

  // get the face surfaces data handle
  DataHandle< CFreal > faceSurf = socket_faceSurf.getDataHandle();

  // resize the datahandle for the face surfaces
  faceSurf.resize(totNbrFaces);

  // get Inner Cells TRS
  SafePtr<TopologicalRegionSet> cellTRS = MeshDataStack::getActive()->getTrs("InnerCells");

  // get face builder and set cells trs
  SafePtr< GeometricEntityPool<FaceToCellGEBuilder> > geoBuilder = getMethodData().getFaceBuilder();
  FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cellTRS;

  // loop over TRSs to compute face surfaces
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

        // compute face surface
        RealVector normal = face.computeAvgCellNormal();
        faceSurf[face.getID()] = 0.0;
        for (CFuint iCoor = 0; iCoor < normal.size(); ++iCoor)
        {
          faceSurf[face.getID()] += normal[iCoor]*normal[iCoor];
        }
        faceSurf[face.getID()] = sqrt(faceSurf[face.getID()]);
// !!!! the following command does not return the correct value in the case of a 3D triangle face!
//         faceSurf[face.getID()] = face.computeVolume();
        cf_assert(faceSurf[face.getID()] > 0.0);

        // release the GeometricEntity
        geoBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::computeCellTransformationMatrices()
{
  CFAUTOTRACE;

   MathTools::MatrixInverterT<2> inverter2;
   MathTools::MatrixInverterT<3> inverter3;

  // Get the datahandle of the transformation matrices for face normals
  DataHandle< RealMatrix > faceNormTransfMatrices = socket_faceNormTransfMatrices.getDataHandle();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get InnerCells TopologicalRegionSet
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the number of cells
  const CFuint nbrCells = cells->getLocalNbGeoEnts();

  // resize faceNormTransfMatrices
  faceNormTransfMatrices.resize(nbrCells);

  // get geometric entity builder
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // Get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // Get element shape
    const CFGeoShape::Type shape = (*elemType)[iElemType].getGeoShape();

    // Get start index of this element type in global element list
    CFuint globalIdx = (*elemType)[iElemType].getStartIdx();

    // Check if geometric shape function is not higher order than 1
    if ((*elemType)[iElemType].getGeoOrder() > 1) throw Common::ShouldNotBeHereException (FromHere(),"Only elements with a P0 or a P1 mapping to a reference element can have a linear transformation to this reference element");

    switch (shape)
    {
      case CFGeoShape::LINE:
      {
        for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
        {
          // Resize the face normal transformation matrix
          faceNormTransfMatrices[globalIdx].resize(1,1);

          // Fill the face normal transformation matrix
          faceNormTransfMatrices[globalIdx](1,1) = 1;
        }
      } break;
      case CFGeoShape::TRIAG:
      {
        for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
        {
          // build the GeometricEntity
          geoData.idx = globalIdx;
          GeometricEntity *const cell = geoBuilder->buildGE();

          // get the nodes
          vector<Node*>* nodes = cell->getNodes();

          // check if element is a triangle
          cf_assert(cell->getShape() == shape);
          cf_assert((*nodes).size() == 3);

          // variabele for regular cell transformation matrix
          RealMatrix transfM(2,2);

          // create cell transformation matrix
          for (CFuint iVec = 0; iVec < 2; ++iVec)
          {
            transfM.setColumn(RealVector((*(*nodes)[iVec+1]) - (*(*nodes)[0])),iVec);
          }

          // resize the face normal transformation matrix
          faceNormTransfMatrices[globalIdx].resize(2,2);

          // invert the cell transformation matrix and put in face normal transformation matrix
          inverter2.invert(transfM,faceNormTransfMatrices[globalIdx]);

          // transpose
          faceNormTransfMatrices[globalIdx].transpose(transfM);
          faceNormTransfMatrices[globalIdx] = transfM;

          // multiply with cell volume
          faceNormTransfMatrices[globalIdx] *= 2*cell->computeVolume();

          //release the GeometricEntity
          geoBuilder->releaseGE();
        }
      } break;
      case CFGeoShape::TETRA:
      {
        for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
        {
          // build the GeometricEntity
          geoData.idx = globalIdx;
          GeometricEntity *const cell = geoBuilder->buildGE();

          // get the nodes
          vector<Node*>* nodes = cell->getNodes();

          // check if element is a triangle
          cf_assert(cell->getShape() == shape);
          cf_assert((*nodes).size() == 4);

          // variabele for regular cell transformation matrix
          RealMatrix transfM(3,3);

          // create cell transformation matrix
          for (CFuint iVec = 0; iVec < 3; ++iVec)
          {
            transfM.setColumn(RealVector((*(*nodes)[iVec+1]) - (*(*nodes)[0])),iVec);
          }

          // resize the face normal transformation matrix
          faceNormTransfMatrices[globalIdx].resize(3,3);

          // invert the cell transformation matrix and put in face normal transformation matrix
          inverter3.invert(transfM,faceNormTransfMatrices[globalIdx]);

          // transpose
          faceNormTransfMatrices[globalIdx].transpose(transfM);
          faceNormTransfMatrices[globalIdx] = transfM;

          // multiply with cell volume
          faceNormTransfMatrices[globalIdx] *= 6*cell->computeVolume();

          //release the GeometricEntity
          geoBuilder->releaseGE();
        }
      } break;
      default:
        throw Common::NotImplementedException (FromHere(),"Only linear, triangular and tetrahedral cells have been implemented for spectral FV!");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::createFaceOrientationStartIndexes()
{
  CFAUTOTRACE;

  // START INDEXES FOR INNER CELLS
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

  // START INDEXES FOR BOUNDARY CELLS
  // get bndFacesStartIdxs from SpectralFVMethodData
  map< std::string , vector< vector< CFuint > > >& bndFacesStartIdxs = getMethodData().getBndFacesStartIdxs();

  // get TRS list
  vector<Common::SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  // number of TRSs
  const CFuint nbTRSs = trsList.size();

  // count the number of TRSs and TRs that are not InnerCells, InnerFaces or PartitionFaces
  CFuint nbBndTRSs = 0;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    if (trsList[iTRS]->hasTag("boundary"))
    {
      ++nbBndTRSs;
    }
  }

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
        throw Common::ShouldNotBeHereException (FromHere(),"Boundary condition corresponding to boundary face not found...");

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

  }  // namespace SpectralFV
}  // namespace COOLFluiD
