#include "Framework/MethodCommandProvider.hh"

#include "DiscontGalerkin/StdSetup.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSetup,DiscontGalerkinSolverData,DiscontGalerkinModule >
  stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  DiscontGalerkinSolverCom(name),
  socket_nstatesProxy("nstatesProxy"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_integrationIndex("integrationIndex"),
  socket_normals("normals"),
  socket_techStatesToNodes("techStatesToNodes"),
  socket_techNodesToStates("techNodesToStates"),
  socket_techNodesCoordinates("techNodesCoordinates")
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
  result.push_back(&socket_integrationIndex);
  result.push_back(&socket_normals);
  result.push_back(&socket_techStatesToNodes);
  result.push_back(&socket_techNodesToStates);
  result.push_back(&socket_techNodesCoordinates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
//   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle < RealVector > techNodesCoordinates = socket_techNodesCoordinates.getDataHandle();
  DataHandle < std::vector< CFuint > > techNodesToStates = socket_techNodesToStates.getDataHandle();
  DataHandle < CFuint > techStatesToNodes = socket_techStatesToNodes.getDataHandle();

  DataHandle<ProxyDofIterator< RealVector >* > nstatesProxy =
    socket_nstatesProxy.getDataHandle();
  DataHandle< std::vector< CFuint > > indexes =
    socket_integrationIndex.getDataHandle();

  DataHandle< RealVector > normals =
    socket_normals.getDataHandle();

  const CFuint nbStates = states.size();
//   const CFuint nbNodes = nodes.size();
  nstatesProxy.resize(1);

  Common::SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");
  const CFuint nbFaces = faces->getLocalNbGeoEnts();

  indexes.resize(nbFaces);
  normals.resize(nbFaces);

//   //   set node to state mapping
//   m_nodeIDToStateID.resize(nbNodes);
//   for (CFuint stateID=0;stateID<nbStates;++stateID) {
//     const CFuint nodeID = states[stateID]->getCoordinates().getLocalID();
//     cf_assert(nodeID < nbStates);
//     m_nodeIDToStateID[nodeID] = stateID;
//   }
//   nstatesProxy[0] =
//     new DofDataHandleIterator< RealVector,State >(states,&m_nodeIDToStateID);

//   m_nodeIDToStateID.resize(nbNodes);
//   for (CFuint stateID=0;stateID<nbStates;++stateID) {
//     if (states[stateID]->getCoordinates().hasLocalID()) {
//       const CFuint nodeID = states[stateID]->getCoordinates().getLocalID();
//       cf_assert(nodeID < nbNodes);
//       m_nodeIDToStateID[nodeID] = stateID;
//     }
//   }
//   nstatesProxy[0] =
//     new DofDataHandleIterator< RealVector,State >(states,&m_nodeIDToStateID);


  // seting of nodes and states for visualization

  //set trs pointer to inner cells trs
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  //get number of inner cells in region inner cells
  const CFuint nbGeos = cells->getLocalNbGeoEnts();

  //get geobuilder to build cells with needed properties (connection, ..)
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoCellBuilder = getMethodData().getStdTrsGeoBuilder();
  CFAUTOTRACE;
  //get structure of data from geobuilder (cells with properties)
  StdTrsGeoBuilder::GeoData& geoCellData = geoCellBuilder->getDataGE();
  //set that we use only data of inner cells trs
  geoCellData.trs = cells;
  std::vector< RealVector > mappedCoords(0);
  //loop over all inner cells
  techStatesToNodes.resize(nbStates);

  RealVector vertice;
  vertice.resize(nbDim);

  const CFuint nbNodes = nodes.size();
  techNodesCoordinates.resize(nbNodes);
  techNodesToStates.resize(nbNodes);
  for(CFuint iNode=0; iNode < nbNodes; ++iNode)
  {
    techNodesCoordinates[iNode].resize(nbDim);
    techNodesCoordinates[iNode] = *nodes[iNode];
  }

  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
  {
    CFLogDebugMax("Cell " << iGeoEnt << "\n");

    // build the GeometricEntity (cell)
    //set index of cell
    geoCellData.idx = iGeoEnt;
    //geo builder make cell
    GeometricEntity& cell = *geoCellBuilder->buildGE();
    //read vector of states on element (cell)
    std::vector<State*>& cellStates = *(cell.getStates());
    //read number of states on element (cell)
    const CFuint nbStatesInCell = cellStates.size();
    //read vector of nodes on element (cell)
    // unused //    std::vector<Node*>&  nodes  = *cell.getNodes();

    if (mappedCoords.size() < nbStatesInCell)
    {
      CFuint i = mappedCoords.size();
      mappedCoords.resize(nbStatesInCell);
      for (;i<nbStatesInCell;++i)
      {
        mappedCoords[i].resize(nbDim);
      }
    }

    cell.getStatesMappedCoordinates(mappedCoords);
    for(CFuint iState=0; iState < nbStatesInCell; ++iState)
    {
      const CFuint stateID = cellStates[iState]->getLocalID();
//       if (nbDim == 2)
//       {
//         vertice = mappedCoords[iState][0]*(*(nodes[0]))
//                 + mappedCoords[iState][1]*(*(nodes[1]))
//                 + (1-mappedCoords[iState][0]-mappedCoords[iState][1])*(*(nodes[2]));
//       }
//       else
//       {
//         vertice = mappedCoords[iState][0]*(*(nodes[0]))
//                 + mappedCoords[iState][1]*(*(nodes[1]))
//                 + mappedCoords[iState][2]*(*(nodes[2]))
//                 + (1-mappedCoords[iState][0]-mappedCoords[iState][1]-mappedCoords[iState][2])*(*(nodes[3]));
//       }

      vertice = states[stateID]->getCoordinates();

      CFuint iNode=0;
      CFuint nbNodes = techNodesCoordinates.size();
      bool find = false;

      while((iNode < nbNodes)&&(!find))
      {
        RealVector hlpVector = techNodesCoordinates[iNode];
        hlpVector-= vertice;
        if ( (hlpVector.norm1()) > 1.0E-5)
        {
          ++iNode;
        }
        else
        {
          find = true;
        }
      }
      techStatesToNodes[stateID] = iNode;
      if (find)
      {
        techNodesToStates[iNode].push_back(stateID);
      }
      else
      {
        techNodesToStates.resize(iNode+1);
        techNodesToStates[iNode].push_back(stateID);
        techNodesCoordinates.resize(iNode+1);
        techNodesCoordinates[iNode].resize(nbDim);
        techNodesCoordinates[iNode] = vertice;
      }
    }
    //release the GeometricEntity
    geoCellBuilder->releaseGE();
  }

//   for(CFuint iNode=0;iNode < nbNodes;++iNode)
//   {
// // if (techNodesToStates[iNode].size() < 3)
// // {
//     CFout << techNodesCoordinates[iNode] << "  /   " << CFendl;
//     for(CFuint iState=0; iState < techNodesToStates[iNode].size();++iState)
//     {
//     CFout << *states[techNodesToStates[iNode][iState]] << " " << CFendl;
//     }
//     CFout << "\n" << CFendl;
// // }
//   }
  //set integrationIndexes and normals to face
  Common::SafePtr<GeometricEntityPool<Framework::FaceToCellGEBuilder> > geoBuilder = getMethodData().getFaceBuilder();

  // get InnerCells TopologicalRegionSet
//   Common::SafePtr<TopologicalRegionSet>
  cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the face builder and set the TRSs
  Framework::FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = faces;
  geoData.isBoundary = false;
  //loop over all inner faces
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    normals[iFace].resize(nbDim);
    indexes[iFace].resize(4);
    CFLogDebugMax("Face " << iFace << "\n");
    //set index of face
    geoData.idx = iFace;
    //geo builder make face
    GeometricEntity& face = *geoBuilder->buildGE();

    //nodes of the actual face
    std::vector<Node*>&  nodes  = *face.getNodes();
    //get the left (neighbouring) cell
    GeometricEntity* cellLeft  = face.getNeighborGeo(0);

    //get nodes and states of the left cell
    std::vector<Node*>&  left_cell_nodes  = *cellLeft->getNodes();
//     std::vector<State*>& left_cell_states = *cellLeft->getStates();

    //get the right (neighbouring) cell
    GeometricEntity* cellRight = face.getNeighborGeo(1);
    //get nodes and states of the right cell
    std::vector<Node*>&  right_cell_nodes  = *cellRight->getNodes();
//     std::vector<State*>& right_cell_states = *cellRight->getStates();

    CFreal detJacobi;
    if (nbDim == 2)
    {
      RealVector hlp = *nodes[1] - *nodes[0];
      detJacobi = sqrt(hlp[0]*hlp[0]+hlp[1]*hlp[1]);
    }
    else
    {
      //     detJacobi = face.computeVolume()*2;
      RealVector hlp1 = *nodes[1] - *nodes[0];
      RealVector hlp2 = *nodes[2] - *nodes[0];
      detJacobi = sqrt((hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1])*(hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1]) + (hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])*(hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])+(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0])*(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0]));//*2.0/2.0;
    }

//*****************************************************************
//*****************************************************************
    //SET Integrator

    //numbers of quadrature points
    CFuint m_nbKvadrPoint = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getIntegratorPattern()[0];

    //get coordinates of quadrature points
    const std::vector<RealVector>& leftCoord = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getQuadraturePointsCoordinates();

    //get coordinates of quadrature points
    const std::vector<RealVector>& rightCoord = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellRight)->getQuadraturePointsCoordinates();

    //compute gradient of shape functions in quadrature points
    std::vector<RealMatrix> leftGradient = cellLeft->computeSolutionShapeFunctionGradientsInMappedCoordinates(leftCoord);


    //compute gradient of shape functions in quadrature points
    std::vector<RealMatrix> rightGradient = cellRight->computeSolutionShapeFunctionGradientsInMappedCoordinates(rightCoord);

//*****************************************************************
//*****************************************************************

    //find local index of face in left and right element
    CFuint m_idxFaceFromLeftCell= 10;
    if (nbDim == 3)
    {
      for (CFuint i=0; i < 4; i++)
        for (CFuint j=0; j < 3; j++)
          if(((*left_cell_nodes[(i)%4]==*nodes[(j)%3])&&(*left_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*left_cell_nodes[(i+2)%4]==*nodes[(j+2)%3]))||((*left_cell_nodes[(i)%4]==*nodes[(j+2)%3])&&(*left_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*left_cell_nodes[(i+2)%4]==*nodes[(j)%3]))) m_idxFaceFromLeftCell=i;
    }
    else
    {
      for (CFuint i=0; i < 3; i++)
        for (CFuint j=0; j < 2; j++)
          if((*left_cell_nodes[(i)%3]==*nodes[(j)%2])&&(*left_cell_nodes[(i+1)%3]==*nodes[(j+1)%2]))
          {
            m_idxFaceFromLeftCell=i;
          }
    }
    cf_assert (m_idxFaceFromLeftCell!=10);
    indexes[iFace][0]=m_idxFaceFromLeftCell;
    CFuint m_idxFaceFromRightCell= 10;
    if (nbDim == 3)
    {
      for (CFuint i=0; i < 4; i++)
        for (CFuint j=0; j < 3; j++)
          if(((*right_cell_nodes[(i)%4]==*nodes[(j)%3])&&(*right_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*right_cell_nodes[(i+2)%4]==*nodes[(j+2)%3]))||
             ((*right_cell_nodes[(i)%4]==*nodes[(j+2)%3])&&(*right_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*right_cell_nodes[(i+2)%4]==*nodes[(j)%3]))) m_idxFaceFromRightCell=i;
    }
    else
    {
      for (CFuint i=0; i < 3; i++)
        for (CFuint j=0; j < 2; j++)
          if((*right_cell_nodes[(i)%3]==*nodes[(j)%2])&&(*right_cell_nodes[(i+1)%3]==*nodes[(j+1)%2]))
          {
            m_idxFaceFromRightCell=i;
          }
    }

if (m_idxFaceFromRightCell==10)
{
  CFout << "Error in face " << iFace << "\n" << CFendl;
  for (CFuint i=0; i < 3; i++)
    CFout << *right_cell_nodes[i] << " " << CFendl;
  CFout << "\n" << CFendl;
  for (CFuint j=0; j < 2; j++)
    CFout << *nodes[j] << " " << CFendl;
  CFout << "\n" << CFendl;
}

    cf_assert (m_idxFaceFromRightCell!=10);
    indexes[iFace][1]=m_idxFaceFromRightCell;

    if (nbDim == 3)
    {
      CFuint rightHlpIndex=10;
      for (CFuint i=0; i < 3; i++)
      {
        RealVector hlpVector=cellLeft->computeCoordFromMappedCoord(leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint + 1]) - cellRight->computeCoordFromMappedCoord(rightCoord[m_idxFaceFromRightCell*m_nbKvadrPoint + 1 + (2-(1+i)%3)]);
        if (hlpVector.norm1()/detJacobi < 0.0001) rightHlpIndex = i;
      }
      assert(rightHlpIndex!=10);
      indexes[iFace][2]=rightHlpIndex;

      CFuint swifted=0;
      RealVector hlpVector=cellLeft->computeCoordFromMappedCoord(leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint + 2]) - cellRight->computeCoordFromMappedCoord(rightCoord[m_idxFaceFromRightCell*m_nbKvadrPoint + 1 + (2-(1+rightHlpIndex + 1)%3)]);
      if (hlpVector.norm1()/detJacobi > 0.0001) swifted = 1;
      indexes[iFace][3]=swifted;
    }
    if (nbDim == 2)
    {
      normals[iFace][0]= (*left_cell_nodes[(m_idxFaceFromLeftCell+1) % 3])[YY] - (*left_cell_nodes[m_idxFaceFromLeftCell])[YY];
      normals[iFace][1]= (*left_cell_nodes[m_idxFaceFromLeftCell])[XX] - (*left_cell_nodes[(m_idxFaceFromLeftCell+1) % 3])[XX];
      normals[iFace] /= sqrt(pow( normals[iFace][0],2)+pow(normals[iFace][1],2));
    }
    else
    {
      RealVector normal = face.computeAvgCellNormal();
      RealVector HlpNormal = *left_cell_nodes[(m_idxFaceFromLeftCell+3)%4] - *left_cell_nodes[(m_idxFaceFromLeftCell)%4];
      if ((normal[0]*HlpNormal[0]+normal[1]*HlpNormal[1]+normal[2]*HlpNormal[2]) > 0)
      {
        normal *=-1;
      }
      normals[iFace]=normal/sqrt(pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2));
    }

    // release the face
    geoBuilder->releaseGE();
  }
// Add integration indexes for boundary faces
  std::vector< Common::SafePtr<TopologicalRegionSet> > m_TrsList = MeshDataStack::getActive()->getTrsList();
  const CFuint number_TRS = m_TrsList.size();
  Framework::GeometricEntity * m_face;

  for(CFuint iTRS = 0; iTRS < number_TRS; ++iTRS)
  {
    if ((*m_TrsList[iTRS]).hasTag("boundary"))
    {
      Common::SafePtr<TopologicalRegionSet> bFaces = m_TrsList[iTRS];//->getCurrentTRS();

      // get the geodata of the face builder and set the TRSs
      geoData.cellsTRS = cells;
      geoData.facesTRS = bFaces;
      geoData.isBoundary = true;

      //get number of faces
      const CFuint nbFaces = bFaces->getLocalNbGeoEnts();
//          CFout << "\n  " << nbFaces << CFendl;
      //loop over all boundary faces
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
      {
        CFLogDebugMax("Face " << iFace << "\n");
        //set index of face
        geoData.idx = iFace;
        //geo builder make face
        m_face = geoBuilder->buildGE();
        //get the (neighbouring) cell
        GeometricEntity* cellLeft  = m_face->getNeighborGeo(0);
        //get nodes of face
        std::vector<Node*>&  nodes  = *m_face->getNodes();
        //get nodes and states of the cell
        std::vector<Node*>&  left_cell_nodes  = *cellLeft->getNodes();
        std::vector<State*>& left_cell_states = *cellLeft->getStates();
        m_face->resizeStates(nodes.size());
        for(CFuint iNode = 0; iNode < nodes.size(); ++iNode)
        {
          CFuint iNodeCell = 0;
          while (left_cell_nodes[iNodeCell] != nodes[iNode])
          {
            ++iNodeCell;
          }
          CFuint iStateCell = 0;
          while (left_cell_states[iStateCell]->getCoordinates() != *nodes[iNode])
          {
            ++iStateCell;
          }
          m_face->setState(iNode, left_cell_states[iStateCell]);
//           CFout << "\n" <<  *left_cell_nodes[iNodeCell] << "  " << *nodes[iNode] << " " << *left_cell_states[iStateCell] << CFendl;
        }

//          std::vector<Node*>&  nodes  = *m_face->getNodes();
// CFout << "  " << nodes.size() << CFendl;
// CFout << "ahoj  " << CFendl;
// if (states[0]->hasLocalID())
// {
//          CFout << "ahoj2  " << CFendl;
// }
//          CFuint hlp = states[0]->getGlobalID();
//          CFout << "ahoj  " << CFendl;
// CFout << "  " << hlp << CFendl;
        geoBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

