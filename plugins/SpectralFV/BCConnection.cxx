#include "Environment/DirPaths.hh"

#include "Framework/BadFormatException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/BaseBndFaceTermComputer.hh"
#include "SpectralFV/BCConnection.hh"
#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SpectralFVElementData.hh"

#include "Common/ConnectivityTable.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider< BCConnection,
                                   SpectralFVMethodData,
                                   BCStateComputer,
                                   SpectralFVModule
                                 > BCConnectionProvider("Connection");

//////////////////////////////////////////////////////////////////////////////

void BCConnection::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector< CFuint > >("ConnectedNodes","List of nodes that are connected (two consecutive nodes are connected)");
  options.addConfigOption< std::string >("ConnNodesFileName","Name of file containing the connected node IDs");
}

//////////////////////////////////////////////////////////////////////////////

BCConnection::BCConnection(const std::string& name) :
  BCStateComputer(name),
  m_states(CFNULL),
  m_cellStatesConn(CFNULL),
  m_statesReconstr(CFNULL),
  m_connectedNodes(),
  m_faceToGCellAndFaceOrient(),
  m_svFaceflxPntsRecCoefs(CFNULL),
  m_ghostCellStates(),
  m_orient(),
  m_connNodesFileName()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_connectedNodes = vector<CFuint>();
  setParameter("ConnectedNodes",&m_connectedNodes);

  m_connNodesFileName = "";
  setParameter("ConnNodesFileName",&m_connNodesFileName);
}

//////////////////////////////////////////////////////////////////////////////

BCConnection::~BCConnection()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCConnection::computeGhostStates(const vector< State* >& intStates,
                                       vector< State* >& ghostStates,
                                       const std::vector< RealVector >& normals,
                                       const std::vector< RealVector >& coords)
{
  m_statesReconstr->reconstructStates(m_ghostCellStates,ghostStates,
                                      (*m_svFaceflxPntsRecCoefs)[m_orient][RIGHT],
                                      m_ghostCellStates.size());
}

//////////////////////////////////////////////////////////////////////////////

void BCConnection::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                          std::vector< std::vector< RealVector* > >& ghostGrads,
                                          const std::vector< RealVector >& normals,
                                          const std::vector< RealVector >& coords)
{
  throw Common::NotImplementedException
      (FromHere(), "Connection BC for the gradients not implemented (add gradients datahandle and set ghost gradients for this).");
  m_statesReconstr->reconstructGradients(m_ghostCellGrads,ghostGrads,
                                         (*m_svFaceflxPntsRecCoefs)[m_orient][RIGHT],
                                         m_ghostCellGrads.size());
}

//////////////////////////////////////////////////////////////////////////////

void BCConnection::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get datahandles to states
  m_states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();

  // get cells-states connectivity
  m_cellStatesConn = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  // get the states reconstructor
  m_statesReconstr = getMethodData().getStatesReconstructor();

  // Get solution reconstruction coefficients
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  std::string computerType = getMethodData().getBndFaceTermComputer()->getType();
  if (computerType == "Std")
  {
    m_svFaceflxPntsRecCoefs = svLocalData[0]->getExtQPntPolyValsPerOrientNoSymm();
  }
  if (computerType == "QuadFree")
  {
    m_svFaceflxPntsRecCoefs = svLocalData[0]->getSolInFaceFluxPntCoefPerOrientNoSymm();
  }
  cf_assert(m_svFaceflxPntsRecCoefs.isNotNull());

  if (m_connectedNodes.size() == 0)
  {
    const std::string fileName =
        Environment::DirPaths::getInstance().getWorkingDir().string() + "/"
        + m_connNodesFileName;
    ifstream connNodesFile;
    connNodesFile.open(fileName.c_str(), ios::in);
    if (connNodesFile.is_open())
    {
      CFuint nbrNodes;
      connNodesFile >> nbrNodes;

      m_connectedNodes.resize(2*nbrNodes);
      for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
      {
        connNodesFile >> m_connectedNodes[2*iNode  ];
        connNodesFile >> m_connectedNodes[2*iNode+1];
      }
    }
    connNodesFile.close();
  }
  if (m_connectedNodes.size() == 0)
  {
    throw BadFormatException (FromHere(),"Connected nodes not specified for connection boundary condition");
  }

  // create the face to cell and face orientation data
  createFaceConnectionData();
}

//////////////////////////////////////////////////////////////////////////////

void BCConnection::createFaceConnectionData()
{
  CFAUTOTRACE;

  // get TRS list
  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();

  // get boundary TRSs
  if (m_trsNames.size() != 2)
  {
    throw BadFormatException (FromHere(),"Connection BC with different number of TRSs than 2");
  }
  vector< SafePtr< TopologicalRegionSet > > m_connTRSs(2);
  const CFuint nbTRSs = trsList.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    if (trsList[iTRS]->getName() == m_trsNames[0])
    {
      if (m_connTRSs[0] == CFNULL)
      {
        m_connTRSs[0] = trsList[iTRS];
      }
      else
      {
        throw BadFormatException (FromHere(),"Two TRSs with the same name found!");
      }
    }
    if (trsList[iTRS]->getName() == m_trsNames[1])
    {
      if (m_connTRSs[1] == CFNULL)
      {
        m_connTRSs[1] = trsList[iTRS];
      }
      else
      {
        throw BadFormatException (FromHere(),"Two TRSs with the same name found!");
      }
    }
  }
  if (m_connTRSs[0] == CFNULL || m_connTRSs[1] == CFNULL)
  {
    throw BadFormatException (FromHere(),"Boundary TRS not found!");
  }

  // get cell-node connectivity
  SafePtr< ConnectivityTable< CFuint > > cellNodeConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get boundary face connectivities
  vector< SafePtr< ConnectivityTable< CFuint > > > bndFaceCellConn (2);
  vector< SafePtr< ConnectivityTable< CFuint > > > bndFaceNodesConn(2);
  vector< SafePtr< vector< CFuint > > >            bndFaceLocalIDs (2);
  for (CFuint iTRS = 0; iTRS < 2; ++iTRS)
  {
    bndFaceCellConn [iTRS] = MeshDataStack::getActive()->getConnectivity(m_trsNames[iTRS]+"-Faces2Cells");
    bndFaceNodesConn[iTRS] = MeshDataStack::getActive()->getConnectivity(m_trsNames[iTRS]+"Nodes");
    bndFaceLocalIDs [iTRS] = m_connTRSs[iTRS]->getGeoEntsLocalIdx();
/*    const CFuint nbrFaces = bndFaceNodesConn[iTRS]->nbRows();
    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      const CFuint nbrFaceNodes = bndFaceNodesConn[iTRS]->nbCols(iFace);
      for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
      {
        CF_DEBUG_OBJ((*bndFaceNodesConn[iTRS])(iFace,iNode));
      }
    }*/
  }

  // get number of faces in the TRSs
  const CFuint nbrFaces = bndFaceCellConn[0]->nbRows();
  if (nbrFaces != bndFaceCellConn [1]->nbRows() ||
      nbrFaces != bndFaceNodesConn[0]->nbRows() ||
      nbrFaces != bndFaceNodesConn[1]->nbRows() ||
      nbrFaces != bndFaceLocalIDs [0]->size()   ||
      nbrFaces != bndFaceLocalIDs [1]->size())
  {
    throw BadFormatException (FromHere(),"Connected TRSs have a different number of faces");
  }

  // Get vector containing the nodes connectivities for each orientation
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  vector< vector < vector < CFuint > > >
      nodeConnPerOrient = *svLocalData[0]->getSVFaceNodeConnPerOrientationNoSymm();
  const CFuint nbrOrient = nodeConnPerOrient.size();

  // loop over TRSs to compute face ghost cell IDs and orientations
  for (CFuint iTRS = 0; iTRS < 2; ++iTRS)
  {
    // connected TRS
    const CFuint iConnTRS = iTRS == 0 ? 1 : 0;
//     CF_DEBUG_OBJ(iConnTRS);

    // loop over faces
    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
//       CF_DEBUG_OBJ(iFace);

      // get the face node IDs, and the matching node IDs (of the connected face)
      const CFuint nbrFaceNodes = bndFaceNodesConn[iTRS]->nbCols(iFace);
      vector< CFuint > currFaceNodes(nbrFaceNodes);
      vector< CFuint > connFaceNodes(nbrFaceNodes);
      for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
      {
        currFaceNodes[iNode] = (*bndFaceNodesConn[iTRS])(iFace,iNode);
        connFaceNodes[iNode] = getMatchingNodeID(currFaceNodes[iNode]);
      }

      // get the connected face local index
      CFuint iConnFace = 0;
      bool faceFound = false;
      for (CFuint iFace2 = 0; iFace2 < nbrFaces && !faceFound; ++iFace2)
      {
        faceFound = true;
        for (CFuint iNode = 0; iNode < nbrFaceNodes && faceFound; ++iNode)
        {
          faceFound = false;
          const CFuint nodeID = connFaceNodes[iNode];
          for (CFuint iNode2 = 0; iNode2 < nbrFaceNodes && !faceFound; ++iNode2)
          {
            faceFound = nodeID == (*bndFaceNodesConn[iConnTRS])(iFace2,iNode2);
          }
        }

        if (faceFound)
        {
          iConnFace = iFace2;
        }
      }
      if (!faceFound)
      {
        throw BadFormatException (FromHere(),"Connected face not found");
      }

      // get internal cell and ghost cell IDs
      const CFuint intCellID   = (*bndFaceCellConn[iTRS    ])(iFace    ,0);
      const CFuint ghostCellID = (*bndFaceCellConn[iConnTRS])(iConnFace,0);

      // set ghost cell ID
      const CFuint faceID = (*bndFaceLocalIDs[iTRS])[iFace];
      m_faceToGCellAndFaceOrient[faceID].first = ghostCellID;

      // get cell node IDs
      const CFuint nbrIntCellNodes   = cellNodeConn->nbCols(intCellID  );
      vector< CFuint > intCellNodes  (nbrIntCellNodes  );
      for (CFuint iNode = 0; iNode < nbrIntCellNodes  ; ++iNode)
      {
        intCellNodes  [iNode] = (*cellNodeConn)(intCellID  ,iNode);
      }
      const CFuint nbrGhostCellNodes = cellNodeConn->nbCols(ghostCellID);
      vector< CFuint > ghostCellNodes(nbrGhostCellNodes);
      for (CFuint iNode = 0; iNode < nbrGhostCellNodes; ++iNode)
      {
        ghostCellNodes[iNode] = (*cellNodeConn)(ghostCellID,iNode);
      }

      // determine face orientation
      bool orientFound = false;
      for (CFuint iOrient = 0; iOrient < nbrOrient && !orientFound; ++iOrient)
      {
        // get cellFaceNodes corresponding to this orientation
        // (assuming all faces have the same number of nodes here!!!)
        vector< CFuint > intCellFaceNodes  (nbrFaceNodes);
        for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
        {
          const CFuint localIdx = nodeConnPerOrient[iOrient][LEFT ][iNode];
          intCellFaceNodes  [iNode] = intCellNodes  [localIdx];
        }
        vector< CFuint > ghostCellFaceNodes(nbrFaceNodes);
        for (CFuint iNode = 0; iNode < nbrFaceNodes; ++iNode)
        {
          const CFuint localIdx = nodeConnPerOrient[iOrient][RIGHT][iNode];
          ghostCellFaceNodes[iNode] = ghostCellNodes[localIdx];
        }

        // check whether nodes correspond to current faces nodes
        orientFound = true;
        for (CFuint iNode = 0; iNode < nbrFaceNodes && orientFound; ++iNode)
        {
          const CFuint nodeID = currFaceNodes[iNode];
          orientFound = false;
          for (CFuint iNode2 = 0; iNode2 < nbrFaceNodes && !orientFound; ++iNode2)
          {
            orientFound = nodeID == intCellFaceNodes  [iNode2];
          }
        }
        for (CFuint iNode = 0; iNode < nbrFaceNodes && orientFound; ++iNode)
        {
          const CFuint nodeID = connFaceNodes[iNode];
          orientFound = false;
          for (CFuint iNode2 = 0; iNode2 < nbrFaceNodes && !orientFound; ++iNode2)
          {
            orientFound = nodeID == ghostCellFaceNodes[iNode2];
          }
        }

        // check if the face nodes coincide
        for (CFuint iNode = 0; iNode < nbrFaceNodes && orientFound; ++iNode)
        {
          orientFound =
              getMatchingNodeID(intCellFaceNodes[iNode]) == ghostCellFaceNodes[iNode];
        }

        if (orientFound)
        {
          m_faceToGCellAndFaceOrient[faceID].second = iOrient;
//           CF_DEBUG_OBJ(iOrient);
        }
      }
      if (!orientFound)
      {
        throw BadFormatException (FromHere(),"Face orientation not found");
      }
    }
  }
//   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////////////

CFuint BCConnection::getMatchingNodeID(const CFuint givenNode)
{
  // number of nodes in the connected node list
  const CFuint nbrConnNodes = m_connectedNodes.size();
  cf_assert(nbrConnNodes % 2 == 0);

  // search matching node
  for (CFuint iNode = 0; iNode < nbrConnNodes; iNode += 2)
  {
    if (givenNode == m_connectedNodes[iNode])
    {
      return m_connectedNodes[iNode+1];
    }
    if (givenNode == m_connectedNodes[iNode+1])
    {
      return m_connectedNodes[iNode];
    }
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void BCConnection::setBndFaceData()
{
  // set ghost cell states
  const CFuint cellID = m_faceToGCellAndFaceOrient[m_faceID].first;
  const CFuint nbrStates = m_cellStatesConn->nbCols(cellID);
  m_ghostCellStates.resize(nbrStates);
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStatesConn)(cellID,iState);
    m_ghostCellStates[iState] = m_states[stateID];
  }

  // set orientation for this face
  m_orient = m_faceToGCellAndFaceOrient[m_faceID].second;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
