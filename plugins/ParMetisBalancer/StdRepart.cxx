#include <boost/filesystem/fstream.hpp>

#include "Common/PE.hh"
#include "Common/NoSuchValueException.hh"
#include "Common/BadValueException.hh"
#include "Common/MPI/MPIException.hh"
#include "Common/ProcessInfo.hh"
#include "Common/CFLog.hh"
#include "Common/OSystem.hh"

#include "Environment/ConcreteProvider.hh"
#include "Environment/DirPaths.hh"


#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodData.hh"
#include "Framework/PathAppender.hh"
#include "Framework/MeshCreator.hh"

#include "ParMetisBalancer/ParMetisBalancer.hh"
#include "ParMetisBalancer/ParMetisBalancerModule.hh"
#include "ParMetisBalancer/StdRepart.hh"


#include "parmetis.h"
#include <mpi.h>
//#include <algorithm>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace ParMetisBalancer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdRepart, ParMetisBalancerData, ParMetisBalancerModule>
stdParMetisBalancerProvider("StdRepart");

//////////////////////////////////////////////////////////////////////////////

StdRepart::StdRepart(const std::string& name) :
  ParMetisBalancerCom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  m_cells(NULL),
  m_nodes(NULL),
  m_states(NULL)
{
  /// Inicializes the command "StdRepart" and sets data socets to be used
}

//////////////////////////////////////////////////////////////////////////////

StdRepart::~StdRepart()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdRepart::configure(Config::ConfigArgs& args)
{
  ParMetisBalancerCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void StdRepart::execute()
{
  // Set handles to Cells, Nodes, States
  m_cells = MeshDataStack::getActive()->getTrs("InnerCells");
  m_nodes = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  m_states= socket_states.getDataHandle();
  // Get info on send/recive nodes, apply part1 coloring
  setupDataStorage();
  // Create continus glogal mapping (requierd by PARMetis)
  globalNumNode();
  // Builds conectivity in CSR format (requierd by PARMetis)
  setupCSR();
  //Parmetis is called, new partitioning is calculated
  callParMetisAdaptiveRepart();
  // update interface data on part2 coloring
  UpdateInterfacePart2();
  // select cells to be send (negotiate cell ownership)
  SelectCellsToSend();
  // select nodes to send
  SelectNodesToSend();
  // build MPI send/recive structures
  PrepereMPIcommStruct();
  // perform MPI communication session
  MPICommunicate();
  // Prepere data for mesh update
  PrepereToUpdate();

  // Raport the results
  CFLogInfo("TecPlotFile write \n");
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  boost::filesystem::path path = "./FSOmesh/bal_test_interf.dat"; // Storage for testing purposes only
  if(dim == 2) DoWriteTec<2>(path,true);
  if(dim == 3) DoWriteTec<3>(path,true);
  
  path = "./FSOmesh/bal_test_nointerf.dat"; // Storage for testing purposes only
  if(dim == 2) DoWriteTec<2>(path,false);
  if(dim == 3) DoWriteTec<3>(path,false);
  
  path="./FSOmesh/bal_test_noremoved.dat"; // Storage for testing purposes only
  if(dim == 2) DoWriteTecNoRemoved<2>(path);
  if(dim == 3) DoWriteTecNoRemoved<3>(path);
  
  path="./FSOmesh/bal_test_recived.dat"; // Storage for testing purposes only
  if(dim == 2) DoWriteTecAfterSendRecive<2>(path);
  if(dim == 3) DoWriteTecAfterSendRecive<3>(path);
  
  // free the alocated memory
  DoClearMemory();
}

//////////////////////////////////////////////////////////////////////////////
void StdRepart::setupDataStorage()
{
  // Set corect sizes for data containers
  dataStorage.ResizeData( m_cells->getNbNodesInTrs() );
  // CommonNodes and part1 coloring are set here
  dataStorage.setupCommonNodes( m_nodes.getGhostSendList(), m_nodes.getGhostReceiveList() );
}

//////////////////////////////////////////////////////////////////////////////
void StdRepart::globalNumNode()
{
  int myNodes = m_nodes.size(); // total number of known m_nodes
  for(CFuint i=0; i<m_nodes.getGhostReceiveList().size(); ++i)
  {
    myNodes -= m_nodes.getGhostReceiveList()[i].size();        // each node that is send from outside is foreign
  }
  
  // MPI Send / Recive data struct
  vector<int> nodenum(PE::GetPE().GetProcessorCount(), myNodes);
  vector<int> in(PE::GetPE().GetProcessorCount(), 0);

  // MPI comunication...
  MPI_Alltoall(&nodenum[0], 1,MPI_INT,&in[0],1,MPI_INT,MPI_COMM_WORLD);
  // ... Done, record.
  dataStorage.vtxdist.resize(PE::GetPE().GetProcessorCount()+1);
  dataStorage.vtxdist[0]=0;
  dataStorage.xadj.resize(myNodes+1);
  dataStorage.xadj[0] = 0;
  for(CFuint i=1,n1=0; i<PE::GetPE().GetProcessorCount()+1; ++i) 
  {
    n1+=in[i-1];
    dataStorage.vtxdist[i]=n1;
  }
  // Assigning global numbers to my m_nodes...
  CFuint tmp_num = 0;
  for(CFuint i=0; i<m_nodes.size(); ++i) if(dataStorage.Part1()[i] == PE::GetPE().GetRank())
  {
    dataStorage.GlobId()[i] = dataStorage.vtxdist[PE::GetPE().GetRank()] + tmp_num;
    ++tmp_num;
  }
  //MPI session -----------
  commMPIstruct<CFuint> myGlobalId; // MPI data struct
  // feed data
  for(CFuint i=0; i<PE::GetPE().GetProcessorCount(); ++i)
    for(CFuint j=0; j<m_nodes.getGhostSendList()[i].size(); ++j)
      myGlobalId.AddSendElem(i, dataStorage.GlobId()[ m_nodes.getGhostSendList()[i][j] ]);
  //Set recive sizes
  //myGlobalId.SetSizeRecivBuf() // -- the recive size is known -> communication unnecessary
  for(CFuint i=0; i<PE::GetPE().GetProcessorCount(); ++i)
    myGlobalId.GetReciData()[i].resize( m_nodes.getGhostReceiveList()[i].size() );
  // MPI communication...
  myGlobalId.MPIcommunicate();
  
  // Apply global names...
  for(CFuint i=0; i<PE::GetPE().GetProcessorCount(); ++i)
    for(CFuint j=0; j<myGlobalId.GetReciData()[i].size(); ++j)
      dataStorage.GlobId()[ m_nodes.getGhostReceiveList()[i][j] ] = myGlobalId.GetReciData()[i][j]; // global names of send m_nodes
  //myGlobalId.MPIraport(); // will cause a lot of output :)
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void StdRepart::setupCSR()
{
  vector<vector<CFuint> > tmp_kon(m_nodes.size()); //for storing conectivity of each node, i'th node has a list o neighbours
  
  // Loop through all cells...
  for(CFuint icell=0; icell < m_cells->getLocalNbGeoEnts(); ++icell) // through all m_cells
  {
    for(CFuint inode=0; inode < m_cells->getNbNodesInGeo(icell); ++inode) // through m_nodes of a cell
    {
      if(dataStorage.Part1()[m_cells->getNodeID(icell,inode)] != PE::GetPE().GetRank())
        continue; // if node is not owned by the proces, than skip it...
      for(CFuint jnode=0; jnode < m_cells->getNbNodesInGeo(icell); ++jnode)
        if(jnode != inode) // only the other m_nodes
        {
          tmp_kon[ m_cells->getNodeID(icell,inode) ].push_back( dataStorage.GlobId()[ m_cells->getNodeID(icell,jnode) ] );
        }
    }
  }
  // need to sort and unique the tmp_kon sub vectors (make them sets)
  vector<CFuint>::iterator it;
  for(CFuint i=0; i<tmp_kon.size(); ++i)
  {
    sort(tmp_kon[i].begin(), tmp_kon[i].end());
    it = unique (tmp_kon[i].begin(), tmp_kon[i].end());
    tmp_kon[i].resize( it - tmp_kon[i].begin() );
  }
  
  /// record into CSR...
  vector<PartitionerData::IndexT>::iterator adj_iter = dataStorage.adjncy.begin();
  vector<PartitionerData::IndexT>::iterator xadj_iter= dataStorage.xadj.begin();
  for(CFuint i=0; i < m_nodes.size(); ++i)
  {
    dataStorage.adjncy.insert( adj_iter, tmp_kon[i].begin(), tmp_kon[i].end() );
    adj_iter = dataStorage.adjncy.end();
    if(tmp_kon[i].size()>0)
    {
      *(xadj_iter+1) = *(xadj_iter) + tmp_kon[i].size();
      xadj_iter++;
      if(xadj_iter == dataStorage.xadj.end())
        throw Common::MPIException (FromHere(),  "Problem z tworzeniem CSR'a - allocation error");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
void StdRepart::callParMetisAdaptiveRepart()
{
  CFAUTOTRACE;
  
  // number of owned m_nodes
  CFuint myNodes = dataStorage.vtxdist[PE::GetPE().GetRank()+1] - dataStorage.vtxdist[PE::GetPE().GetRank()];

  MPI_Comm comm = PE::GetPE().GetCommunicator();   // Do wywolania ParMETIS'a
  
  // ParMETIS flags and options  -- see ParMETIS documentation
  PartitionerData::IndexT wgtflag = 0;                 // 0 for no weights(vwgt and adjwgt=NULL), 2 weight on vertices only(adjwgt=NULL)
  PartitionerData::IndexT numflag = 0;                 // 0 for C style array numbering
  PartitionerData::IndexT ncon    = 0;                 // no of weights for each vertex
  PartitionerData::IndexT nparts  = PE::GetPE().GetProcessorCount(); // number of subdomains - processes
  
  PartitionerData::RealT itr = 1000;
  
  PartitionerData::IndexT options[4];
  options[0] = 0;
  
  PartitionerData::IndexT edgecut =0;                  //number of egdes that are cutted

  PartitionerData::RealT *tpwgts = NULL;            //array for weights
  PartitionerData::RealT *ubvec = NULL;             //array of size ncon for imbalance tolerance

  PartitionerData::IndexT  *vwgt=NULL;
  PartitionerData::IndexT  *adjwgt=NULL;
  PartitionerData::IndexT  *vsize=NULL; //weight arrays
  PartitionerData::IndexT* part = new PartitionerData::IndexT[myNodes];
  for(CFuint i=0; i<myNodes; ++i)
  {
    part[i] = 0;
  }
  
  ///HACK:: introduce waights for testing -- to be removed later!! 
    wgtflag=2;  // 0 for no weights(vwgt and adjwgt=NULL), 2 weight on vertices only(adjwgt=NULL)
    ncon   =1;  // no of weights for each vertex
    
    tpwgts = new PartitionerData::RealT[ncon*nparts];
    for(int i=0; i<(ncon*nparts); ++i) tpwgts[i] = 1./nparts;
    
    ubvec = new PartitionerData::RealT[ncon];
    for(int i=0; i<ncon; ++i) ubvec[i] = 1.05;
    
    vwgt = new PartitionerData::IndexT[myNodes];
    for(CFuint i=0; i<(myNodes); ++i) vwgt[i]=1;

      if(PE::GetPE().GetRank()==0)
        for(CFuint i=0; i<(myNodes); ++i) vwgt[i]=3;
      if(PE::GetPE().GetRank()==1)
        for(CFuint i=0; i<(myNodes); ++i) vwgt[i]=10;
      if(PE::GetPE().GetRank()==2)
        for(CFuint i=0; i<(myNodes); ++i) vwgt[i]=5;
  /// /////////////////////////////////////////////////////////////

  CFLogDebugMin( "Calling ParMetis::AdaptiveRepart()\n");
  Common::Stopwatch<Common::WallTime> MetisTimer;

  ParMETIS_V3_AdaptiveRepart (&dataStorage.vtxdist[0], &dataStorage.xadj[0], &dataStorage.adjncy[0], vwgt, vsize, adjwgt, 
        &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, &itr, options, &edgecut, part, &comm);
  
  CFuint i1=0;
  for(CFuint i=0; i<m_nodes.size(); ++i) if( dataStorage.Part1()[i] == PE::GetPE().GetRank() )
  {
    //cout<<part[i1]<<" ";
    dataStorage.Part2()[i]=part[i1]; //ustawienie nowego podzialu w wezlach pomijajac obce wezly
    i1+=1;
  }
  //cout<<" myNodes:"<<PE::GetPE().GetRank()<<" "<<myNodes<<" "<<i1<<endl;
  delete [] part;
  delete [] vwgt;
  delete [] tpwgts;
  delete [] ubvec;

  MetisTimer.stop ();
  CFLog(NOTICE, "ParMetis::AdaptiveRepart() took " << MetisTimer << "\n");
}

//////////////////////////////////////////////////////////////////////////////
void StdRepart::UpdateInterfacePart2()
{
  //MPI session
  commMPIstruct<CFuint> interfPartMPI; // MPI data struct
  // feed data
  for(CFuint i=0; i<PE::GetPE().GetProcessorCount(); ++i)
    for(CFuint j=0; j<m_nodes.getGhostSendList()[i].size(); ++j)
      interfPartMPI.AddSendElem(i, dataStorage.Part2()[ m_nodes.getGhostSendList()[i][j] ]);
  //Set recive sizes
  //myGlobalId.SetSizeRecivBuf() // -- the recive size is known -> communication unnecessary
  for(CFuint i=0; i<PE::GetPE().GetProcessorCount(); ++i)
    interfPartMPI.GetReciData()[i].resize( m_nodes.getGhostReceiveList()[i].size() );
  // MPI communication...
  interfPartMPI.MPIcommunicate();
  
  // Apply part2 to interface nodes
  for(CFuint i=0; i<PE::GetPE().GetProcessorCount(); ++i)
    for(CFuint j=0; j<m_nodes.getGhostReceiveList()[i].size(); ++j)
      dataStorage.Part2()[ m_nodes.getGhostReceiveList()[i][j] ] = interfPartMPI.GetReciData()[i][j]; // global names of send m_nodes
}

//////////////////////////////////////////////////////////////////////////////
/// Make "shoping lists" of entities to be send/removed from procesors
void StdRepart::SelectCellsToSend()
{
  /// Processes are only concerned for cells thay own,
  // proces owns a cell which is either composed of nodes with its part, or for which it is the oldest (has lowest rank)
  // I use the following approch to decide whather a cell is to be send/ droped
  // 1. Cells are sent to a process to which all its nodes are signed
  // 2. If nodes composing a cell belong to more than one process, the cell is signed to the oldest one (by rank).
  // 3. If a cell is to be send someware it is necessary to check if the process is not allready in posesion of the cell.
  //    For this I use dataStorage.commNodes container. If all nodes are in a container, I am sure the cell is known.
  //!        it is actually a piece of bull... since there might be a cell which has all nodes in commnodes, and still is not known!!
  // 4. So called triple pionts will be sended anyway, as not all nodes known to the process might be specified in commNodes.
  // 5. It also applies to others cells sometime! To corect for that a process will have to scan through its "ghostCells" list
  //    to avoid repeting cells
  // 6. It should not happen that more than one process sends the same cell, as a process is only allowed to send a cell it owns,
  //    and cell ownership is (I hope :-) well defined.
  /// Above should, assure proper mesh reconstruction!
  
  for(CFuint icell=0; icell < m_cells->getLocalNbGeoEnts(); ++icell) // through all m_cells - decisions are made for each cell
  {
    bool my_to_decide = true; // is process an owner?
    CFuint destination_proces = PE::GetPE().GetProcessorCount(); // to the youngest process, (highest rank number)
    for(CFuint inode=0; inode < m_cells->getNbNodesInGeo(icell); ++inode) // through m_nodes composing a cell
    {
      if( dataStorage.Part1()[m_cells->getNodeID(icell,inode)] < PE::GetPE().GetRank() ) // if there is a node from older process
      {           //  node part                          <          //  processor rank  => the cell is own by someone else!!
        dataStorage.GhostCell().push_back(icell); // record cell as a foreign cell
        my_to_decide = false;                     // mark cell as foreign
        break;                                    // go out of for loop into cell loop
      }
      destination_proces = min(destination_proces, dataStorage.Part2()[m_cells->getNodeID(icell,inode)]); // where to?
    }
    if( my_to_decide && destination_proces != PE::GetPE().GetRank() )
    {
      dataStorage.CellsToSend()[destination_proces].push_back(icell);  // sign to send
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
/// Makes a list of nodes that are to be sended to other processes
void StdRepart::SelectNodesToSend()
{ // Based upon CellsToSend() list, nodes are to be selected to go along
  // It might however happen that a node necessary to create a cell does not belong to a process which sends a cell.
  // Since there is no assurance that the process owning a node will send it as well (the cell which it has might go someware else!),
  // A process sends all nodes necessary for the cell creation, commNodes container is chacked if the target process is aware of the node.
  // It might still happen that a node allready known to a process is sended, or that more than one process sends it,
  // due to that it will be nessecary to correct for that.
  
  for(CFuint iproc=0; iproc<PE::GetPE().GetProcessorCount(); ++iproc) // through all processes
  if(dataStorage.CellsToSend()[iproc].size() > 0)  // only if CellsToSend() not empty
  {
    for(CFuint iscell=0; iscell<dataStorage.CellsToSend()[iproc].size(); ++iscell) // through CellsToSend() collection for iproc
    {
      for(CFuint inode=0; inode < m_cells->getNbNodesInGeo( dataStorage.CellsToSend()[iproc][iscell] ); ++inode) // through nodes of the sended cell
      {
        if( !dataStorage.IsCommon(m_cells->getNodeID(dataStorage.CellsToSend()[iproc][iscell] ,inode), iproc) ) // if the node is not known to the process
          dataStorage.NodesToSend()[iproc].push_back( m_cells->getNodeID(dataStorage.CellsToSend()[iproc][iscell], inode) ); // record for sending
      }
    }
    dataStorage.UniqeSendNodes(iproc);
  }
}

//////////////////////////////////////////////////////////////////////////////
/// Preperes send recive structures for future MPI comuniation, send struct are filled with data
void StdRepart::PrepereMPIcommStruct()
{
  for(CFuint proc=0; proc<PE::GetPE().GetProcessorCount(); ++proc) if(proc != PE::GetPE().GetRank()) // through all but me
  {
    // nodes & states
    for(CFuint inode=0; inode < dataStorage.NodesToSend()[proc].size(); ++inode)
    {
      // four intiger values goes with the node
      dataStorage.m_intData[0].AddSendElem( proc, m_nodes[ dataStorage.NodesToSend()[proc][inode] ]->getGlobalID() );
      dataStorage.m_intData[0].AddSendElem( proc, dataStorage.GlobId()[ dataStorage.NodesToSend()[proc][inode] ] );
      dataStorage.m_intData[0].AddSendElem( proc, dataStorage.Part1()    [ dataStorage.NodesToSend()[proc][inode] ] );
      dataStorage.m_intData[0].AddSendElem( proc, dataStorage.Part2()    [ dataStorage.NodesToSend()[proc][inode] ] );
      // xyz of a node
      for(CFuint j=0; j<(*m_nodes[0]).size(); ++j)
        dataStorage.m_douData[0].AddSendElem( proc, (*m_nodes[ dataStorage.NodesToSend()[proc][inode] ])[j] );
      // data of a state
      for(CFuint j=0; j<(*m_states[0]).size(); ++j)
        dataStorage.m_douData[2].AddSendElem( proc, (*m_states[ dataStorage.NodesToSend()[proc][inode] ])[j] );
    }
    // cells
    for(CFuint icell=0; icell < dataStorage.CellsToSend()[proc].size(); ++icell)
    {
      for(CFuint inode=0; inode < m_cells->getNbNodesInGeo( dataStorage.CellsToSend()[proc][icell] ); ++inode) // through nodes of the sended cell
      {
        dataStorage.m_intData[1].AddSendElem( proc, dataStorage.GlobId()[ m_cells->getNodeID(dataStorage.CellsToSend()[proc][icell], inode) ] );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
/// Processors exchange data
void StdRepart::MPICommunicate()
{
  dataStorage.ExecuteMPIcomm();
}

//////////////////////////////////////////////////////////////////////////////
void StdRepart::PrepereToUpdate()
{
  dataStorage.MakeListsOfKnownNodes();// List of nodes that could be doubled is made
  
  // assume all sendCells for deletion - have been sent and are not mine
  for( CFuint p=0; p<PE::GetPE().GetProcessorCount(); ++p ) if(p != PE::GetPE().GetRank() )
    dataStorage.DelCellCol().insert(dataStorage.DelCellCol().end(), dataStorage.CellsToSend()[p].begin(), dataStorage.CellsToSend()[p].end());
  
  // Do a ghost cell container for lookup, cheack weather a cell is to be removed
  for(CFuint icell=0; icell<dataStorage.GhostCell().size(); ++icell) // through all cells in ghost cell container
  {
    vector<CFuint> singleCell(m_cells->getNbNodesInGeo( dataStorage.GhostCell()[ icell ] ) ); //temporary for nodes
    bool toBeDeleted = true; // is the cell to be removed?
    for(CFuint inode=0; inode < m_cells->getNbNodesInGeo( dataStorage.GhostCell()[ icell ] ); ++inode) // through m_nodes composing a cell
    {
      if( dataStorage.Part2()[m_cells->getNodeID(dataStorage.GhostCell()[icell],inode)] == PE::GetPE().GetRank() ) // if cell has a node that is mine
      {           //  node part2                         <          //  processor rank  => the cell is my
        toBeDeleted = false; // it could be mine so do not delete
      }
      else if( dataStorage.Part2()[m_cells->getNodeID(dataStorage.GhostCell()[icell],inode)] < PE::GetPE().GetRank() )// if cell has a node that is from older process
      {
        toBeDeleted = true; //it belongs to older process so it does not belong to me!!
        break;
      }
      singleCell[inode]=dataStorage.GlobId()[ m_cells->getNodeID(dataStorage.GhostCell()[icell],inode) ]; // global names of nodes stored
    }
    
    if(toBeDeleted == false) // the cell is my, so it could have been sent to me
    {
      sort(singleCell.begin(), singleCell.end() ); // sort so imput is faster
      set<CFuint> cellSet;  // set so push_back works
      cellSet.insert(singleCell.begin(), singleCell.end() ); // fill set
      dataStorage.ghostCellsContainer.push_back( cellSet ); //  add to container
    }
    else if(toBeDeleted == true) // no one should have sent this one to me
      dataStorage.DelCellCol().push_back( dataStorage.GhostCell()[ icell ] ); // add to be deleted - local name
  }
  sort( dataStorage.ghostCellsContainer.begin(), dataStorage.ghostCellsContainer.end() ); // sort so binary find works

  sort( dataStorage.DelCellCol().begin(), dataStorage.DelCellCol().end() ); // sort DelCellCol() collection so binary_find works
  vector<CFuint>::iterator it;
  it = unique (dataStorage.DelCellCol().begin(), dataStorage.DelCellCol().end());     // make the save list uniqe
  dataStorage.DelCellCol().resize( it - dataStorage.DelCellCol().begin() );

  /// New Mesh elements are created...
  dataStorage.MakeNewMeshElements( PhysicalModelStack::getActive()->getDim(), (*m_states[0]).size()); // dimension and no of variables in state
  SelectNodesToDelete();
}

//////////////////////////////////////////////////////////////////////////////
/// Nodes to be removed are selected here
void StdRepart::SelectNodesToDelete()
{
  /// Select nodes that are necessary to compose saved cells
  // all nodes composing seved cells and the one recived should be saved
  //sort( dataStorage.DelCellCol().begin(), dataStorage.DelCellCol().end() ); // sort DelCellCol() collection so binary_find works
  vector<CFuint> toBeSavedNodesLocalId; // container to store saved nodes
  for(CFuint icell=0; icell < m_cells->getLocalNbGeoEnts(); ++icell) // through all m_cells - decisions are made for each cell
  if( !binary_search(dataStorage.DelCellCol().begin(), dataStorage.DelCellCol().end(), icell) ) // if icell is not to be deleted -> save it's nodes
  {
    for(CFuint inode=0; inode < m_cells->getNbNodesInGeo(icell); ++inode) // through m_nodes composing a saved cell
    {
      toBeSavedNodesLocalId.push_back(m_cells->getNodeID(icell,inode)); // save node local ID
    }
  }
  sort( toBeSavedNodesLocalId.begin(), toBeSavedNodesLocalId.end() );
  vector<CFuint>::iterator it;
  it = unique (toBeSavedNodesLocalId.begin(), toBeSavedNodesLocalId.end());     // make the save list uniqe
  toBeSavedNodesLocalId.resize( it - toBeSavedNodesLocalId.begin() );
  /// Select nodes that are necessary to be saved, so recived cells have something to be made of - global numbers
  //vector<CFuint> toBeSavedNodesGlobalId;  // - declared as a class member
  for(CFuint icell=0; icell < dataStorage.newCellVec.size(); ++icell) // through all new_cells
  {
    //vector<CFuint> tmpNodes;
    //tmpNodes = dataStorage.newCellVec[icell].GetNodes();
    toBeSavedNodesGlobalId.insert(toBeSavedNodesGlobalId.end(), dataStorage.newCellVec[icell].GetNodes().begin(), dataStorage.newCellVec[icell].GetNodes().end());
  }
  sort( toBeSavedNodesGlobalId.begin(), toBeSavedNodesGlobalId.end() );
  //vector<CFuint>::iterator it;
  it = unique (toBeSavedNodesGlobalId.begin(), toBeSavedNodesGlobalId.end());     // make the save list uniqe
  toBeSavedNodesGlobalId.resize( it - toBeSavedNodesGlobalId.begin() );
  
  for(CFuint inode=0; inode<m_nodes.size(); ++inode) // through nodes
  {
    if( binary_search(toBeSavedNodesLocalId.begin(), toBeSavedNodesLocalId.end(), inode) )  // if node found in to be saved do not at to delete collection
      continue;
    if( binary_search(toBeSavedNodesGlobalId.begin(), toBeSavedNodesGlobalId.end(), dataStorage.GlobId()[inode]) ) // if found in global
      continue;
    // if not found before add to delete node collection
    dataStorage.DelNodeCol().push_back( inode ); // sorted!
  }
}

//////////////////////////////////////////////////////////////////////////////
/// Data Storage memory is dealocated here
void StdRepart::DoClearMemory()
{
  // dataStorage owns all data
  dataStorage.DoClearMemory();
}

//////////////////////////////////////////////////////////////////////////////
template <CFuint DIM> /// Test printout
void StdRepart::DoWriteTecNoRemoved(boost::filesystem::path& path)
{
  vector<CFint> orgNamesToNewNames(m_nodes.size()); //maps orginal nodal localId to new localId; -1 -> node removed
  CFuint j=0;
  for(CFuint inode=0; inode<m_nodes.size(); ++inode) // through nodes
  {
    if( !binary_search(dataStorage.DelNodeCol().begin(), dataStorage.DelNodeCol().end(), inode) )
    {
      orgNamesToNewNames[inode] = j;
      ++j;
    }
    else
      orgNamesToNewNames[inode] = -10;
  }

  boost::filesystem::ofstream of;
  path = Framework::PathAppender::getInstance().appendParallel( path );
  of.open(path);  // file opened and ready to write!
  
  // output file header
  of<<"TITLE      =  Unstructured grid data"<<endl;
  if(DIM == 2)
    of<<"VARIABLES  =  \"x0\" \"x1\" \"CFglobId\" \"PARglobId\" \"part1\" \"part2\" ";
  else if(DIM == 3)
    of<<"VARIABLES  =  \"x0\" \"x1\" \"x2\" \"CFglobId\" \"PARglobId\" \"part1\" \"part2\" ";
  for(CFuint i=0; i<(*m_states[1]).size(); ++i)
    of<<"\"state"<<i<<"\" ";
  of<<endl;
  if(j == 0)
  {
    of<<"ZONE   T=\"ZONE0 Triag\", N="<<1<<", E="<<1<<", F=FEPOINT, ET=TRIANGLE"<<endl;
    for(CFuint i=0; i<DIM+4+(*m_states[1]).size(); ++i)
      of<<i<<" ";
    of<<endl;
    for(CFuint i=0; i<DIM+1; ++i)
      of<<"1 ";
    of<<endl;
  }
  else
  {
    of<<"ZONE   T=\"ZONE0 Triag\", N="<<j<<", E="<<m_cells->getLocalNbGeoEnts()-dataStorage.DelCellCol().size();
    if(DIM == 2)
      of<<", F=FEPOINT, ET=TRIANGLE"<<endl;
    else if(DIM == 3)
      of<<", F=FEPOINT, ET=TETRAHEDRON"<<endl;
  
    // m_nodes
    for(CFuint i=0; i<m_cells->getNbNodesInTrs(); ++i)
    if(!binary_search(dataStorage.DelNodeCol().begin(), dataStorage.DelNodeCol().end(), i))// only nodes not deleted are printed
    {
      of<<*m_nodes[i];
      of<<" "<<m_nodes[i]->getGlobalID()<<" "<<dataStorage.GlobId()[i]<<" "<<dataStorage.Part1()[i]<<" "<<dataStorage.Part2()[i]<<" "<<*m_states[i]<<endl;
    }
  
    // m_cells
    for(CFuint icell=0; icell < m_cells->getLocalNbGeoEnts(); ++icell) // through all m_cells
    if ( !binary_search(dataStorage.DelCellCol().begin(), dataStorage.DelCellCol().end(), icell) )//only cells not deleted are printed
    {
      for(CFuint inode=0; inode < m_cells->getNbNodesInGeo(icell); ++inode)
        of<<orgNamesToNewNames[m_cells->getNodeID(icell,inode)]+1<<" "; // TECPlot has different enumeration!!
      of<<endl;
    }
  }
  // close the stream to a file!
  of.close();
}

//////////////////////////////////////////////////////////////////////////////
template <CFuint DIM> /// Test printout - writes data that have been recived from others
void StdRepart::DoWriteTecAfterSendRecive(boost::filesystem::path& path)
{
  // create local map. Global ID -> Local Id for the nodes that are to be used to represent recived elements
  map<CFuint, CFuint> NodesMap;  // maps GlobalId to LocalId for orginal nodes
  map<CFuint, CFuint>::iterator it;     // to look through
  
  CFuint iPrintNode = 0;
  for(CFuint inode=0; inode < m_nodes.size(); ++inode) // through orginal nodes
  if( binary_search( toBeSavedNodesGlobalId.begin(), toBeSavedNodesGlobalId.end(), dataStorage.GlobId()[inode] ) ) // if global id of a node node will be needed
  {
    NodesMap.insert(NodesMap.end(), pair<CFuint, CFuint>(dataStorage.GlobId()[inode], iPrintNode+1) );
    ++iPrintNode;
  }
//   for(CFuint inode=0; inode<dataStorage.newNodeVec.size(); ++inode) // through recived nodes
//   {
//     NodesMap.insert(NodesMap.end(), pair<CFuint, CFuint>(dataStorage.newNodeVec[inode].GetGlobalId(), iPrintNode+1) );
//     ++iPrintNode;
//   }
  // maps to navigate through nodes are done...
  
  boost::filesystem::ofstream of;
  path = Framework::PathAppender::getInstance().appendParallel( path );
  of.open(path);  // file opened and ready to write!
  
  // output file header
  of<<"TITLE      =  Unstructured grid data"<<endl;
  if(DIM == 2)
    of<<"VARIABLES  =  \"x0\" \"x1\" \"CFglobId\" \"PARglobId\" \"part1\" \"part2\""<<" ";
  else if(DIM == 3)
    of<<"VARIABLES  =  \"x0\" \"x1\" \"x2\" \"CFglobId\" \"PARglobId\" \"part1\" \"part2\""<<" ";
  for(CFuint i=0; i<(*m_states[1]).size(); ++i)
    of<<"\"state"<<i<<"\" ";
  of<<endl;
  if(dataStorage.newCellVec.empty() )
  {
    of<<"ZONE   T=\"ZONE0 Triag\", N="<<1<<", E="<<1<<", F=FEPOINT, ET=TRIANGLE"<<endl;
    for(CFuint i=0; i<DIM+4+(*m_states[1]).size(); ++i)
      of<<i<<" ";
    of<<endl;
    for(CFuint i=0; i<DIM+1; ++i)
      of<<"1 ";
    of<<endl;
  }
  else
  {
    of<<"ZONE   T=\"ZONE0 Triag\", N="<<iPrintNode+dataStorage.newNodeVec.size();
    of<<", E="<<dataStorage.newCellVec.size();
    if(DIM == 2)
      of<<", F=FEPOINT, ET=TRIANGLE"<<endl;
    else if(DIM == 3)
      of<<", F=FEPOINT, ET=TETRAHEDRON"<<endl;
    // the saved ones
    //cout<<" aa "<<PE::GetPE().GetRank()<<" aa "<<iPrintNode<<endl;
    for(CFuint i=0; i < m_nodes.size(); ++i) // through oryginal nodes
    if( binary_search( toBeSavedNodesGlobalId.begin(), toBeSavedNodesGlobalId.end(), dataStorage.GlobId()[i] ) ) // if global id of a node node will be needed
    {
      of<<*m_nodes[i]<<" "<<m_nodes[i]->getGlobalID()<<" "<<dataStorage.GlobId()[i]<<" "<<dataStorage.Part1()[i]<<" "<<dataStorage.Part2()[i];
      of<<" "<<*m_states[i]<<endl;
      //cout<<i<<" ";
    }
    // recived nodes
    //cout<<endl<<" aa "<<PE::GetPE().GetRank()<<" aa "<<dataStorage.newNodeVec.size()<<endl;
    for(CFuint i=0; i<dataStorage.newNodeVec.size(); ++i)
    {
      for(CFuint j=0; j<dataStorage.newNodeVec[i].GetPosition().size(); ++j)
      {
        of<<dataStorage.newNodeVec[i].GetPosition()[j]<<" ";
        //cout<<dataStorage.newNodeVec[i].GetPosition()[j]<<" ";
      }
      //cout<<endl;
      of<<dataStorage.newNodeVec[i].GetCFGlobId()<<" "<<dataStorage.newNodeVec[i].GetGlobalId()<<" ";
      of<<dataStorage.newNodeVec[i].GetPart1()<<" "<<dataStorage.newNodeVec[i].GetPart2()<<" ";
      for(CFuint j=0; j<dataStorage.newStateVec[i].GetData().size(); ++j)
      {
        of<<dataStorage.newStateVec[i].GetData()[j]<<" ";
        //cout<<dataStorage.newNodeVec[i].GetPosition()[j]<<" ";
      }
      of<<endl;
      NodesMap.insert(NodesMap.end(), pair<CFuint, CFuint>(dataStorage.newNodeVec[i].GetGlobalId(), iPrintNode+1) );
      ++iPrintNode;
      //cout<<iPrintNode<<" ";
    }
  
  
    //cout<<endl;
  
    // m_cells
    for(CFuint icell=0; icell < dataStorage.newCellVec.size(); ++icell) // through all m_cells
    {
      for(CFuint inode=0; inode<dataStorage.newCellVec[icell].GetNodes().size(); ++inode)
      {
        it = NodesMap.find( dataStorage.newCellVec[icell].GetNodes()[inode] );
        if(it == NodesMap.end() ) throw Common::ParallelException (FromHere(),"Cell Adresses an unknowen node!!!");
        of<<(*it).second<<" ";
      }
      of<<endl;
    }
  }
  // close the stream to a file!
  of.close();
}

//////////////////////////////////////////////////////////////////////////////
template <CFuint DIM> /// Test printout
void StdRepart::DoWriteTec(boost::filesystem::path& path, bool interface_view)
{
  if(DIM !=2 && DIM !=3)     throw Common::NoSuchValueException(FromHere(), "No such dimension!!");

  //CFLogInfo("Writing TecPlotFile\n");

  //output file, modify file name
  boost::filesystem::ofstream of;
  path = Framework::PathAppender::getInstance().appendParallel( path );
  of.open(path);  // file opened and ready to write!

  /// I want the following structure:
  // cordinate1 cordinate2 ... global_Id owner_process_old owner_process_new some adittional data

  // output file header
  of<<"TITLE      =  Unstructured grid data"<<endl;
  if(DIM == 2)
    of<<"VARIABLES  =  \"x0\" \"x1\" \"CFglobId\" \"PARglobId\" \"part1\" \"part2\" ";
  else if(DIM == 3)
    of<<"VARIABLES  =  \"x0\" \"x1\" \"x2\" \"CFglobId\" \"PARglobId\" \"part1\" \"part2\" ";
  for(CFuint i=0; i<(*m_states[1]).size(); ++i)
    of<<"\"state"<<i<<"\" ";
  of<<endl;
  of<<"ZONE   T=\"ZONE0 Triag\", N="<<m_cells->getNbNodesInTrs()<<", E=";
  if(!interface_view)
    of<<m_cells->getLocalNbGeoEnts()-dataStorage.GhostCell().size();
  else
    of<<m_cells->getLocalNbGeoEnts();
  if(DIM == 2)
    of<<", F=FEPOINT, ET=TRIANGLE"<<endl;
  else if(DIM == 3)
    of<<", F=FEPOINT, ET=TETRAHEDRON"<<endl;

  // m_nodes
  for(CFuint i=0; i<m_cells->getNbNodesInTrs(); ++i)
  {
      of<<*m_nodes[i];
    of<<" "<<m_nodes[i]->getGlobalID()<<" "<<dataStorage.GlobId()[i]<<" "<<dataStorage.Part1()[i]<<" "<<dataStorage.Part2()[i]<<" "<<*m_states[i]<<endl;
  }

  // m_cells
  CFuint ighost = 0;
  for(CFuint icell=0; icell < m_cells->getLocalNbGeoEnts(); ++icell) // through all m_cells
  {
    if(!interface_view)
      if( dataStorage.GhostCell().size() > 0) if( icell == dataStorage.ghostCells[ighost] )
      {
        ++ighost; //ghostCells container is sorted, I do not need to search through it
        continue;
      }
    for(CFuint inode=0; inode < m_cells->getNbNodesInGeo(icell); ++inode)
      of<<m_cells->getNodeID(icell,inode)+1<<" "; // TECPlot has different enumeration!!
    of<<endl;
  }

  // close the stream to a file!
  of.close();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdRepart::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParMetisBalancer

} // namespace COOLFluiD

