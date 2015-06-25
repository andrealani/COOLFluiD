#include "DataStorage.hh"

#include <algorithm>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

    namespace ParMetisBalancer {
    
//////////////////////////////////////////////////////////////////////////////
// Setup
void DataStorage::setupCommonNodes(const vector<vector<CFuint> >& sNodes, const vector<vector<CFuint> >& rNodes)
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  cf_assert(commNodes.size() == PE::GetPE().GetProcessorCount(nsp));
  // Setup CommonNodes container
  for(CFuint i=0; i<commNodes.size(); ++i)
  {
    commNodes[i].insert(commNodes[i].end(), sNodes[i].begin(), sNodes[i].end());
    commNodes[i].insert(commNodes[i].end(), rNodes[i].begin(), rNodes[i].end());
    sort( commNodes[i].begin(), commNodes[i].end() ); // sort the entries for binary_search
  }
  // Setup part1 coloring
  cf_assert(part1.size() > 0 );
  for(CFuint i=0; i<rNodes.size(); ++i)
  {
    for(CFuint j=0; j<rNodes[i].size(); ++j)
    {
      part1[ rNodes[i][j] ] = i;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
bool DataStorage::IsCommon(CFuint inode, CFuint procrank)
{
  return binary_search( commNodes[procrank].begin(), commNodes[procrank].end(), inode);
}

//////////////////////////////////////////////////////////////////////////////
void DataStorage::ResizeData(const CFuint noNodes)
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  commNodes.resize( PE::GetPE().GetProcessorCount(nsp) ); // CommonNodes
  part1.resize( noNodes, PE::GetPE().GetRank(nsp) ); // Most of nodes belong to process
  part2.resize( noNodes, 0); // this will be set after ParMETIS is run
  globNumNode.resize(noNodes);
  cellSend.resize( PE::GetPE().GetProcessorCount(nsp) );
  nodeSend.resize( PE::GetPE().GetProcessorCount(nsp) );
  m_intData.resize(3); // for nodes [0], cells [1], states [2] - int
  m_douData.resize(3); // for nodes [0], cells [1], states [2] - double
}

//////////////////////////////////////////////////////////////////////////////
void DataStorage::UniqeSendNodes(const CFuint destProc)
{
  sort( nodeSend[destProc].begin(), nodeSend[destProc].end() ); // sort
  vector<CFuint>::iterator it;
  it = unique (nodeSend[destProc].begin(), nodeSend[destProc].end());     // uniqe
  nodeSend[destProc].resize( it - nodeSend[destProc].begin() );
}

//////////////////////////////////////////////////////////////////////////////
void DataStorage::ExecuteMPIcomm()
{
  for(uint i=0; i<m_intData.size(); ++i)
  {
    m_intData[i].SetSizeRecivBuf();
    m_intData[i].MPIcommunicate();
  }
  for(uint i=0; i<m_douData.size(); ++i)
  {
    m_douData[i].SetSizeRecivBuf();
    m_douData[i].MPIcommunicate();
  }
}

//////////////////////////////////////////////////////////////////////////////
void DataStorage::MakeListsOfKnownNodes()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  
  for( CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p ) {
    if(p != PE::GetPE().GetRank(nsp) )
      knownNodes.insert(knownNodes.end(), CommonNodes()[p].begin(), CommonNodes()[p].end());
  }
  
  // rewrite KnownNodes() into global names
  for(CFuint i=0; i<knownNodes.size(); ++i)
    knownNodes[i] = GlobId()[ knownNodes[i] ]; // put glob num instead of local num
  sort( knownNodes.begin(), knownNodes.end() ); // sort for future find
}

///////////////////////////////////////////////////////////////////////////////
void DataStorage::AddNewCell(const vector<CFuint>& nCell, const CFuint dim)
{
  cf_assert(nCell.size() == dim);
  
  set<CFuint> testCell;
  testCell.insert( nCell.begin(), nCell.end() );
  if( !binary_search( ghostCellsContainer.begin(), ghostCellsContainer.end(), testCell ) )  // if not found, add
  {
    newCell TmpCell( nCell );
    newCellVec.push_back( TmpCell ); /// NewCell is described by GlobId!!
  }
}

//////////////////////////////////////////////////////////////////////////////
void DataStorage::AddNewNode(const CFuint cfGlobId, const CFuint globId, const CFuint part1, const CFuint part2, const vector<double>& position, const vector<double>& stateVec)
{
  cf_assert(1);
  
  if( !binary_search( knownNodes.begin(), knownNodes.end(), globId ) )
  {
    newNode  newTmpNode(cfGlobId, globId, part1, part2, position);
    newNodeVec.push_back( newTmpNode );
    AddNewState(stateVec);
    
    knownNodes.push_back( globId ); // add as known
    sort( knownNodes.begin(), knownNodes.end() ); //sort
  }
}

//////////////////////////////////////////////////////////////////////////////
void DataStorage::AddNewState(const vector<double>& stateVec)
{
  newState newTmpState( stateVec );
  newStateVec.push_back(newTmpState);
}

//////////////////////////////////////////////////////////////////////////////
void DataStorage::MakeNewMeshElements(const CFuint dim, const CFuint statesize)
{
  
  CFuint nodesPerCell;  // how many nodes makes one cell - SIMPLEX!!
  if(dim == 2)  nodesPerCell = 3;
  else if(dim == 3) nodesPerCell = 4;
  else if(dim != 2 || dim !=3) throw Common::ParallelException (FromHere(),"Dimension Problem while creating new cells!!");

  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  // add new cells
  for(CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p) if( p != PE::GetPE().GetRank(nsp) )// through processes that communicated
  {
    // for now cell has no double type data association
    cf_assert( m_douData[1].GetReciData()[p].empty() );
    
    while( m_intData[1].GetReciData()[p].size() != 0 ) // if the cell recive container is not empty
    {
      vector<CFuint> tmpCell(nodesPerCell);
      for(CFuint i=0; i<nodesPerCell; ++i)
      {
        cf_assert( !m_intData[1].GetReciData()[p].empty() );  //should not be empty!!
        tmpCell[nodesPerCell-i-1] = m_intData[1].GetReciData()[p].back();  //! copy into tmpCell (reverse !! -> nodesPerCell-i-1)
        m_intData[1].GetReciData()[p].pop_back(); // delete last entry
      }
      AddNewCell( tmpCell, dim );
    }
  }
  
  // add new nodes
  for(CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p) if( p != PE::GetPE().GetRank(nsp) )// through processes that communicated
  {
    // since pop_back is used I am coping in reverse direction than it is stored!!
    CFuint Part2; // part2CFuint CFGlobId; // coolfluid global id
    CFuint Part1;
    CFuint GlobalId; // global id which I am using
    CFuint CFGlobId;
    vector<CFdouble> m_xyz(dim); // position of a node
    vector<CFdouble> stateVec(statesize); // state vars
    while( m_intData[0].GetReciData()[p].size() != 0 && m_douData[0].GetReciData()[p].size() != 0 ) // if the node recive container is not empty
    {
      Part2    = m_intData[0].GetReciData()[p].back(); m_intData[0].GetReciData()[p].pop_back();
      Part1    = m_intData[0].GetReciData()[p].back(); m_intData[0].GetReciData()[p].pop_back();
      GlobalId = m_intData[0].GetReciData()[p].back(); m_intData[0].GetReciData()[p].pop_back();
      cf_assert( !m_intData[0].GetReciData()[p].empty() );
      CFGlobId = m_intData[0].GetReciData()[p].back(); m_intData[0].GetReciData()[p].pop_back();
      for(CFuint i=0; i<dim; ++i) // xyz has a size of dim
      {
        cf_assert( !m_douData[0].GetReciData()[p].empty() );
        m_xyz[dim-i-1] = m_douData[0].GetReciData()[p].back();  // copy into position (reverse !! -> dim-i-1)
        m_douData[0].GetReciData()[p].pop_back(); // delete last entry
      }
      for(CFuint i=0; i<statesize; ++i) // state variables
      {
        cf_assert( !m_douData[2].GetReciData()[p].empty() );
        stateVec[statesize-i-1] = m_douData[2].GetReciData()[p].back();// copy into position (reverse !! -> dim-i-1)
        m_douData[2].GetReciData()[p].pop_back(); // delete last entry
      }
      AddNewNode(CFGlobId, GlobalId, Part1, Part2, m_xyz, stateVec);
    }
    cf_assert( m_intData[0].GetReciData()[p].empty() );
    cf_assert( m_douData[0].GetReciData()[p].empty() );
  }
}

//////////////////////////////////////////////////////////////////////////////
void DataStorage::DoClearMemory()
{
  part1.clear();
  part2.clear();
  globNumNode.clear();
  cellSend.clear();
  nodeSend.clear();
  ghostCells.clear();
  cellDell.clear();
  nodeDell.clear();
  newNodeVec.clear();
  newCellVec.clear();
  newStateVec.clear();
  m_douData.clear();
  m_intData.clear();
  commNodes.clear();
  knownNodes.clear();
  ghostCellsContainer.clear();
  xadj.clear();
  adjncy.clear();
  vtxdist.clear();
}
    }
}
