#ifndef COOLFluiD_RadiativeTransfer_ParallelVector_hh
#define COOLFluiD_RadiativeTransfer_ParallelVector_hh

#include "Framework/SocketBundleSetter.hh"
#include "Framework/MeshData.hh"
#include <vector>
#include "LagrangianSolver/LagrangianSolverModule.hh"
#include "Common/MPI/MPIStructDef.hh"

namespace COOLFluiD{

namespace LagrangianSolver {


template<typename T>
class ParallelVector : public Framework::SocketBundleSetter{
public:
  
  ParallelVector() : Framework::SocketBundleSetter(), m_isSetup(false) {}
  
  ~ParallelVector() {}
  
  void sincronizeAdd();

  void sincronizeAssign();

  void getSharedEntries();
  
  void resize(CFuint nbPoints, T value = static_cast<T>( CFNULL ) ){ m_data.resize(nbPoints, value); } 
  T& operator[](CFuint idx){return (m_data[idx]); }
  
  const T& operator[](CFuint idx) const { return (m_data[idx]); }
  
  ParallelVector<T>& operator = (T val) 
  {
    for(CFuint i=0; i<m_data.size();++i ){
      m_data[i] = val;
    }
    return *this;
  } 
  
  CFuint size(){ return m_data.size(); }
  
private:

  std::vector<T> m_data;

  std::vector<CFuint> m_sendVectorOverlapLocalIDs;
  std::vector<CFuint> m_recvVectorOverlapLocalIDs;

  std::vector<int> m_sendVectorOverlapRankDisps;
  std::vector<int> m_recvVectorOverlapRankDisps;
  
  std::vector<int> m_sendVectorOverlapRankCounts;
  std::vector<int> m_recvVectorOverlapRankCounts;
  
  CFuint m_myProcessRank;
  MPI_Comm m_comm;
  CFuint m_nbProcesses;

  bool m_isSetup;
};


template<typename T>
void ParallelVector<T>::sincronizeAdd()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  if( Common::PE::GetPE().GetProcessorCount(nsp) == 1 ) {
    return;
  }
  
  if ( !m_isSetup ){
    getSharedEntries(); 
  }

  std::vector<CFuint> sendBuffer( m_sendVectorOverlapLocalIDs.size() );
  std::vector<CFuint> recvBuffer( m_recvVectorOverlapLocalIDs.size() );

  for( CFuint i =0; i < m_sendVectorOverlapLocalIDs.size(); ++i){
    sendBuffer[i] = m_data[ m_sendVectorOverlapLocalIDs[i] ];
  }
  
  MPI_Alltoallv(&sendBuffer[0], &m_sendVectorOverlapRankCounts[0], &m_sendVectorOverlapRankDisps[0], 
		Common::MPIStructDef::getMPIType(&sendBuffer[0]),
		&recvBuffer[0], &m_recvVectorOverlapRankCounts[0], &m_recvVectorOverlapRankDisps[0], 
		Common::MPIStructDef::getMPIType(&recvBuffer[0]), m_comm);
  
  for(CFuint i=0; i< m_recvVectorOverlapLocalIDs.size(); ++i){
    m_data[ m_recvVectorOverlapLocalIDs[i] ] += recvBuffer[i];
  }
  
}

template<typename T>
void ParallelVector<T>::sincronizeAssign()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  if(  Common::PE::GetPE().GetProcessorCount(nsp) == 1 ){
    return;
  }

 if ( !m_isSetup ){
     getSharedEntries(); 
  }
  std::vector<CFuint> sendBuffer( m_sendVectorOverlapLocalIDs.size() );
  std::vector<CFuint> recvBuffer( m_recvVectorOverlapLocalIDs.size() );

  for( CFuint i =0; i < m_sendVectorOverlapLocalIDs.size(); ++i){
    sendBuffer[i] = m_data[ m_sendVectorOverlapLocalIDs[i] ];
  }

  MPI_Alltoallv(&sendBuffer[0], &m_sendVectorOverlapRankCounts[0], &m_sendVectorOverlapRankDisps[0], 
		Common::MPIStructDef::getMPIType(&sendBuffer[0]),
		&recvBuffer[0], &m_recvVectorOverlapRankCounts[0], &m_recvVectorOverlapRankDisps[0],
		Common::MPIStructDef::getMPIType(&recvBuffer[0]), m_comm);
  
  for(CFuint i=0; i< m_recvVectorOverlapLocalIDs.size(); ++i){
    m_data[ m_recvVectorOverlapLocalIDs[i] ] = recvBuffer[i];
  }
}

template<typename T>
void ParallelVector<T>::getSharedEntries()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  if(  Common::PE::GetPE().GetProcessorCount(nsp) == 1) {
    return;
  }
  
  m_isSetup = true;
  // MPI parameters
  m_myProcessRank = Common::PE::GetPE().GetRank(nsp);
  m_nbProcesses   = Common::PE::GetPE().GetProcessorCount(nsp);
  m_comm          = Common::PE::GetPE().GetCommunicator(nsp);
  
  m_sendVectorOverlapRankDisps.resize(m_nbProcesses,0);
  m_recvVectorOverlapRankDisps.resize(m_nbProcesses,0);

  m_sendVectorOverlapRankCounts.resize(m_nbProcesses,0);
  m_recvVectorOverlapRankCounts.resize(m_nbProcesses,0);

  std::map<CFuint, CFuint> globalLocalGhostMap;
  std::map<CFuint, CFuint> globalLocalOwnerMap;

  Framework::CellTrsGeoBuilder::GeoData &cellsData = m_cellBuilder.getDataGE();
  const CFuint nbCells = cellsData.trs->getLocalNbGeoEnts();

  for(CFuint i=0; i<nbCells; ++i){
    cellsData.idx = i;
    Framework::GeometricEntity *const cell = m_cellBuilder.buildGE();

    if (cell->getState(0)->isParUpdatable() ){
      globalLocalGhostMap.insert( std::pair<CFuint,CFuint>( cell->getState(0)->getGlobalID(),
                                                       cell->getState(0)->getLocalID()    )  );
    }
    else{
      globalLocalOwnerMap.insert( std::pair<CFuint,CFuint>( cell->getState(0)->getGlobalID(),
                                                       cell->getState(0)->getLocalID()    )  );

    }
    m_cellBuilder.releaseGE();
  }

  CFuint maxNbCells = 0;
  CFuint localNbCells = globalLocalGhostMap.size() ;
  MPI_Allreduce(&localNbCells, &maxNbCells, 1, 
		Common::MPIStructDef::getMPIType(&localNbCells), MPI_MAX, m_comm);
  
  std::vector<CFuint> bCastGlobalCellIDs; bCastGlobalCellIDs.resize(maxNbCells+1);
  std::vector<CFuint> replyGlobalCellIDs; replyGlobalCellIDs.reserve(maxNbCells+1);
  std::vector<CFuint> replyProcessorRank; replyProcessorRank.reserve(maxNbCells+1);
  std::vector<CFuint> gatherGlobalCellIDs; gatherGlobalCellIDs.reserve(maxNbCells+1);
  
  m_recvVectorOverlapLocalIDs.reserve(maxNbCells);

  std::map<CFuint,CFuint>::iterator it;

  for (CFuint processRank= 0; processRank < m_nbProcesses; ++ processRank) {

    if (processRank == m_myProcessRank) {

      // first entry is the size of the overlap to be broadcast
      bCastGlobalCellIDs[0]=globalLocalGhostMap.size();

      CFuint idx=1;
      for (it = globalLocalGhostMap.begin() ; it != globalLocalGhostMap.end(); ++it , ++idx){
        bCastGlobalCellIDs[idx]=(it->first);
      }
    }
    
    MPI_Bcast(&bCastGlobalCellIDs[0], bCastGlobalCellIDs.size(),
	      Common::MPIStructDef::getMPIType(&bCastGlobalCellIDs[0]), processRank, m_comm);
    
    replyGlobalCellIDs.clear();

    if ( processRank != m_myProcessRank) {

      const CFuint nbC = bCastGlobalCellIDs[0]+1;

      m_recvVectorOverlapRankDisps[processRank] = m_recvVectorOverlapLocalIDs.end() - m_recvVectorOverlapLocalIDs.begin();

      for (CFuint c = 1; c < nbC; ++c) {
        it = globalLocalOwnerMap.find(bCastGlobalCellIDs[c]) ;

        if ( it != globalLocalOwnerMap.end() ) {
          replyGlobalCellIDs.push_back( bCastGlobalCellIDs[c] );
          m_recvVectorOverlapLocalIDs.push_back( it->second ) ;
          m_recvVectorOverlapRankCounts[processRank]++;
        }
      }
    }


    std::vector<int> gatherCounts(m_nbProcesses);
    std::vector<int> gatherDisps(m_nbProcesses);

    int nbReply = int(replyGlobalCellIDs.size());
    
    MPI_Gather(&nbReply, 1, Common::MPIStructDef::getMPIType(&nbReply), &gatherCounts[0], 
	       1, Common::MPIStructDef::getMPIType(&nbReply), processRank, m_comm);
    
    CFint totalSendCount=0;
    for(CFuint i=0;i<gatherCounts.size(); ++i){
      gatherDisps[i] = totalSendCount;
      totalSendCount += gatherCounts[i];
    }
    
    gatherGlobalCellIDs.resize( totalSendCount );
    MPI_Gatherv( &replyGlobalCellIDs[0], int(replyGlobalCellIDs.size()), 
                 Common::MPIStructDef::getMPIType(&replyGlobalCellIDs[0]), 
                 &gatherGlobalCellIDs[0], &gatherCounts[0], &gatherDisps[0], 
                 Common::MPIStructDef::getMPIType(&gatherGlobalCellIDs[0]), 
		             int(processRank), m_comm );
    
    if(m_myProcessRank == processRank){
      
      m_sendVectorOverlapLocalIDs.resize( gatherGlobalCellIDs.size() );

      for(CFuint i=0 ; i< gatherGlobalCellIDs.size(); ++i ){
        it = globalLocalGhostMap.find(gatherGlobalCellIDs[i]) ;

        cf_assert(it != globalLocalGhostMap.end());

        m_sendVectorOverlapLocalIDs[i]=it->second;
      }
      m_sendVectorOverlapRankDisps  = gatherDisps;
      m_sendVectorOverlapRankCounts = gatherCounts;

    }
  }

}


}

}


#endif
