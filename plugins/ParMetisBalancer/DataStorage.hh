#ifndef COOLFluiD_IO_ParMetisBalancer_DataStorage_hh
#define COOLFluiD_IO_ParMetisBalancer_DataStorage_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include <parmetis.h>
#include "Common/PE.hh"
#include "Common/CFAssert.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Common/ParallelException.hh"
#include "Framework/PartitionerData.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

#include <iostream>

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace ParMetisBalancer {
    class StdRepart;
    template <typename TYPE>
    class commMPIstruct;
    class newNode;
    class newCell;
    class newState;

  /**
   * This clas is used as a data container for repartitioner
   * It stores data that are used inside the partitioner
   *
   */
class DataStorage {
friend class COOLFluiD::ParMetisBalancer::StdRepart;

public:
  DataStorage() {}
  ~DataStorage() {}

//protected: // methods to acces data containers
  const vector<CFuint>& Part1()       {return part1;} // only to read
  vector<CFuint>& Part2()       {return part2;}
  vector<CFuint>& GlobId()      {return globNumNode;}
  vector<CFuint>& DelCellCol()  {return cellDell;}
  vector<CFuint>& DelNodeCol()  {return nodeDell;}
  vector<CFuint>& GhostCell()   {return ghostCells;}
  //vector<CFuint>& KnownNodes()  {return knownNodes;}
  vector<vector<CFuint> >& CellsToSend() {return cellSend;}
  vector<vector<CFuint> >& NodesToSend() {return nodeSend;}
  const vector<vector<CFuint> >& CommonNodes() {return commNodes;} // Only look through
  

public: // methods
  /// Common Nodes/ ghost cells
  void setupCommonNodes(const vector<vector<CFuint> >& sNodes, const vector<vector<CFuint> >& rNodes);
  bool IsCommon(CFuint inode, CFuint procrank);
  
  /// Resize
  void ResizeData(const CFuint noNodes);
  
  ///Set/Get Node/Cell Send
  void   UniqeSendNodes(const CFuint destProc); // Make a list of send nodes uniqe
  
  /// containers for MPI send/recive operations are being setup here
  //void InitMPIcommStruct();
  
  /// Fills MPI send struct with data
//   void FillIntSendStruct(const vector<vector<CFuint> >& nodeData, const vector<vector<CFuint> >& cellData, const vector<vector<CFuint> >& sData); // int
//   void FillDouSendStruct(const vector<vector<CFdouble> >& nodeData, const vector<vector<CFdouble> >& cellData, const vector<vector<CFdouble> >& sData); // double
  
  /// MPI comunication is caried out
  void ExecuteMPIcomm();
  
  /// makes lists of nodes against wchich send ones will be chacked
  void MakeListsOfKnownNodes();
  
  /// translates MPI recive vectors into NewNode/NewCell colections
  void MakeNewMeshElements(const CFuint dim, const CFuint statesize);
  
  ///Dealocates memory
  void DoClearMemory();

protected: // data members accesible from outside

  // specify to where the node belonges - old part
  vector<CFuint> part1; // Part1
  // specify to where it belonges now - new part - filled by ParMETIS!
  vector<CFuint> part2;  // Part2
  // global numbers
  vector<CFuint> globNumNode; // global node mapping
  // list of cells to be send
  vector<vector<CFuint> > cellSend; // Cells to send
  // list of nodes to be send
  vector<vector<CFuint> > nodeSend; // nodes to send
  // cells for which other processes decide
  vector<CFuint> ghostCells;
  // list of cells/nodes to be removed
  vector<CFuint> cellDell;
  vector<CFuint> nodeDell;

protected: // Mesh elements to be added
  void AddNewNode(const CFuint cfGlobId, const CFuint globId, const CFuint part1, const CFuint part2, const vector<double>& position, const vector<double>& stateVec);
  void AddNewCell(const vector<CFuint>& newCell, const CFuint dim);
  void AddNewState(const vector<double>& stateVec);
  vector<newNode> newNodeVec;
  vector<newCell> newCellVec;
  vector<newState>newStateVec;
  
protected: // MPI communication structures. [0] Node, [1] cell, [2] state

  vector<commMPIstruct<CFuint> >    m_intData;  // int type data for node [0] and cell [1]
  vector<commMPIstruct<CFdouble> >  m_douData;  // double type data for node [0] and cell [1]
  
protected:
  // Sorted list of nodes known to other processes (there might be more of those at triple points)
  vector< vector<CFuint> > commNodes;
  // known nodes
  vector<CFuint> knownNodes;
  // known cells
  vector<set<CFuint> > ghostCellsContainer;
  
protected: // CSR format data for ParMETIS
  vector<Framework::PartitionerData::IndexT>     xadj;         //CSR format: navigation helper for adjncy
  vector<Framework::PartitionerData::IndexT>     adjncy;       //CSR format: conectivity of nodes (global indexes)
  vector<Framework::PartitionerData::IndexT>     vtxdist;      //CSR amount of nodes processes own
}; // Data Storage end


class newNode{
  public:
    newNode() {throw Common::ParallelException (FromHere(),"Not to be used!!");}; // not to be used
    newNode(const CFuint cfGlobId, const CFuint globId, const CFuint part1, const CFuint part2, const vector<double>& position) :
      CFGlobId(cfGlobId), GlobalId(globId), Part1(part1), Part2(part2), m_xyz(position) {}
    ~newNode() {};
  public: // acess
    CFuint GetCFGlobId() {return CFGlobId;}
    CFuint GetGlobalId() {return GlobalId;}
    CFuint GetPart1() {return Part1;}
    CFuint GetPart2() {return Part2;}
    const vector<double>& GetPosition() {return m_xyz;}
    
  protected: // members
    CFuint CFGlobId; // coolfluid global id
    CFuint GlobalId; // global id which I am using
    CFuint Part1;
    CFuint Part2;
    vector<double> m_xyz; // position of a node
};
// this class stores data which compose a new cell
class newCell{
  public:
    newCell() {throw Common::ParallelException (FromHere(),"Not to be used!!");};
    newCell(const vector<CFuint>& nodeList) : nodes(nodeList) {}
    ~newCell() {};
  
  public: // acess
    const vector<CFuint>& GetNodes() {return nodes;}
    
  protected: // members
    vector<CFuint> nodes; // list of nodes composing this cell
};
// this class stores data which compose a new state
class newState{
  public:
    newState() {throw Common::ParallelException (FromHere(),"Not to be used!!");}; // not to be used
    newState(const vector<CFdouble>& imput) : stateDataVec(imput) {}
    ~newState() {};
  public: // acess
    const vector<double>& GetData() {return stateDataVec;}
    
  protected: // members
    vector<double> stateDataVec;
};

/** This class is used to carry MPI data exchange,it,
  * 1) holds data of appropriate type in m_MPISend
  * 2) comunicates the volume of communication to other processes
  * 3) does MPI data exchange
  * Data containers ( MPISend and MPIreci ) are accesible outside through reference (think of something else?)
  */
template <typename TYPE>
class commMPIstruct{

public:
  commMPIstruct() { Init(); }; // MPIStruct is initialized at construction
  ~commMPIstruct() {};

public: // methods to acess conteiners:
  vector<vector<TYPE> >& GetReciData() {return m_MPIReci;}
  void AddSendElem(const CFuint proc, const TYPE value) {m_MPISend[proc].push_back( value );}

protected:
  // initialize
  void Init();

public: // methods for MPI communication

  // sizes for recive bufers are comunicated to appropriate proccesses
  void SetSizeRecivBuf();
  // MPI comunication session
  void MPIcommunicate();
  // Raprorts what has been send/recived TODO: will be removed -- does a lot of output
  void MPIraport();
  
private: // members
  vector<vector<TYPE> > m_MPISend; // bufer for send info
  vector<vector<TYPE> > m_MPIReci; // bufer for recive info
};  // commMPIstruct end

template <typename TYPE>
inline void commMPIstruct<TYPE>::Init()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  m_MPISend.resize(PE::GetPE().GetProcessorCount(nsp));
  m_MPIReci.resize(PE::GetPE().GetProcessorCount(nsp));
}

template <typename TYPE>
inline void commMPIstruct<TYPE>::SetSizeRecivBuf()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  vector<CFuint> out(PE::GetPE().GetProcessorCount(nsp), 0);
  vector<CFuint> in (PE::GetPE().GetProcessorCount(nsp), 0);
  
  for(CFuint i=0; i<PE::GetPE().GetProcessorCount(nsp); ++i)
    out[i] = m_MPISend[i].size();
  
  MPI_Alltoall(&(out[0]) , 1, MPI_INT, &(in[0]) , 1, MPI_INT, PE::GetPE().GetCommunicator(nsp));
  
  for(CFuint i=0; i<PE::GetPE().GetProcessorCount(nsp); ++i)
    m_MPIReci[i].resize( in[i] );
}

template <typename TYPE>
inline void commMPIstruct<TYPE>::MPIcommunicate()
{
  throw Common::ParallelException (FromHere(),"Not implemented for any TYPE");
}

template <>
inline void commMPIstruct<CFuint>::MPIcommunicate()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  vector<MPI_Request>	tabreqs(  2*PE::GetPE().GetProcessorCount(nsp) );
  vector<MPI_Status>	tabstats( 2*PE::GetPE().GetProcessorCount(nsp) );
  int tag = 1;
  int licznik=0;
  //----------------- recive
  for(CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p) 
    if(p!=PE::GetPE().GetRank(nsp) && m_MPIReci[p].size() > 0)
    { 
      MPI_Irecv( &m_MPIReci[p][0], m_MPIReci[p].size(), 
	Common::MPIStructDef::getMPIType(&m_MPIReci[p][0]), p, tag, 
		 PE::GetPE().GetCommunicator(nsp), &tabreqs[licznik]);
      ++licznik;
    }
  // ----------------  send
  for(CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p) 
    if(p!=PE::GetPE().GetRank(nsp) && m_MPISend[p].size() > 0)
      {
	//cout<<mProcId<<": pSize["<<p<<"]="<<pSize[p]<<" pSizeIn[p]"<<pSizeIn[p]<<endl;
      MPI_Isend(&m_MPISend[p][0], m_MPISend[p].size(), 
	Common::MPIStructDef::getMPIType(&m_MPISend[p][0]), p, tag, PE::GetPE().GetCommunicator(nsp), &tabreqs[licznik]); //wysylam id
      ++licznik;
    }
  MPI_Waitall( licznik, &tabreqs[0], &tabstats[0]);
}

template <>
inline void commMPIstruct<CFdouble>::MPIcommunicate()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  vector<MPI_Request>	tabreqs(  2*PE::GetPE().GetProcessorCount(nsp) );
  vector<MPI_Status>	tabstats( 2*PE::GetPE().GetProcessorCount(nsp) );
  int tag = 1;
  int licznik=0;
  //----------------- recive
  for(CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p) 
    if(p!=PE::GetPE().GetRank(nsp) && m_MPIReci[p].size() > 0)
      { 
      MPI_Irecv( &m_MPIReci[p][0], m_MPIReci[p].size(), MPI_DOUBLE, p, tag, PE::GetPE().GetCommunicator(nsp), &tabreqs[licznik]);
      ++licznik;
    }
  // ----------------  send
  for(CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p) 
    if(p!=PE::GetPE().GetRank(nsp) && m_MPISend[p].size() > 0)
    {
      //cout<<mProcId<<": pSize["<<p<<"]="<<pSize[p]<<" pSizeIn[p]"<<pSizeIn[p]<<endl;
      MPI_Isend(&m_MPISend[p][0], m_MPISend[p].size(), MPI_DOUBLE, p, tag, PE::GetPE().GetCommunicator(nsp), &tabreqs[licznik]); //wysylam id
      ++licznik;
    }
  MPI_Waitall( licznik, &tabreqs[0], &tabstats[0]);
}

template <typename TYPE>
inline void commMPIstruct<TYPE>::MPIraport()
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  for(CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p)
  {
    cout<<"To "<<p<<" I have sent "<<m_MPISend[p].size()<<" elements: "<<endl;
    for(CFuint i=0; i<m_MPISend[p].size(); ++i)
      cout<<m_MPISend[p][i]<<" ";
    cout<<endl;
  }
  
  for(CFuint p=0; p<PE::GetPE().GetProcessorCount(nsp); ++p)
    {
    cout<<"From "<<p<<" I have recived "<<m_MPIReci[p].size()<<" elements: "<<endl;
    for(CFuint i=0; i<m_MPIReci[p].size(); ++i)
      cout<<m_MPIReci[p][i]<<" ";
    cout<<endl;
    }
}
  }
}
#endif

//////////////////////////////////////////////////////////////////////////////
