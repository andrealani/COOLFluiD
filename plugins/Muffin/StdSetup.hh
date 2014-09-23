#ifndef COOLFluiD_Muffin_StdSetup_hh
#define COOLFluiD_Muffin_StdSetup_hh

#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {


/// This is a standard command to setup the Muffin method
class StdSetup : public MuffinCom {

 public:  // functions

  /// Constructor
  explicit StdSetup(const std::string& name);

  /// Destructor
  ~StdSetup();

  /// Execute processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_nodes);
    r.push_back(&s_states);
    r.push_back(&s_bstatesneighbors);
    return r;
  }

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > r;
    r.push_back(&s_nstatesproxy);
    r.push_back(&s_faceneighcell);
    r.push_back(&s_mn_volume);
    r.push_back(&s_mn_bnormal);
    r.push_back(&s_mn_bnarea);
    r.push_back(&s_mn_walldistance);
    r.push_back(&s_mn_wallnode);
    r.push_back(&s_mn_priority);
    return r;
  }

  /// Node to state ID map
  std::vector< CFuint > m_nodeIDToStateID;

  /// Set nodal states proxy
  void setNStatesProxy();

  /// Set faces neighboring cells, for the given boundaries TRS's
  void setFaceNeighCell(const std::vector< Common::SafePtr< Framework::TopologicalRegionSet > >& btrs);

  /// Set node-wise closest distance to node on a wall, and node on wall's closest node, for the given walls TRS's
  void setWallDistance(const std::vector< Common::SafePtr< Framework::TopologicalRegionSet > >& walls);

  /// Set node-wise boundary normalized normals, for the given boundaries TRS's
  void setBNodeNormals(const std::vector< Common::SafePtr< Framework::TopologicalRegionSet > >& btrs);

  /// Set node-wise volume
  void setNodeVolume();

  /// Set initial solution
  void setInitialSolution();


 private:  // sockets

  /// Socket to access nodes
  Framework::DataSocketSink< Framework::Node*,Framework::GLOBAL > s_nodes;

  /// Socket to access states
  Framework::DataSocketSink< Framework::State*,Framework::GLOBAL > s_states;

  /// Socket to access neighbor states
  Framework::DataSocketSink< std::valarray< Framework::State* > > s_bstatesneighbors;

  /// Socket to provide proxy to use the data handle of nodal states independently from storage type being RealVector or State*
  Framework::DataSocketSource< Framework::ProxyDofIterator< RealVector >* > s_nstatesproxy;

  /// Socket to provide mapping surface face to corresponding boundary element
  Framework::DataSocketSource< std::pair< CFuint,CFuint > > s_faceneighcell;

  /// Socket to provide node-wise volume
  Framework::DataSocketSource< CFreal > s_mn_volume;

  /// Socket to provide node-wise boundary normalized normals
  Framework::DataSocketSource< RealVector > s_mn_bnormal;

  /// Socket to provide node-wise dual mesh node area
  Framework::DataSocketSource< CFreal > s_mn_bnarea;

  /// Socket to provide node-wise closest distance to node on a wall
  Framework::DataSocketSource< CFreal > s_mn_walldistance;

  /// Socket to provide node-wise node on wall's closest node
  Framework::DataSocketSource< CFuint > s_mn_wallnode;

  /// Socket to provide node-wise boundary condition application priority
  Framework::DataSocketSource< CFuint > s_mn_priority;

};  // class StdSetup


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_StdSetup_hh

