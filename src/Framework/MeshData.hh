// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshData_hh
#define COOLFluiD_Framework_MeshData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SelfRegistPtr.hh"
#include "Common/SafePtr.hh"
#include "Common/CFMap.hh"
#include "Config/ConfigObject.hh"

#include "Framework/MeshStatistics.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Framework/ElementTypeData.hh"
#include "Framework/Storage.hh"
#include "Framework/NamespaceGroup.hh"
#include "Framework/NamespaceStack.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common { 
    template <typename TYPE> class GeneralStorage; 
    class FactoryRegistry;
  }

  namespace Framework {

    class Node;
    class State;
    class MapGeoToTrsAndIdx;
    class MeshDataStack;
    class DomainModel;
    
//////////////////////////////////////////////////////////////////////////////

/// This class offers a facade for the Data realted to the Mesh,
/// to TopologicalRegionSet's and to Mesh ConnectivityTable's.
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Dries Kimpe
class Framework_API MeshData :
    public NamespaceGroup,
    public Config::ConfigObject,
    public Common::NonCopyable<MeshData>
{

friend class MeshDataStack;

public: // typedefs

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Typedef for ConnectivityTable's in MeshData
  typedef Common::ConnectivityTable<CFuint> ConnTable;

public: // methods

  /// Configures this MeshData.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Instructs MeshData to prepare for loading the mesh.
  /// It places the empty storage of states, nodes and elements in the DataStorage.
  void allocateSockets();

  /// Instructs MeshData to clean up after the mesh has been used.
  /// It deallocates the storage of states, nodes and elements in the DataStorage.
  void deallocateSockets();

  ///Deallocates the Nodes and the States
  void deallocateNodesStates();

  /// Instructs MeshData to allocate some default connectivity in the storage
  void allocateConnectivity();

  /// Instructs MeshData to clean up the connectivity storage
  void deallocateConnectivity();

  /// Deallocate function to allow deallocation of memory whenever
  /// needed during the simulation. Deallocation would be
  /// impossible otherwise, since MeshData is a singleton
  /// and the destructor cannot explicitly be called.
  void deallocate();

  /// Reallocate function to allow reallocation of memory
  /// (following a deallocation) whenever needed during
  /// the simulation. Reallocation would be impossible
  /// otherwise, since MeshData is a singleton
  /// and the constructor cannot explicitly be called.
  void reallocate();

  /// Set the factory registry
  void setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr);
  
  /// Get the factory registry
  Common::SafePtr<Common::FactoryRegistry> getFactoryRegistry();
  
  /// Access the mesh statistics
  MeshStatistics& Statistics() { return m_statistics; }
  
  /// Set the number of overlap layers
  void setNbOverlapLayers(CFuint nbOverlapLayers) {m_nbOverLayers = nbOverlapLayers;}
  
  /// Get the number of overlap layers
  CFuint getNbOverlapLayers() const {return m_nbOverLayers;}
  
  /// Get the provided sockets by the MeshData
  std::vector<Common::SafePtr<BaseDataSocketSource> > getProvidedSockets();

  /// @returns the DomainModel of this data
  Common::SafePtr< DomainModel > getDomainModel();

  /// @return the data storage
  Common::SafePtr< DataStorage > getDataStorage();

  /// @return the elementType's data
  Common::SafePtr< std::vector<ElementTypeData> > getElementTypeData();

  /// creates the element type data for the TRS
  void createElementTypeData(const std::string trsName);

  /// @return the elementType's data for this TRS
  Common::SafePtr< std::vector<ElementTypeData> > getElementTypeData(const std::string trsName);

  /// @return the ConnectivityTable corresponding to the given name
  /// @param name the name to identify the connectivity
  Common::SafePtr<ConnTable> getConnectivity(const std::string& name);

  /// Store the given ConnectivityTable
  /// @param name the name to identify the connectivity
  /// @param conn the pointer to the connvectivity
  void storeConnectivity(const std::string& name, ConnTable *const conn);

  /// Removes and deletes the given ConnectivityTable from Connectivity Storage
  /// @param name the name to identify the connectivity
  void deleteConnectivity(const std::string& name);

  /// Remove only the given ConnectivityTable from Connectivity Storage
  /// @param name the name to identify the connectivity
  void removeConnectivity(const std::string& name);

  /// @return the MapGeoToTrsAndIdx corresponding to the given name
  /// @param name the name to identify the connectivity
  Common::SafePtr<MapGeoToTrsAndIdx> getMapGeoToTrs(const std::string& name);

  /// Store the given MapGeoToTrsAndIdx
  /// @param name the name to identify the MapGeoToTrsAndIdx
  /// @param conn the pointer to the  MapGeoToTrsAndIdx
  void storeMapGeoToTrs(const std::string& name, MapGeoToTrsAndIdx *const mapG);

  /// Remove only the given MapGeoToTrsAndIdx from MapGeoToTrsAndIdx Storage
  /// @param name the name to identify the connectivity
  void removeMapGeoToTrs(const std::string& name);

  /// Returns a vector with pointers to all the TRSs defined in this MeshData
  /// @post Do not assume the order of the TRSs inside this vector will remain constant.
  /// @return the list of TopologicalRegionSet's
  std::vector< Common::SafePtr<TopologicalRegionSet> > getTrsList();

  /// Returns a vector with pointers to all the TRSs defined in this MeshData satisfying all
  /// the tags given in the tagList
  /// @param tagList the list of the tags that the TRS need to have
  /// @post Do not assume the order of the TRSs inside this vector will remain constant.
  /// @return the list of TopologicalRegionSet's which satisfy all the tags
  std::vector< Common::SafePtr<TopologicalRegionSet> >
    getFilteredTrsList(const std::vector<std::string>& tagList);

  /// Returns a vector with pointers to all the TRSs defined in this MeshData satisfying
  /// the tags given
  /// @param tag the tag that the TRS need to have
  /// @post Do not assume the order of the TRSs inside this vector will remain constant.
  /// @return the list of TopologicalRegionSet's which satisfy the tag
  std::vector< Common::SafePtr<TopologicalRegionSet> >
    getFilteredTrsList(const std::string& tag);

  /// Adds the TopologicalRegionSet to the storage
  /// @param trs the TRS pointer to store
  void addTrs(TopologicalRegionSet* trs);

  /// @return the TopologicalRegionSet by its name
  Common::SafePtr<TopologicalRegionSet> getTrs(const std::string& name);

  /// Get the local number of nodes coming from the "nodes" socket
  CFuint getNbNodes()
  {
    const CFuint nbNodes = socket_nodes.getDataHandle().size();
    return nbNodes;
  }

  /// Get the local number of states coming from the "states" socket
  CFuint getNbStates()
  {
    const CFuint nbStates = socket_states.getDataHandle().size();
    return nbStates;
  }

  /// Get the states datasocket
  Framework::DataSocketSink< Framework::State*, Framework::GLOBAL> getStateDataSocketSink();

  /// Get the nodes datasocket
  Framework::DataSocketSink< Framework::Node*, Framework::GLOBAL> getNodeDataSocketSink();

  /// Set the total number of nodes for all processors
  void setTotalNodeCount (CFuint Nc)
  {
      m_totalNodes = Nc;
  }

  /// Set the total number of state for all processors
  void setTotalStateCount (CFuint Nc)
  {
      m_totalStates = Nc;
  }

  /// Get the total number of nodes for all processors
  CFuint getTotalNodeCount() const
  {
      return m_totalNodes;
  }

  /// Get the total number of states for all processors
  CFuint getTotalStateCount() const
  {
    return m_totalStates;
  }

  /// Set the total number of elements for all processors
  void setTotalElementCount(const std::vector<CFuint>& totalElementCount)
  {
    cf_assert(totalElementCount.size() > 0);
    m_totalElements = totalElementCount;
  }

  /// Get the total number of elements for all processors
  std::vector<CFuint> & getTotalElements()
  {
    return m_totalElements;
  }
    
  /// Get the array storing the global IDs of elements
  Common::SafePtr<std::vector<CFuint> > getGlobalElementIDs()
  {
    return &m_globalElementIDs;
  }

  /// Get the array storing the global IDs of nodes
  Common::SafePtr<std::vector<CFuint> > getGlobalNodeIDs()
  {
    return &m_globalNodeIDs;
  }
  
  /// Get the array storing the global IDs of states
  Common::SafePtr<std::vector<CFuint> > getGlobalStateIDs()
  {
    return &m_globalStateIDs;
  }

  /// Get the total number of mesh element types for all processors
  Common::SafePtr<std::vector<std::vector<std::vector<CFuint> > > >
  getGlobalTRSGeoIDs()
  {
    return &m_globalTRSGeoIDs;
  }
  
  /// Get the info about the TRS's in all processors
  std::vector<std::vector<CFuint> > & getTotalTRSInfo ()
  {
      return m_totalTRSInfo;
  }

  /// Get the names of the TRS's in all processors
  std::vector<std::string> & getTotalTRSNames ()
  {
    return m_totalTRSMap;
  }

  /// Get the TRS names
  std::vector<std::string>& getTRSNameList ()
  {
    return m_trsList;
  }

  /// Default destructor
  ~MeshData();

private: // methods

  /// Constructor
  MeshData(const std::string& name);

private: // member data
  
  /// factory registry to allow polymorphic creation of objects
  Common::SafePtr<Common::FactoryRegistry> m_fr;
  
  /// list of the names of the TRS
  MeshStatistics m_statistics;

  /// local and parallel data held and owned by the MeshData
  DataStorage* m_dataStorage;

  /// the model describing the Domain of computation
  Common::SelfRegistPtr < DomainModel > m_domainmodel;

  /// storage of the connectivities
  Common::GeneralStorage<ConnTable> m_connectivityStorage;

  /// storage of the connectivities
  Common::GeneralStorage<TopologicalRegionSet> m_trsStorage;

  /// storage of the mapping from geometric entity to
  /// TRS data
  Common::GeneralStorage<MapGeoToTrsAndIdx> m_mapGeoToTrsStorage;

  /// the vector of data of ElementTypes
  std::vector<ElementTypeData>* m_elementType;

  Common::CFMap<std::string, std::vector<ElementTypeData>* > m_groupElementTypeMap;
  
  /// number of overlap layers for parallel communication
  CFuint m_nbOverLayers;  
  
  /// flag telling if the mesh is allocated
  bool m_allocated;

  /// flag indicating if nodes and states have same connectivity
  bool m_sameNodeStateConnectivity;
  
  // Global info
  
  /// The global number of states in the mesh
  CFuint m_totalStates;

  /// The global number of nodes in the mesh
  CFuint m_totalNodes;

  /// For each element type, the global number of elements
  std::vector<CFuint> m_totalElements;

  /// global element IDs
  std::vector<CFuint> m_globalElementIDs;

  /// global node IDs
  std::vector<CFuint> m_globalNodeIDs;
  
  /// global state IDs
  std::vector<CFuint> m_globalStateIDs;
  
  /// global IDs of the GeometricEntity's in the TRS
  std::vector<std::vector<std::vector<CFuint> > > m_globalTRSGeoIDs;
  
  /// For each TRS, TR, the global number of elements
  std::vector<std::vector<CFuint> > m_totalTRSInfo;
  
  /// Name of the TRSs stored in _totalTRSInfo
  std::vector<std::string> m_totalTRSMap;

  /// socket for State's
  Framework::DataSocketSource<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for Node's
  Framework::DataSocketSource<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// list of the names of the TRS
  std::vector<std::string> m_trsList;
  
  /// string for configuring the domain model
  std::string m_domainmodel_str;
  
}; // end of class MeshData

//////////////////////////////////////////////////////////////////////////////

class Framework_API MeshDataStack : public NamespaceStack<MeshData> {
public:

  /// Returns the instance of the Active MeshData
  /// which is the one on top of the stack
  /// @return SafePtr to the active MeshData
  static Common::SafePtr<MeshData> getActive();

  /// Returns the instance of this meshDataStack
  /// This is the access point to the Singleton
  /// @return the instance of the singleton
  static MeshDataStack& getInstance();

protected: // helper functions from NamespaceStack

  /// Gets the name of the MeshData from the Namespace
  /// @param nsp the Namespace from where to get the object name
  std::string getObjectName(const Common::SafePtr<Namespace>& nsp) const;
  
  /// Creates a MeshData with the supplied name
  /// @param name of the MeshData
  MeshData * createObject(const std::string& name);

}; // end of class MeshDataStack;

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshData_hh
