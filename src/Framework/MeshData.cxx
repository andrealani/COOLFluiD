// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/GeneralStorage.hh"

#include "Framework/State.hh"
#include "Framework/GeometricEntityRegister.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/MeshData.hh"
#include "Framework/DomainModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void MeshData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("listTRS","List of TRS's to be built.");
  options.addConfigOption< std::string >("DomainModel","Type of domain model to describe the computational domain.");
  options.addConfigOption< NamespaceGroup::StorageType >("Namespaces","List the Namespaces which will be present in this MeshData");
  options.addConfigOption< bool >("sameNodeStateConnectivity","Option to assume the Node and State connectivity to be the same.");
}

//////////////////////////////////////////////////////////////////////////////

MeshData::MeshData(const std::string& name) :
  NamespaceGroup(name),
  ConfigObject(name),
  m_fr(CFNULL),
  m_statistics(),
  m_dataStorage(new DataStorage()),
  m_domainmodel(),
  m_connectivityStorage(),
  m_trsStorage(),
  m_mapGeoToTrsStorage(),
  m_elementType(new vector<ElementTypeData>()),
  m_groupElementTypeMap(),
  m_nbOverLayers(1),
  m_allocated(true),
  m_sameNodeStateConnectivity(false),
  m_totalStates(0),
  m_totalNodes(0),
  m_totalElements(),
  m_globalElementIDs(), 
  m_globalNodeIDs(),
  m_globalStateIDs(),
  m_globalTRSGeoIDs(),
  m_totalTRSInfo(),
  m_totalTRSMap(),
  socket_states("states"),
  socket_nodes("nodes")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  // config options
  setParameter("Namespaces",&_namespaceList);

  m_trsList = vector<std::string>();
  setParameter("listTRS",&m_trsList);

  m_sameNodeStateConnectivity = false;
  setParameter("sameNodeStateConnectivity",&m_sameNodeStateConnectivity);

  m_domainmodel_str = "Null";
  setParameter("DomainModel",&m_domainmodel_str);
}

//////////////////////////////////////////////////////////////////////////////

MeshData::~MeshData()
{
  CFAUTOTRACE;

  deallocate();
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ConfigObject::configure(args);

  socket_states.setParentNamespace(getPrimaryNamespace());
  socket_nodes.setParentNamespace(getPrimaryNamespace());

  SafePtr<DomainModel::PROVIDER> prov_dm = CFNULL;
  try {
    prov_dm = FACTORY_GET_PROVIDER(getFactoryRegistry(), DomainModel, m_domainmodel_str);
  }
  catch (Common::NoSuchValueException& e) {
    CFLog (DEBUG_MIN, e.what() << "\n" );
    prov_dm = FACTORY_GET_PROVIDER(getFactoryRegistry(), DomainModel, "Null");
  }
  cf_assert(prov_dm.isNotNull());
 
  Common::SelfRegistPtr < DomainModel >* dm = prov_dm->createPtr(m_domainmodel_str);   

  m_domainmodel = *dm;

  cf_assert(m_domainmodel.isNotNull());
  m_domainmodel->setFactoryRegistry(getFactoryRegistry());
 
  configureNested ( m_domainmodel.getPtr(), args );

  GeometricEntityRegister::getInstance().setFactoryRegistry(getFactoryRegistry());
  
  delete dm;
}

//////////////////////////////////////////////////////////////////////////////
 
Common::SafePtr< DomainModel > MeshData::getDomainModel()
{
  cf_assert(m_domainmodel.isNotNull());
  return m_domainmodel.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr< DataStorage > MeshData::getDataStorage()
{
  cf_assert(m_dataStorage != CFNULL);
  return m_dataStorage;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr< std::vector<ElementTypeData> > MeshData::getElementTypeData()
{
  cf_assert(m_elementType != CFNULL);
  return m_elementType;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr< std::vector<ElementTypeData> >
MeshData::getElementTypeData(const std::string trsName)
{
  return m_groupElementTypeMap.find(trsName);
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<TopologicalRegionSet> MeshData::getTrs(const std::string& name)
{
  return m_trsStorage.getEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr<TopologicalRegionSet> > MeshData::getTrsList()
{
  std::vector< SafePtr<TopologicalRegionSet> > all;
  std::transform(m_trsStorage.begin(),
                 m_trsStorage.end(),
                 back_inserter(all),
                 GeneralStorage<TopologicalRegionSet>::extract);
  return all;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr<TopologicalRegionSet> > MeshData::getFilteredTrsList(const std::string& tag)
{
  std::vector< SafePtr<TopologicalRegionSet> > all = getTrsList();
  std::vector< SafePtr<TopologicalRegionSet> > filtered;

  const CFuint nbTrs = all.size();
  for(CFuint iTrs=0; iTrs < nbTrs; ++iTrs)
  {
    if(all[iTrs]->hasTag(tag)) filtered.push_back(all[iTrs]);
  }

  return filtered;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr<TopologicalRegionSet> > MeshData::getFilteredTrsList(const std::vector<std::string>& tagList)
{
  std::vector< SafePtr<TopologicalRegionSet> > all = getTrsList();
  std::vector< SafePtr<TopologicalRegionSet> > filtered;

  const CFuint nbTrs = all.size();
  const CFuint nbTags = tagList.size();
  for(CFuint iTrs=0; iTrs < nbTrs; ++iTrs)
  {
    bool hasAllTags = true;
    for(CFuint iTag=0; iTag < nbTags; ++iTag)
    {
      if(!all[iTrs]->hasTag(tagList[iTag])) hasAllTags = false;
    }
    if(hasAllTags) filtered.push_back(all[iTrs]);
  }

  return filtered;
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::deallocate()
{
  deallocateSockets();
    
  if (m_allocated)
  {
    deletePtr(m_dataStorage);
    deletePtr(m_elementType);
        
    for(CFuint i=0; i<m_groupElementTypeMap.size(); ++i){
      deletePtr(m_groupElementTypeMap[i]);
    }
    m_groupElementTypeMap.clear();
    
    // the TRSs are owned by this storage
    m_trsStorage.deleteAllEntries();  
        
    // these maps are not owned by this storage
    m_mapGeoToTrsStorage.removeAllEntries();
    
    Statistics().reset();
    
    m_allocated = false;
  }
  
  /// GeometricEntityRegistry should be a member of MeshData
  /// Not a Singleton
  GeometricEntityRegister::getInstance().clear(); 
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::reallocate()
{
  CFAUTOTRACE;

  if (!m_allocated) {
    m_dataStorage = new DataStorage();
    m_elementType = new vector<ElementTypeData>();

    //reset the map
    for(CFuint i=0; i<m_groupElementTypeMap.size(); ++i){
      deletePtr(m_groupElementTypeMap[i]);
    }
    m_groupElementTypeMap.clear();

    Statistics().reset();

    m_allocated = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::createElementTypeData(const std::string trsName)
{
  if(!m_groupElementTypeMap.exists(trsName)) m_groupElementTypeMap.insert(trsName, new std::vector<ElementTypeData>);
  m_groupElementTypeMap.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::allocateSockets()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "MeshData::allocateSockets() => start\n");
  
  // AL: the namespace could be an input of this function
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  // allocate the State's and the Node's
  socket_states.allocate(this->getDataStorage(), nsp);
  socket_nodes.allocate(this->getDataStorage(), nsp);
  
  CFLog(VERBOSE, "MeshData::allocateSockets() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::deallocateSockets()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "MeshData::deallocateSockets() => start\n");
  
  // deallocate the State's and the Node's
  socket_states.deallocate();
  socket_nodes.deallocate();
  
  CFLog(VERBOSE, "MeshData::deallocateSockets() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::deallocateNodesStates()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  /// @todo there is still something strange here...
  /// If you delete first the nodes, then the states, sometimes the states
  /// will try to delete a node that has already been deleted!

  // delete the states
  for(CFuint i = 0; i < states.size(); i++)
  {
    deletePtr(states[i]);
  }

  // delete the nodes
  for(CFuint i = 0; i < nodes.size(); i++)
  {
    deletePtr(nodes[i]);
  }

  // deallocateSockets();

  /// @todo the use of IndexList has to be re-discussed !!!!!!!
  IndexList<Node>::getList().reset();
  IndexList<State>::getList().reset();
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::allocateConnectivity()
{
  CFAUTOTRACE;

  ConnectivityTable<CFuint>* cellNodes = new ConnectivityTable<CFuint>();

  // if the mesh has only isoparametric elements point to the same table
  // this provides a  big saving in memory

  ConnTable* cellStates = (m_sameNodeStateConnectivity) ? cellNodes : new ConnectivityTable<CFuint>();

  // store them by their name
  storeConnectivity("cellNodes_InnerCells", cellNodes);
  storeConnectivity("cellStates_InnerCells", cellStates);
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::deallocateConnectivity()
{
  CFAUTOTRACE;

  std::vector<std::string> tagList(2);
  tagList[0] = "inner";
  tagList[1] = "cell";

  std::vector< Common::SafePtr<TopologicalRegionSet> > innerCellsList =
    MeshDataStack::getActive()->getFilteredTrsList(tagList);

  for(CFuint iGroup=0; iGroup < innerCellsList.size(); ++iGroup){

    deleteConnectivity("cellNodes_" + innerCellsList[iGroup]->getName());

    if (m_sameNodeStateConnectivity) {
      removeConnectivity("cellStates_" + innerCellsList[iGroup]->getName());
    }
    else {
      deleteConnectivity("cellStates_" + innerCellsList[iGroup]->getName());
    }
  }
  // the connectivity tables are not owned by this storage
  m_connectivityStorage.deleteAllEntries();

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MeshData::getProvidedSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MeshData::ConnTable>
MeshData::getConnectivity(const std::string& name)
{
  return m_connectivityStorage.getEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::storeConnectivity(const std::string& name, MeshData::ConnTable
*const conn)
{
  // cleanup the entry if it already exists
  if (m_connectivityStorage.checkEntry(name)) {
    m_connectivityStorage.deleteEntry(name);
  }
  m_connectivityStorage.addEntry(name,conn);
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::removeConnectivity(const std::string& name)
{
  m_connectivityStorage.removeEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::deleteConnectivity(const std::string& name)
{
  m_connectivityStorage.deleteEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::addTrs(TopologicalRegionSet * trs)
{
  // cleanup the entry if it already exists
  if (m_trsStorage.checkEntry(trs->getName())) {
    m_trsStorage.deleteEntry(trs->getName());
  }
  m_trsStorage.addEntry(trs->getName(),trs);
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::storeMapGeoToTrs(const std::string& name,
                                MapGeoToTrsAndIdx *const mapG)
{
  // cleanup the entry if it already exists
  if (m_mapGeoToTrsStorage.checkEntry(name)) {
    m_mapGeoToTrsStorage.deleteEntry(name);
  }
  m_mapGeoToTrsStorage.addEntry(name,mapG);
}
    
//////////////////////////////////////////////////////////////////////////////

void MeshData::removeMapGeoToTrs(const std::string& name)
{
  m_mapGeoToTrsStorage.removeEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MapGeoToTrsAndIdx>
MeshData::getMapGeoToTrs(const std::string& name)
{
  return m_mapGeoToTrsStorage.getEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

MeshDataStack& MeshDataStack::getInstance()
{
  static MeshDataStack aMeshDataStack;
  return aMeshDataStack;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MeshData> MeshDataStack::getActive()
{

  cf_assert(getInstance().isEnabled());
  // returns the first on the stack
  return getInstance().top();
}

//////////////////////////////////////////////////////////////////////////////

MeshData* MeshDataStack::createObject(const std::string& name)
{
  MeshData * ptr = new MeshData(name);
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

std::string
MeshDataStack::getObjectName(const Common::SafePtr<Namespace>& nsp) const
{
  return nsp->getMeshDataName();
}

//////////////////////////////////////////////////////////////////////////////

DataSocketSink<Framework::Node*, Framework::GLOBAL> MeshData::getNodeDataSocketSink()
{
  return socket_nodes;
}

//////////////////////////////////////////////////////////////////////////////

DataSocketSink<Framework::State*, Framework::GLOBAL> MeshData::getStateDataSocketSink()
{
  return socket_states;
}

//////////////////////////////////////////////////////////////////////////////

void MeshData::setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr)
{
  m_fr = fr;
}
    
//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Common::FactoryRegistry> MeshData::getFactoryRegistry() 
{
#ifdef CF_HAVE_SINGLE_EXEC
  cf_assert(m_fr != CFNULL);
#endif
  return m_fr;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

