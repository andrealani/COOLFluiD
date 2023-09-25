// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_TopologicalRegionSet_hh
#define COOLFluiD_Framework_TopologicalRegionSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NamedObject.hh"
#include "Common/TaggedObject.hh"

#include "Framework/TopologicalRegion.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a topological region set, a collection of
/// super patches.
/// It stores a list of TopologicalRegion's, local and global index
/// lists and connectivity tables (by acquaintaince) for the geo-state
/// and geo-node connectivity.
/// TopologicalRegionSet provides a unified point of access and indexing
/// of the full underlying storage (global view of data related to the
/// aggregated TopologicalRegion's).
/// Also known as the acronym TRS.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API TopologicalRegionSet :
    public Common::NamedObject,
    public Common::TaggedObject {

public:

  /// Constructor
  /// @param nameTRS enum type holding the name of the TopologicalRegionSet
  /// @param listTR  list of the topological regions
  /// @pre  listTR != CFNULL
  TopologicalRegionSet (const std::string& nameTRS, std::vector<TopologicalRegion*>* listTR);

  /// Destructor
  ~TopologicalRegionSet();

  /// Get the number of nodes in the given GeometricEntity
  CFuint getNbNodesInGeo(CFuint iGeo) const { return m_geo2nodes->nbCols(iGeo); }

  /// Get the number of states in the given GeometricEntity
  CFuint getNbStatesInGeo(CFuint iGeo) const { return m_geo2states->nbCols(iGeo); }

  /// Get the nodeID
  CFuint getNodeID(CFuint iGeo, CFuint iNode) const
  {
    cf_assert(iNode < m_geo2nodes->nbCols(iGeo));
    return (*m_geo2nodes)(iGeo,iNode);
  }

  /// Get the stateID
  CFuint getStateID(CFuint iGeo, CFuint iState) const
  {
    cf_assert(iState < m_geo2states->nbCols(iGeo));
    return (*m_geo2states)(iGeo,iState);
  }

  /// Set the nodeID
  void setNodeID ( CFuint iGeo, CFuint iNode, CFuint nodeID)
  {
    cf_assert(iNode < m_geo2nodes->nbCols(iGeo));
    (*m_geo2nodes)(iGeo,iNode) = nodeID;
  }

  /// Set the stateID
  void setStateID ( CFuint iGeo, CFuint iState, CFuint stateID)
  {
    cf_assert(iState < m_geo2states->nbCols(iGeo));
    (*m_geo2states)(iGeo,iState) = stateID;
  }

  /// Get the number of local geometrical entities
  CFuint getLocalNbGeoEnts() const {  return m_geoEntLocalIdx->size();  }

  /// Get the number of global geometrical entities
  CFuint getGlobalNbGeoEnts() const  {  return m_globalNbGeoEnts; }

  /// Get the geometric type encapsulating the knowledge
  /// of the geometrical and solution-related shape of the
  /// given GeometricEntity
  CFuint getGeoType(CFuint iGeo) const
  {
    cf_assert(iGeo < m_geoTypes->size());
    return (*m_geoTypes)[iGeo];
  }

  /// Sets the geometric entity types for all the GeometricEntities in this TRS
  void setGeoTypes(std::vector<CFuint> *const geoTypes)  
  { deletePtr(m_geoTypes); m_geoTypes = geoTypes; }
  
  /// Sets the global number of GeoEnts in this TR.
  /// Will coincide with the number of present GeoEnts if it is not parallel.
  void setGlobalNbGeoEnts(CFuint nb) {  m_globalNbGeoEnts = nb; }

  /// Get the connectivity GeometricEntity to Node's
  Common::SafePtr<Common::ConnectivityTable<CFuint> > getGeo2NodesConn() const { return m_geo2nodes; }

  /// Get the connectivity GeometricEntity to State's
  Common::SafePtr<Common::ConnectivityTable<CFuint> > getGeo2StatesConn() const { return m_geo2states; }

  /// Set the connectivity GeometricEntity to Node's
  void setGeo2NodesConn(Common::SafePtr<Common::ConnectivityTable<CFuint> > geo2nodes);

  /// Set the connectivity GeometricEntity to State's
  void setGeo2StatesConn(Common::SafePtr<Common::ConnectivityTable<CFuint> > geo2states);

  /// Create the list of states
  void createStatesList();

  /// Create the list of nodes
  void createNodesList();

  /// Gets the list of TRs
  /// @return list of TRs
  Common::SafePtr<std::vector<TopologicalRegion*> >
  getTopologicalRegionList() const
  {
    cf_assert(m_listTR != CFNULL);
    return m_listTR;
  }

  /// Gets the number of the TRs
  /// @return the number of the TRs
  CFuint getNbTRs() const
  {
    cf_assert(m_listTR != CFNULL);
    return m_listTR->size();
  }

  /// Gets the pointer to a specified TR
  /// @parameter iTR  ID of the requested TR
  /// @return the pointer to the specified TR
  Common::SafePtr<TopologicalRegion> operator[] (CFuint iTR) const  {  return getTopologicalRegion(iTR); }

  /// Gets the pointer to a specified TR
  /// @parameter iTR  ID of the requested TR
  /// @return the pointer to the specified TR
  Common::SafePtr<TopologicalRegion> getTopologicalRegion(CFuint iTR) const
  {
    cf_assert(m_listTR != CFNULL);
    cf_assert(iTR < getNbTRs());

    TopologicalRegion* ptr = (*m_listTR)[iTR];
    cf_assert(ptr != CFNULL);
    return ptr;
  }

  /// Gets the local ID of the GeometricEntity iGeo in TR iTR
  CFuint getLocalGeoID(CFuint iGeo) const {  return (*m_geoEntLocalIdx)[iGeo]; }

  /// Gets the global ID of the GeometricEntity iGeo
  CFuint getGlobalGeoID(CFuint iGeo) const  {  return (*m_geoEntGlobalIdx)[iGeo];  }

  /// Sets the global ID of the GeometricEntity iGeo
  void setGlobalGeoID(CFuint iGeo, CFuint globalGeoID)
  {
    cf_assert(iGeo < m_geoEntGlobalIdx->size());
    (*m_geoEntGlobalIdx)[iGeo] = globalGeoID;
  }

  /// Gets the global index list of the GeometricEntity's in this TRS
  /// @post return CFNULL if the indexes match the position of the
  /// GeoEnts in the GeomEntList
  Common::SafePtr<std::vector<CFuint> > getGeoEntsGlobalIdx()  { return m_geoEntGlobalIdx; }

  /// Get the list of geometric type encapsulating the knowledge
  /// of the geometrical and solution-related shape of the
  /// of the GeometricEntity's in this TRS
  Common::SafePtr<std::vector<CFuint> > getGeoTypes()  { return m_geoTypes; }

  /// Sets the global index list of checking the GeometricEntity's in this TRS
  /// @pre m_geoEntGlobalIdx == CFNULL
  /// @pre setGlobalNbGeoEnts() must be called first to set global number
  ///                           of GeoEnts
  /// @param geoEntIdx the list of global indexes
  /// @post return CFNULL if the indexes match the position of
  /// the GeoEnts in the GeomEntList
  void setGeoEntsGlobalIdx(std::vector<CFuint> *const geoEntIdx) 
  { deletePtr(m_geoEntGlobalIdx); m_geoEntGlobalIdx = geoEntIdx; }

  /// Gets the local index list of the GeometricEntity's in this TRS
  /// @post return CFNULL if the indexes match the position of the
  /// GeoEnts in the GeomEntList
  Common::SafePtr<std::vector<CFuint> > getGeoEntsLocalIdx() { return m_geoEntLocalIdx; }

  /// Sets the local index list of checking the GeometricEntity's in this TRS
  /// @pre m_geoEntLocalIdx == CFNULL
  /// @pre setLocalNbGeoEnts() must be called first to set local number
  ///                           of GeoEnts
  /// @param geoEntIdx the list of local indexes
  /// @post return CFNULL if the indexes match the position of
  /// the GeoEnts in the GeomEntList
  void setGeoEntsLocalIdx(std::vector<CFuint> *const geoEntIdx) 
  { deletePtr(m_geoEntLocalIdx); m_geoEntLocalIdx = geoEntIdx; }
  
  /// Get the list of all the states in this topological region set.
  /// @return the set of all the states (without duplications) in this TRS
  Common::SafePtr<std::vector<CFuint> > getStatesInTrs()
  {
    if (!m_hasStatesList) {
      createStatesList();
    }
    return &m_statesList;
  }

  /// Get the number of states in the current TRS
  CFuint getNbStatesInTrs() const
  {
    cf_assert(m_hasStatesList);
    return m_statesList.size();
  }

  /// Get the list of all the nodes in this topological region set.
  /// @return the set of all the nodes (without duplications) in this TRS
  Common::SafePtr<std::vector<CFuint> > getNodesInTrs()
  {
    if (!m_hasNodesList) {
      createNodesList();
    }
    return &m_nodesList;
  }

  /// Get the number of nodes in the current TRS
  CFuint getNbNodesInTrs()
  {
    if (!m_hasNodesList) {
      createNodesList();
    }
    return m_nodesList.size();
  }

private: //  data

  /// Copy constructor is declared private and is not defined to prevent
  /// default copy construction and assignment (this follows the Scott
  /// Meyers guideline in "Effective C++")
  TopologicalRegionSet(const TopologicalRegionSet& init);

  /// List of TopologicalRegion's belonging to this topological
  std::vector<TopologicalRegion*>* m_listTR;

  /// Map of the local indexes of the GeometricEntities in all TRs
  /// They are stored here and numbered sequentially looping over the TRs
  std::vector<CFuint>* m_geoEntLocalIdx;

  /// Map of the global indexes of the GeometricEntities in all TRs
  /// They are stored here and numbered sequentially looping over the TRs
  std::vector<CFuint>* m_geoEntGlobalIdx;

  /// global number og GeometricEntity's in this TRS
  CFuint m_globalNbGeoEnts;

  /// array storing the geometrical types of all the GeometricEntities in this TRS
  ///@todo use SharedPtr here
  std::vector<CFuint>* m_geoTypes;

  /// geo-states connectivity
  Common::SafePtr<Common::ConnectivityTable<CFuint> > m_geo2states;

  /// geo-nodes connectivity
  Common::SafePtr<Common::ConnectivityTable<CFuint> > m_geo2nodes;

  /// List of the states in all the TopologicalRegion's belonging
  /// to this topological region
  std::vector<CFuint> m_statesList;

  /// List of the nodes in all the TopologicalRegion's belonging
  /// to this topological region
  std::vector<CFuint> m_nodesList;

  /// flag telling if this TRS has a list of States
  bool m_hasStatesList;

  /// flag telling if this TRS has a list of Nodes
  bool m_hasNodesList;

}; // end of class TopologicalRegionSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_TopologicalRegionSet_hh
