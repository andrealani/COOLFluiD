// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_TopologicalRegion_hh
#define COOLFluiD_Framework_TopologicalRegion_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/ConnectivityTable.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the concept of TopologicalRegion
/// It stores pointers to the lists of local and global IDs
/// for the corresponding GeometricEntity's. These lists are
/// aggregated by the corresponding TopologicalRegionSet.
/// There is also acquaintance of the GeometricEntity IDs to State IDs
/// and GeometricEntity IDs to Node IDs connectivities
/// Also know as the acronym TR.
/// @author Andrea Lani
/// @author Tiago quintino
class Framework_API TopologicalRegion {
public:

  /// Constructor
  TopologicalRegion();

  /// Destructor
  ~TopologicalRegion();

  /// Get the number of nodes in the given GeometricEntity
  CFuint getNbNodesInGeo(CFuint iGeo) const
  {
    return m_geo2nodes->nbCols(m_startGeoIdx + iGeo);
  }

  /// Get the number of states in the given GeometricEntity
  CFuint getNbStatesInGeo(CFuint iGeo) const
  {
    return m_geo2states->nbCols(m_startGeoIdx + iGeo);
  }

  /// Get the nodeID
  CFuint getNodeID(CFuint iGeo, CFuint iNode) const
  {
    cf_assert(iNode < m_geo2nodes->nbCols(m_startGeoIdx + iGeo));
    return (*m_geo2nodes)(m_startGeoIdx + iGeo, iNode);
  }

  /// Get the stateID
  CFuint getStateID(CFuint iGeo, CFuint iState) const
  {
    cf_assert(iState < m_geo2states->nbCols(m_startGeoIdx + iGeo));
    return (*m_geo2states)(m_startGeoIdx + iGeo, iState);
  }

  /// @return the idx of GeometricEntity's in the current TR
  ///         corresponding to the sequential numbering of the
  ///         GeometricEntity's in
  ///         the currengt TRS
  /// @param local ID of the GeometricEntity in the current TR
  CFuint getGeoIDInTrs(CFuint iGeo) const
  {
    return iGeo + m_startGeoIdx;
  }

  /// Gets the local ID of the GeometricEntity iGeo
  CFuint getLocalGeoID(CFuint iGeo) const
  {
    return m_localIdx[iGeo];
  }

  /// Gets the global ID of the GeometricEntity iGeo
  CFuint getGlobalGeoID(CFuint iGeo) const
  {
    return m_globalIdx[iGeo];
  }

  /// Get the number of local geometrical entities
  CFuint getLocalNbGeoEnts() const
  {
    return m_localNbGeoEnts;
  }

  /// Sets the local number of GeoEnts in this TR.
  void setLocalNbGeoEnts(CFuint nb)
  {
    m_localNbGeoEnts = nb;
  }

  /// Gets the global number of GeoEnts in this TR.
  /// Will coincide with the number of present GeoEnts if it is not parallel.
  CFuint getGlobalNbGeoEnts() const
  {
    return m_globalNbGeoEnts;
  }

  /// Sets the global number of GeoEnts in this TR.
  /// Will coincide with the number of present GeoEnts if it is not parallel.
  void setGlobalNbGeoEnts(CFuint nb)
  {
    m_globalNbGeoEnts = nb;
  }

  /// Set the connectivity GeometricEntity to Node's
  void setGeo2NodesConn
  (Common::SafePtr<Common::ConnectivityTable<CFuint> > geo2nodes)
  {
    m_geo2nodes = geo2nodes;
  }

  /// Set the connectivity GeometricEntity to State's
  void setGeo2StatesConn
  (Common::SafePtr<Common::ConnectivityTable<CFuint> > geo2states)
  {
    m_geo2states = geo2states;
  }

  /// Puts the states in TR in the supplied vector
  /// Vector might be filled, that state pointers
  /// already present will no be removed.
  /// @param statesInTR the vector where to put the states
  /// @note the algorithm implemented could be improved
  void putStatesInTR(std::vector<CFuint>& statesInTR) const;

  /// Puts the nodes in TR in the supplied vector
  /// Vector might be filled, that node pointers already
  /// present will no be removed.
  /// @param nodesInTR the vector where to put the nodes
  /// @note the algorithm implemented could be improved
  void putNodesInTR(std::vector<CFuint>& nodesInTR) const;

  /// Sets the global index list of checking the GeometricEntity's in this TR
  /// By default_globalIdx == CFNULL
  /// @pre setGlobalNbGeoEnts() must be called first to set global number
  ///                           of GeoEnts
  /// @param geoEntIdx the list of global indexes
  /// @post return CFNULL if the indexes match the position of
  /// the GeoEnts in the GeomEntList
  void setGeoEntsGlobalIdx(CFuint *const geoEntIdx)
  {
    m_globalIdx = geoEntIdx;
  }

  /// Sets the list of the local indexes ofthe GeometricEntity's in this TR
  /// @pre setLocalNbGeoEnts() must be called first to set local number
  ///                           of GeoEnts
  /// @param geoEntIdx  the list of local indexes
  /// @param startGeoIdx   the starting idx of GeometricEntity's in the current TR
  ///                   in the sequential numbering of the TRS
  void setGeoEntsLocalIdx(CFuint *const geoEntIdx, CFuint startGeoIdx)
  {
    m_localIdx = geoEntIdx;
    m_startGeoIdx = startGeoIdx;
  }

private: // data

  /// start idx for the GeometricEntity's in this TR starting from
  /// 0 (first geo in first TR of this TRS)
  CFuint m_startGeoIdx;

  /// the local number of GeometricEntity's in this TR
  CFuint m_localNbGeoEnts;

  /// the global number of GeometricEntity's in this TR
  CFuint m_globalNbGeoEnts;

  /// pointer to the beginning of the local geo indexes
  /// storage in TRS corresponding  to this TR
  CFuint* m_localIdx;

  /// pointer to the beginning of the global geo indexes
  /// storage in TRS corresponding  to this TR
  CFuint* m_globalIdx;

  /// geo-states connectivity
  Common::SafePtr<Common::ConnectivityTable<CFuint> > m_geo2states;

  /// geo-nodes connectivity
  Common::SafePtr<Common::ConnectivityTable<CFuint> > m_geo2nodes;

}; // end of class TopologicalRegion

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_TopologicalRegion_hh
