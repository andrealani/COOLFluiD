// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshStatistics_hh
#define COOLFluiD_Framework_MeshStatistics_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class stores the statistics of the mesh
/// @todo we need to create a step after the loading of the mesh
///       where all these values are update automatically such that we dont
///       depende on the concrete mesh data builder to do it.
/// @todo should move the global information here also
/// @author Tiago Quintino
class Framework_API MeshStatistics {

public: // methods

  /// Constructor
  MeshStatistics();

  /// Resets the statisitics to zero
  void reset();

  /// Get the max number of states in cell
  CFuint getMaxNbStatesInCell()
  {
    return m_maxNbStatesInCell;
  }

  /// Get the max number of nodes in cell
  CFuint getMaxNbNodesInCell()
  {
    return m_maxNbNodesInCell;
  }

  /// Get the max number of faces in cell
  CFuint getMaxNbFacesInCell()
  {
    return m_maxNbFacesInCell;
  }

  /// Set the max number of states in cell
  void setMaxNbStatesInCell(const CFuint maxNbStatesInCell);

  /// Set the max number of nodes in cell
  void setMaxNbNodesInCell(const CFuint maxNbNodesInCell);

  /// Set the max number of faces in cell
  void setMaxNbFacesInCell(const CFuint maxNbFacesInCell);

  /// Set the number of faces in the all mesh
  void setNbFaces(const CFuint nbFaces)
  {
    m_nbFaces = nbFaces;
  }

  /// Get the number of faces in the mesh
  CFuint getNbFaces()
  {
    return m_nbFaces;
  }

private: // member data

  /// the number of faces in the all mesh
  CFuint m_nbFaces;

  /// the maximum number of states in the
  CFuint m_maxNbStatesInCell;

  /// the maximum number of nodes in the element
  CFuint m_maxNbNodesInCell;

  /// the maximum number of faces in a cell
  CFuint m_maxNbFacesInCell;

}; // end of class MeshStatistics

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshStatistics_hh
