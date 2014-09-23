// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FaceToCellGEBuilder_hh
#define COOLFluiD_Framework_FaceToCellGEBuilder_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/NonCopyable.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/Storage.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

  class TopologicalRegionSet;

//////////////////////////////////////////////////////////////////////////////

/// This class builds the Face's used by the Framework method.
/// Each face has a pointer to the neighbor cell.
/// @author Kris Van Den Abeele
/// @author Tiago Quintino
class Framework_API FaceToCellGEBuilder {
public:

  /// This nested struct represents groups the data needed by this builder.
  /// It must be non copyable to force the client code to use by reference
  /// the data aggregated by this FaceToCellGEBuilder
  /// @see TopologicalRegionSet
  /// @author Tiago Quintino
  struct GeoData : public Common::NonCopyable<FaceToCellGEBuilder> {

    /// Default constructor
    GeoData() {}

    /// pointer to inner faces TRS
    Common::SafePtr< Framework::TopologicalRegionSet > facesTRS;

    /// pointer to inner cells TRS
    Common::SafePtr< Framework::TopologicalRegionSet > cellsTRS;

    /// geo index in inner faces TRS
    CFuint idx;

    /// flag indicating if the TRS is a boundary
    bool isBoundary;

  }; // end GeoData

  /// Constructor
  FaceToCellGEBuilder();

  /// Destructor
  ~FaceToCellGEBuilder();

  /// Set up the pool
  void setup();

  /// Get the data of the GeometricEntity builder.
  /// This allows the client code to set the data and then
  /// let the builder work on its own updated data.
  FaceToCellGEBuilder::GeoData& getDataGE() { return m_data; }

  /// Build the GeometricEntity corresponding to the given local ID
  /// in the correspondng TopologicalRegionSet
  Framework::GeometricEntity* buildGE();

  /// Release the build GeometricEntity's to make them again
  /// available for creation of new ones
  void releaseGE();

private: // helper functions

  void assembleNeighbourCell(const CFuint side);

private: // data

  /// data of this builder
  FaceToCellGEBuilder::GeoData  m_data;

  /// pointer to the connectivity faces-cells
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_FaceToCells;

  /// pointer to the connectivity cells-states
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_CellToStates;

  /// pointer to the connectivity cells-nodes
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_CellToNodes;

  /// handle to the State's storage
  Framework::DataHandle< Framework::State*, Framework::GLOBAL> m_states;

  /// handle to the Node's storage
  Framework::DataHandle< Framework::Node*, Framework::GLOBAL> m_nodes;

  /// GeometricEntity's pool for Cells ordered by GeoType
  /// every Face will have two cells which may or may not be
  /// of same type, so we need to create two cells of each type
  std::map<CFuint,
           std::pair<Framework::GeometricEntity*,Framework::GeometricEntity*>
          > m_PoolCells;

  /// GeometricEntity's pool for Faces ordered by GeoType
  std::map<CFuint,Framework::GeometricEntity*> m_PoolFaces;

  /// list of the built GeometricEntity's of CFGeoEnt::Type CFGeoEnt::CELL
  std::vector<Framework::GeometricEntity*> m_BuiltGeoCells;

  /// list of the built GeometricEntity's of CFGeoEnt::Type CFGeoEnt::FACE
  std::vector<Framework::GeometricEntity*> m_BuiltGeoFaces;

  ///flag to known if the builder has already been setup
  bool m_isSetup;

  /// temporary face for working
  Framework::GeometricEntity * m_tmpFace;


}; // end of class FaceToCellGEBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FaceToCellGEBuilder_hh
