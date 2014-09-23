// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FaceCellTrsGeoBuilder_hh
#define COOLFluiD_Framework_FaceCellTrsGeoBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CellTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a builder of a standard GeometricEntity,
/// with Node's and State's starting from some TRS-based data. It is
/// not supposed to be overridden.
/// This builder allows to create, use and release only one GeometricEntity
/// at a time.
/// @see GeometricEntity
/// @see Node
/// @see State
/// @see ConnectivityTable
/// @see TopologicalRegionSet
/// @see TopologicalRegion
/// @see MeshData
/// @author Andrea Lani
class Framework_API FaceCellTrsGeoBuilder : public CellTrsGeoBuilder {
public:

  /// This nested struct represents groups the data needed by this builder.
  /// It must be non copyable to force the client code to use by reference
  /// the data aggregated by this FaceCellTrsGeoBuilder
  /// @see TopologicalRegionSet
  /// @author Andrea Lani
  struct GeoData : public Common::NonCopyable<GeoData> {

    /// Default constructor
    GeoData() {}

    /// pointer to TRS of cells
    Common::SafePtr<Framework::TopologicalRegionSet> cells;
    
    /// pointer to TRS of faces
    Common::SafePtr<Framework::TopologicalRegionSet> faces;
    
    /// flag telling if the face is on the boundary
    bool isBFace;
    
    /// flag telling if all cells have to be built
    bool allCells;
    
    /// face index in face TRS
    CFuint idx;
  };
  
  /// Constructor
  FaceCellTrsGeoBuilder();

  /// Destructor
  ~FaceCellTrsGeoBuilder();

  /// Sets the cell flags DataSocket
  void setCellFlagSocket(Framework::DataSocketSink<bool> cellFlagSocket)
  {
    socket_cellFlag = cellFlagSocket;
  }
  
  /// Get the data of the GeometricEntity builder.
  /// This allows the client code to set the data and then
  /// let the builder work on its own updated data.
  FaceCellTrsGeoBuilder::GeoData& getDataGE()
  {
    return m_fcdata;
  }
  
  /// Build the GeometricEntity corresponding to the given local ID
  /// in the correspondng TopologicalRegionSet
  Framework::GeometricEntity* buildGE();
  
private:
  
  /// socket for cell flags
  Framework::DataSocketSink<bool> socket_cellFlag;
  
  /// data of this builder
  FaceCellTrsGeoBuilder::GeoData  m_fcdata;
  
}; // end of class FaceCellTrsGeoBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FaceCellTrsGeoBuilder_hh
