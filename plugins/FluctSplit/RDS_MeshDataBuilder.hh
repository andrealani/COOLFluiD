#ifndef COOLFluiD_Numerics_FluctSplit_RDS_MeshDataBuilder_hh
#define COOLFluiD_Numerics_FluctSplit_RDS_MeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class builds numerics dependent data inside MeshData for FEM schemes
/// on 2D and 3D meshes.
/// @see MeshDataBuilder
/// @author Tiago Quintino
/// @author Andrea Lani
class FluctSplit_API RDS_MeshDataBuilder : public Framework::MeshDataBuilder {

public: // functions

  /// Constructor
  /// @param name name of the builder used for configuration
  RDS_MeshDataBuilder(const std::string& name);

  /// Destructor
  ~RDS_MeshDataBuilder();

  /// Releases temporary memory used in building the mesh
  virtual void releaseMemory();

protected: // functions

  /// Set the max number of states in cell
  virtual void setMaxNbStatesInCell();

  /// Set the max number of nodes in cell
  virtual void setMaxNbNodesInCell();

  /// Set the max number of faces in cell
  virtual void setMaxNbFacesInCell();

  /// Create and set the mapping between faces and TRSs
  virtual void setMapGeoToTrs();

private: // data

}; // end of class RDS_MeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RDS_MeshDataBuilder_hh
