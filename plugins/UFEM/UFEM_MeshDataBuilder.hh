#ifndef COOLFluiD_UFEM_UFEM_MeshDataBuilder_hh
#define COOLFluiD_UFEM_UFEM_MeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshDataBuilder.hh"
#include "UFEM/UFEMSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/// This class builds numerics dependent data inside MeshData for UFEM schemes
/// on 2D and 3D meshes.
/// @see MeshDataBuilder
/// @author Tiago Quintino
class UFEM_API UFEM_MeshDataBuilder : public Framework::MeshDataBuilder {

public: // functions

  /// Constructor
  /// @param name name of the builder used for configuration
  UFEM_MeshDataBuilder(const std::string& name);

  /// Destructor
  ~UFEM_MeshDataBuilder();

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

}; // end of class UFEM_MeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_UFEM_MeshDataBuilder_hh
