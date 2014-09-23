#ifndef COOLFluiD_Muffin_MuffinMeshDataBuilder_hh
#define COOLFluiD_Muffin_MuffinMeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

/// This class builds Muffin solver data inside MeshData
/// @see MeshDataBuilder
class MuffinMeshDataBuilder : public Framework::MeshDataBuilder {

public: // functions

  /// Constructor
  MuffinMeshDataBuilder(const std::string& name);

  /// Destructor
  ~MuffinMeshDataBuilder();

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

};  // end of class MuffinMeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

  }  // end of namespace Muffin
}  // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Muffin_MuffinMeshDataBuilder_hh

