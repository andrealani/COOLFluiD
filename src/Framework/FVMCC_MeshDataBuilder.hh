#ifndef COOLFluiD_Framework_FVMCC_MeshDataBuilder_hh
#define COOLFluiD_Framework_FVMCC_MeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMultiMap.hh"

#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class builds numerics dependent data inside MeshData for cell
 * center Finite Volume schemes on 1D, 2D and 3D meshes.
 *
 * @see MeshDataBuilder
 *
 * @author Andrea Lani
 */
class Framework_API FVMCC_MeshDataBuilder : public Framework::MeshDataBuilder {

public: // typedefs

  typedef Common::ConnectivityTable<CFuint> ConnTable;

public: // functions

  /**
   * Constructor
   *
   * @param name name of the builder used for configuration
   */
  FVMCC_MeshDataBuilder(const std::string& name);

  /**
   * Destructor
   */
  ~FVMCC_MeshDataBuilder();

  /**
   * Releases temporary memory used in building the mesh
   */
  virtual void releaseMemory();

  /**
   * Creates a the TopologicalRegionSet and put it in the MeshDataBuilder
   *
   * @param nbGeomEntsInTr  list of the nb of geometric entities in each tr
   * @param name  name of the TRS
   * @param trGeoCon connectivity for GeometricEntity's present in the TRS
   * @param isWritable  flag telling if the TRS can be written in output file
   */
  virtual Common::SafePtr<Framework::TopologicalRegionSet> createTopologicalRegionSet
  (const std::vector<CFuint>& nbGeomEntsInTr,
   const std::string& name,
   const Framework::TRGeoConn& trGeoCon,
   const CFuint iTRS);

protected: // functions

  /**
   * Set the max number of states in cell
   */
  void setMaxNbStatesInCell();

  /**
   * Set the max number of nodes in cell
   */
  void setMaxNbNodesInCell();

  /**
   * Set the max number of faces in cell
   */
  void setMaxNbFacesInCell();

  /**
   * Creates all the TopologicalRegionSet's from CFmeshData
   */
  void createTopologicalRegionSets();

private: // helper functions

  /**
   * Create the cell-faces connectivity
   */
  void createCellFaces();

  /**
   * Renumber local cells so that their IDs are equal to the
   * local state IDs
   * Make a consistency check on cells: verify that
   * global stateID == global cellID in each cell
   */
  void renumberCells();

  /**
   * Create the InnerFaces TRS
   */
  void createInnerFacesTRS();

  /**
   * Create the Boundary faces TRS
   */
  void createBoundaryFacesTRS();

  /**
   * Create the partition boundary faces TRS
   */
  void createPartitionBoundaryFacesTRS();

  /**
   * Create and set the mapping between faces and TRSs
   */
  void setMapGeoToTrs();

private: // data

  /// total number of faces
  CFuint m_nbFaces;

  /// geometric entity types of internal faces
  std::vector<CFuint>* m_inGeoTypes;

  /// local geo IDs of internal faces
  std::vector<CFuint>* m_inLocalGeoIDs;

  /// geometric entity types of boundary faces
  std::vector<CFuint> m_bGeoType;

  /// local geo IDs of boundary faces
  std::vector<CFuint> m_bLocalGeoIDs;

  /// face node connectivity in element type
  std::vector<Common::Table<CFuint>*> m_faceNodeElement;

  /// number of faces per element
  std::valarray<CFuint> m_nbFacesPerElem;

  /// total number of faces
  std::vector<bool> m_isBFace;

  /// geometric entity type IDs
  std::vector<CFuint> m_geoTypeIDs;

  /// number of nodes in inner faces
  std::valarray<CFuint> m_nbInFacesNodes;

  /// number of nodes in boundary faces
  std::valarray<CFuint> m_nbBFacesNodes;

  /// boundary face-state connectivity
  std::vector<CFuint> m_bFaceStateID;

  /// boundary face-node connectivity
  ConnTable* m_bFaceNodes;

  /// mapping between stateID of the face and face idx
  Common::CFMultiMap<CFuint, CFuint>  m_mapStateIdToFaceIdx;

  /// flag telling if the face is a partition face
  std::valarray<bool> m_isPartitionFace;

}; // end of class FVMCC_MeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FVMCC_MeshDataBuilder_hh
