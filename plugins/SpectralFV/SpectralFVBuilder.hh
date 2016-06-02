#ifndef COOLFluiD_SpectralFV_SpectralFVBuilder_hh
#define COOLFluiD_SpectralFV_SpectralFVBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class builds data inside MeshData
 * @see MeshDataBuilder
 */
class SpectralFVBuilder : public Framework::MeshDataBuilder {

public: // typedefs

  typedef Common::ConnectivityTable<CFuint> ConnTable;

public: // functions

  /**
   * Constructor
   *
   * @param name name of the builder used for configuration
   */
  SpectralFVBuilder(const std::string& name);

  /**
   * Destructor
   */
 virtual ~SpectralFVBuilder();

  /**
   * Releases temporary memory used in building the mesh
   */
  virtual void releaseMemory();

  /// Configures this object from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Computes the data for the element types
   * This function overides the one in MeshDataBuilder
   */
  virtual void computeGeoTypeInfo();

protected: // functions

  /// Set the max number of states in cell
  virtual void setMaxNbStatesInCell();

  /// Set the max number of nodes in cell
  virtual void setMaxNbNodesInCell();

  /// Set the max number of faces in cell
  virtual void setMaxNbFacesInCell();

  /// Create and set the mapping between faces and TRSs
  virtual void setMapGeoToTrs();

  /**
   * Creates all the TopologicalRegionSet's in the MeshData
   */
  virtual void createTopologicalRegionSets();

private: // helper functions

  /**
   * Create the cell-faces connectivity
   */
  void createCellFaces();

  /**
   * Create the InnerFaces TRS
   */
  void createInnerFacesTRS();

  /**
   * Reorder the faces in the InnerFacesTRS. Group the faces according to their connectivity to cells.
   */
  void reorderInnerFacesTRS();

  /**
   * Returns a vector containing all the possible node connectivities between two cells (LINE, TRIAG, TETRA)
   */
  std::vector< std::vector < std::vector < CFuint > > > getNodeConnPerOrientation(CFuint dim);

  /**
   * Create the Boundary faces TRS
   */
  void createBoundaryFacesTRS();

  /**
   * Creates a the TopologicalRegionSet and put it in the MeshDataBuilder
   *
   * @param nbGeomEntsInTr  list of the nb of geometric entities in each tr
   * @param name  name of the TRS
   * @param trGeoCon connectivity for GeometricEntity's present in the TRS
   * @param isWritable  flag telling if the TRS can be written in output file
   */
  Common::SafePtr<Framework::TopologicalRegionSet>
  createTopologicalRegionSet
  (const std::vector<CFuint>& nbGeomEntsInTr,
   const std::string& name,
   const Framework::TRGeoConn& trGeoCon,
   const CFuint iTRS);

  /**
   * Create the partition boundary faces TRS
   */
  void createPartitionBoundaryFacesTRS();

  /**
   * Reorder the faces in the BoundaryFacesTRS. Group the faces according to their connectivity to cells.
   */
  void reorderBoundaryFacesTRS();

  /**
   * Returns a vector containing all the face orientations (LINE, TRIAG, TETRA)
   */
  std::vector< std::vector < CFuint > > getBFaceOrientations(CFuint dim);

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

  /// boolean telling whether face is boundary face
  std::vector<bool> m_isBFace;

  /// geometric entity type IDs
  std::vector<CFuint> m_geoTypeIDs;

  /// number of nodes in inner faces
  std::valarray<CFuint> m_nbInFacesNodes;

  /// number of nodes in boundary faces
  std::valarray<CFuint> m_nbBFacesNodes;

  /// boundary face-cell connectivity
  std::vector<CFuint> m_bFaceCell;

  /// boundary face-node connectivity
  ConnTable * m_bFaceNodes;

  /// flag telling if the face is a partition face
  std::valarray<bool> m_isPartitionFace;

  /// max spectral finite volume polynomial order
  CFuint m_maxNbrStatesInCell;

};  // end of class SpectralFVBuilder

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_SpectralFV_SpectralFVBuilder_hh

