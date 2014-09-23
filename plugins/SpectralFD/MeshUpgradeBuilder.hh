#ifndef COOLFluiD_SpectralFD_MeshUpgradeBuilder_hh
#define COOLFluiD_SpectralFD_MeshUpgradeBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/SpectralFDBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class builds data inside MeshData.
 * It assumes a Finite Volume mesh has been read and it upgrades
 * the elements to SpectralFD elements.
 * It can also upgrade the geometrical order of the mesh.
 */
class MeshUpgradeBuilder : public SpectralFD::SpectralFDBuilder {

public: // static functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

public: // functions

  /**
   * Constructor
   *
   * @param name name of the builder used for configuration
   */
  MeshUpgradeBuilder(const std::string& name);

  /**
   * Destructor
   */
  ~MeshUpgradeBuilder();

  /**
   * Releases temporary memory used in building the mesh
   */
  virtual void releaseMemory();

  /// Configures this object from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Creates all the TopologicalRegionSet's in the MeshData
   * This function overides the one in MeshDataBuilder
   */
  virtual void createTopologicalRegionSets();

  /**
   * Computes the data for the element types
   * This function overides the one in MeshDataBuilder
   */
  virtual void computeGeoTypeInfo();

protected: // functions

  /**
   * Recreate the cell-states connectivity
   * for the new SpectralFD elements
   */
  virtual void upgradeStateConnectivity();

  /**
   * Recreates the states to match the new SpectralFD elements
   */
  virtual void recreateStates();

  /**
   * Recreate cell-nodes and face-nodes connectivity
   */
  virtual void upgradeNodeConnectivity();

  /**
   * Recreates the nodes to match the new geometrically high-order elements
   */
  virtual void recreateNodes();

  /**
   * Returns the number of states in a spectral difference cell of given shape and order
   * @param geoShape the shape of the spectral difference cell
   * @param polyOrder the order of the polynomial interpolation in the spectral difference cell
   * @return the number of solution points in the given spectral difference cell type
   */
  CFuint getNbrOfStatesInSDCellType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /**
   * Returns the number of nodes in a cell of given shape and order
   * @param geoShape the shape of the cell
   * @param polyOrder the order of the geometrical polynomial interpolation in the cell
   * @return the number of nodes in the given cell type
   */
  CFuint getNbrOfNodesInCellType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /**
   * Returns the number of internal nodes in a cell of given shape and order
   * @param geoShape the shape of the cell
   * @param polyOrder the order of the geometrical polynomial interpolation in the cell
   * @return the number of internal nodes in the given cell type
   */
  CFuint getNbrOfInternalNodesInCellType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /**
   * Returns the number of nodes in a face of given shape and order
   * @param geoShape the shape of the face
   * @param polyOrder the order of the geometrical polynomial interpolation in the face
   * @return the number of nodes in the given face type
   */
  CFuint getNbrOfNodesInFaceType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /**
   * Returns the number of `internal' nodes in a face of given shape and order
   * @param geoShape the shape of the face
   * @param polyOrder the order of the geometrical polynomial interpolation in the face
   * @return the number of `internal' nodes in the given face type
   */
  CFuint getNbrOfInternalNodesInFaceType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

protected: // data

  /// spectral finite difference polynomial order
  CFPolyOrder::Type m_solPolyOrder;

  /// string holding the spectral finite difference polynomial order
  std::string m_solPolyOrderStr;

  /// geometrical polynomial order
  CFPolyOrder::Type m_geoPolyOrder;

  /// string holding the geometrical polynomial order
  std::string m_geoPolyOrderStr;

  /// previous geometrical polynomial order
  CFPolyOrder::Type m_prevGeoPolyOrder;

  /// boundary face connectivity
  std::vector< std::vector<CFuint*> > m_bndFacesNodes;

  /// total number of nodes (including new high-order nodes
  CFuint m_totNbrNodes;

};  // end of class MeshUpgradeBuilder

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_SpectralFD_MeshUpgradeBuilder_hh

