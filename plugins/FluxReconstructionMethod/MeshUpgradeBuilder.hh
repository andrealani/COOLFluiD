#ifndef COOLFluiD_FluxReconstructionMethod_MeshUpgradeBuilder_hh
#define COOLFluiD_FluxReconstructionMethod_MeshUpgradeBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/FluxReconstructionBuilder.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class builds data inside MeshData.
 * It assumes a Finite Volume mesh has been read and it upgrades
 * the elements to FluxReconstructionMethod elements.
 * It can also upgrade the geometrical order of the mesh.
 */
class MeshUpgradeBuilder : public FluxReconstructionMethod::FluxReconstructionBuilder {

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
   * for the new FluxReconstruction elements
   */
  virtual void upgradeStateConnectivity();

  /**
   * Recreates the states to match the new FluxReconstruction elements
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
   * Divide the elements in equal parts to form new elements
   */
  virtual void divideElements();

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
  
  /**
   * Get the mapped coordinates of the new nodes for dividing the elements
   */
  std::vector< RealVector > getNewNodesMappedCoords(CFGeoShape::Type shape,CFuint solOrder, CFuint cellIdx,  const Framework::MeshData::ConnTable nodesConn);

protected: // data

  /// spectral finite difference polynomial order
  CFPolyOrder::Type m_solPolyOrder;

  /// string holding the spectral finite difference polynomial order
  std::string m_solPolyOrderStr;
  
  /// number of times the cells will be divided in equal parts to form new cells
  CFuint m_elementDivision;

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
  
  /// stores whether an old state is updatable
  std::vector< bool > m_updatables;
  
  /// stores the globalIDs of the old states
  std::vector< CFuint > m_globalIDs;
  
  /// stores the localID of the element of a new state
  std::vector< CFuint > m_elemIDOfState;
  
  /// stores the local ID of state within its element for a new state
  std::vector< CFuint > m_elemLocalIDOfState;
  
  /// stores the old local ID of the first state in an element
  std::vector< CFuint > m_elemFirstStateLocalID;
  
  /// connectivity pattern
  std::valarray<CFuint>  m_pattern;
  

};  // end of class MeshUpgradeBuilder

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_MeshUpgradeBuilder_hh

