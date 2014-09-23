#ifndef COOLFluiD_Numerics_FiniteElement_FEM_HighOrderMeshUpdater_hh
#define COOLFluiD_Numerics_FiniteElement_FEM_HighOrderMeshUpdater_hh

//////////////////////////////////////////////////////////////////////////////

#include "FEM_MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class builds numerics dependent data inside MeshData for FEM schemes
 * on 2D and 3D meshes.
 * Updating meshes with P1P1 elements to P1PX meshes
 *
 * @see MeshDataBuilder
 *
 * @author Thomas Wuilbaut
 *
 */
class FEM_HighOrderMeshUpdater : public FEM_MeshDataBuilder {

public: // typedefs

  typedef std::pair< std::vector<CFuint*> , std::vector<CFuint*> >  BoundaryFaceConnect;

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   *
   * @param name name of the builder used for configuration
   */
  FEM_HighOrderMeshUpdater(const std::string& name);

  /**
   * Destructor
   */
  ~FEM_HighOrderMeshUpdater();

  /**
   * Releases temporary memory used in building the mesh
   */
  virtual void releaseMemory();

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
   * for the new SpectralFV elements
   */
  virtual void upgradeStateConnectivity();

  /**
   * Recreates the states to match the new SpectralFV elements
   */
  virtual void recreateStates();


private: // functions

  /**
   * Update the TRS Data with the new states
   */
  void updateTRSData();

  /**
   * Sets the non nodal states of the elements
   * @param nextStateID stateID of the first new state
   */
  void setExtraStates(CFuint& nextStateID);

  /**
   * Returns the number of states in an element of given shape and order
   * @param geoShape the shape of the element
   * @param polyOrder the order of the polynomial interpolation in the element
   * @return the number of states in the given element type
   */
  CFuint getNbStatesInNewType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /**
   * Returns the number of interior states in an element of given shape and order
   * @param geoShape the shape of the element
   * @param polyOrder the order of the polynomial interpolation in the element
   * @return the number of states in the given element type
   */
  CFuint getNbInteriorStatesInNewType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /**
   * Returns the number of states in a boundary element
   * @param oldNbStates the number of states of the original boundary face
   * @param polyOrder the order of the polynomial interpolation in the face
   * @return the number of states in the given face
   */
  CFuint getNewNbStatesInBoundaryGeo(CFuint oldNbStates, CFPolyOrder::Type polyOrder);


protected: // data

  /// spectral finite volume polynomial order
  CFPolyOrder::Type m_newPolyOrder;

  /// Polynomial orderit is the what you set with the ootin, the we convert it to PolyOrder
  CFuint m_newPolyOrder_int;

  /// string holding the spectral finite volume polynomial order
  std::string m_newPolyOrderStr;

  std::vector<BoundaryFaceConnect> _boundaryFaces;

  CFuint _totalNewNbStates;

}; // end of class FEM_HighOrderMeshUpdater

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_FEM_HighOrderMeshUpdater_hh
