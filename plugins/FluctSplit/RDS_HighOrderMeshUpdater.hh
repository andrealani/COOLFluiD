#ifndef COOLFluiD_Numerics_FluctSplit_RDS_HighOrderMeshUpdater_hh
#define COOLFluiD_Numerics_FluctSplit_RDS_HighOrderMeshUpdater_hh

//////////////////////////////////////////////////////////////////////////////

#include "RDS_MeshDataBuilder.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//   namespace Framework { class Node; }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class builds numerics dependent data inside MeshData for RDS schemes
/// on 2D and 3D meshes.
/// Updating meshes with P1P1 elements to PmPn meshes
/// @see MeshDataBuilder
/// @author Thomas Wuilbaut
class FluctSplit_API RDS_HighOrderMeshUpdater : public RDS_MeshDataBuilder {

public: // typedefs

  typedef std::vector<CFuint*> BoundaryFaceConnect;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  /// @param name name of the builder used for configuration
  RDS_HighOrderMeshUpdater(const std::string& name);

  /// Destructor
  ~RDS_HighOrderMeshUpdater();

  /// Configures this object from the supplied arguments
  /// @param args configuration parameters
  virtual void configure ( Config::ConfigArgs& args );

  /// Releases temporary memory used in building the mesh
  virtual void releaseMemory();

  /// Creates all the TopologicalRegionSet's in the MeshData
  /// This function overides the one in MeshDataBuilder
  virtual void createTopologicalRegionSets();

  /// Computes the data for the element types
  /// This function overides the one in MeshDataBuilder
  virtual void computeGeoTypeInfo();

protected: // functions

  /// Recreate the cell-nodes connectivity

  virtual void upgradeNodalConnectivity();

  /// Recreate the cell-states connectivity
  /// for the new SpectralFV elements
  virtual void upgradeStateConnectivity();

  /// Recreates the nodes to match the new SpectralFV elements
  virtual void recreateNodes();


  /// Recreates the states to match the new SpectralFV elements
  virtual void recreateStates();


private: // functions

  /// Update the TRS Data with the new states
  void updateTRSData();

  /// Sets the nodes of the elements that are not the vertices of the elements
  /// @param nextNodeID nodeID of the first new node
  void setExtraNodes(CFuint& nextNodeID);


  /// Sets the non nodal states of the elements
  /// @param nextStateID stateID of the first new state
  void setExtraStates(CFuint& nextStateID);

  /// Returns the number of degrees of freedom in an element of given shape and order
  /// @param geoShape the shape of the element
  /// @param polyOrder the order of the polynomial interpolation in the element
  /// @return the number of DOF in the given element type
  CFuint getNbDofInNewType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /// Returns the number of interior degrees of freedom in an element of given shape and order
  /// @param geoShape the shape of the element
  /// @param polyOrder the order of the polynomial interpolation in the element
  /// @return the number of DOF in the given element type
  CFuint getNbInteriorDofInNewType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /// Returns the number of states in a boundary element
  /// @param oldNbStates the number of states of the original boundary face
  /// @param polyOrder the order of the polynomial interpolation in the face
  /// @return the number of states in the given face
  CFuint getNewNbDofInBoundaryGeo(CFuint oldNbDof, CFPolyOrder::Type polyOrder);


  /// Vector to store the number of P1 nodes per cell

  std::vector<int> nbP1CellNodes;

  /// Vector to store the number of P1 states per cell

  std::vector<int> nbP1CellStates;

  /// Vector to store the number of P1 nodes per face
//    std::vector<int> nbP1FaceNodes;


protected: // data

  ///Information about local connectivity for all element types in the mesh
  ///Information is stored twice: first for P1 elements
  ///Then for Pn elements, where n is the order to which we want to update to
  ///see functions setExtraNodes and setExtraStates

  std::vector<Table<CFuint>*> m_faceNodeP1Element;
  std::vector<Table<CFuint>*> m_faceNodePnElement;

  std::vector<Table<CFuint>*> m_faceStateP1Element;
  std::vector<Table<CFuint>*> m_faceStatePnElement;

  /// spectral finite volume polynomial order
  CFPolyOrder::Type m_newGeoPolyOrder;
  CFPolyOrder::Type m_newSolPolyOrder;

  /// Polynomial order_int is what you set with the option in CFcase file, the we convert it to PolyOrder
  CFuint m_newGeoPolyOrder_int;
  CFuint m_newSolPolyOrder_int;


  /// Should geometry/solution space be updated?
  bool m_updateGeometry;
  bool m_updateSolution;

  /// string holding the spectral finite volume polynomial order
  std::string m_newPolyOrderStr;

  std::vector<BoundaryFaceConnect> _boundaryFacesNodes;
  std::vector<BoundaryFaceConnect> _boundaryFacesStates;

  CFuint _totalNewNbNodes;
  CFuint _totalNewNbStates;

  /// nodal coordinates
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

}; // end of class RDS_HighOrderMeshUpdater

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RDS_HighOrderMeshUpdater_hh
