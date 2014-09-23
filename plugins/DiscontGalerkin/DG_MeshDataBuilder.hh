#ifndef COOLFluiD_DiscontGalerkin_DG_MeshDataBuilder_hh
#define COOLFluiD_DiscontGalerkin_DG_MeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class builds numerics dependent data inside MeshData for FEM schemes
 * on 2D and 3D meshes.
 *
 * @see MeshDataBuilder
 *
 * @author Tiago Quintino
 * @author Martin Holik
 * @author Vit Dolejsi
 *
 */
class DG_MeshDataBuilder : public Framework::MeshDataBuilder {

public: // typedefs

  typedef Common::ConnectivityTable<CFuint> ConnTable;

public: // functions

  /**
   * Constructor
   *
   * @param name name of the builder used for configuration
   */
  DG_MeshDataBuilder(const std::string& name);

  /**
   * Destructor
   */
  ~DG_MeshDataBuilder();

  /**
   * Releases temporary memory used in building the mesh
   */
  virtual void releaseMemory();


protected: // functions

  /**
   * Creates all the TopologicalRegionSet's in the MeshData
   */
  virtual void createTopologicalRegionSets();

  /**
   * Set the max number of states in cell
   */
  virtual void setMaxNbStatesInCell();

  /**
   * Set the max number of nodes in cell
   */
  virtual void setMaxNbNodesInCell();

  /**
   * Set the max number of faces in cell
   */
  virtual void setMaxNbFacesInCell();

  /**
   * Create and set the mapping between faces and TRSs
   */
  virtual void setMapGeoToTrs();

private: // nested helper class

  class LessThan
  {
  public:

    bool operator() (const std::pair<CFuint,CFuint>& p1,
                     const std::pair<CFuint,CFuint>& p2) const
    {
      return (p1.first < p2.first) ? true : false;
    }

  }; // end class LessThan

  class Equal
  {
  public:

    bool operator() (const std::pair<CFuint,CFuint>& p1,
                     const std::pair<CFuint,CFuint>& p2) const
    {
      return p1.first == p2.first;
    }

  }; // end class LessThan

private: // helper functions

  /**
   * Create the cell-faces connectivity
   */
  void createCellFaces();

  /**
   * Create the InnerFaces TRS
   */
  void createInnerFacesTRS();


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

  /// face to node local connectivity for each element type
  std::vector<Common::Table<CFuint>*> m_localFace2Node;

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

  /// boundary face-cell connectivity
  std::vector<CFint> m_bFaceCell;

  /// boundary face-node connectivity
  ConnTable * m_bFaceNodes;

}; // end of class DG_MeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

    } // namespace DiscontGalerkin

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_DG_MeshDataBuilder_hh
