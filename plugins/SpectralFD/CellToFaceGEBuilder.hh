#ifndef COOLFluiD_Framework_CellToFaceGEBuilder_hh
#define COOLFluiD_Framework_CellToFaceGEBuilder_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/NonCopyable.hh"

#include "Framework/CFSide.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/Storage.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class TopologicalRegionSet; }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class builds Cell's. Each Cell has pointers to the neighbouring Face's.
 *
 * @author Kris Van Den Abeele
 */
class CellToFaceGEBuilder {
public:

  /**
   * This nested struct represents groups the data needed by this builder.
   * It must be non copyable to force the client code to use by reference
   * the data aggregated by this CellToFaceGEBuilder
   *
   * @see TopologicalRegionSet
   *
   * @author Tiago Quintino
   */
  struct GeoData : public Common::NonCopyable<CellToFaceGEBuilder> {

    /**
     * Default constructor
     */
    GeoData() {}

    /// pointer to inner cells TRS
    Common::SafePtr< Framework::TopologicalRegionSet > trs;

    /// geo index in inner cells TRS
    CFuint idx;

  }; // end GeoData

  /**
   * Constructor
   */
  CellToFaceGEBuilder();

  /**
   * Destructor
   */
  ~CellToFaceGEBuilder();

  /**
   * Set up the pool
   */
  void setup();

  /**
   * Get the data of the GeometricEntity builder.
   * This allows the client code to set the data and then
   * let the builder work on its own updated data.
   */
  CellToFaceGEBuilder::GeoData& getDataGE()
  {
    return m_data;
  }

  /**
   * @return m_isFaceOnBoundary
   */
  Common::SafePtr< std::vector< bool > > getIsFaceOnBoundary()
  {
    return &m_isFaceOnBoundary;
  }

  /**
   * @return m_neighbrCellSide
   */
  Common::SafePtr< std::vector< CFuint > > getNeighbrCellSide()
  {
    return &m_neighbrCellSide;
  }

  /**
   * @return m_currentCellSide
   */
  Common::SafePtr< std::vector< CFuint > > getCurrentCellSide()
  {
    return &m_currentCellSide;
  }

  /**
   * @return m_faceOrient
   */
  Common::SafePtr< std::vector< CFuint > > getFaceOrient()
  {
    return &m_faceOrient;
  }

  /**
   * @return m_faceBCIdx
   */
  Common::SafePtr< std::vector< CFuint > > getFaceBCIdx()
  {
    return &m_faceBCIdx;
  }

  /**
   * Build the GeometricEntity corresponding to the given local ID
   * in the correspondng TopologicalRegionSet
   */
  Framework::GeometricEntity* buildGE();

  /**
   * Release the build GeometricEntity's to make them again
   * available for creation of new ones
   */
  void releaseGE();

private: // helper functions

  void assembleNeighbourFace(const CFuint cellGeoType, const CFuint faceLocalID);

private: // data

  /// data of this builder
  CellToFaceGEBuilder::GeoData  m_data;

  /// pointer to the connectivity cells-states
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_cellToStates;

  /// pointer to the connectivity cells-nodes
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_cellToNodes;

  /// pointer to the connectivity cells-faces
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_cellToFaces;

  /// pointer to the connectivity internal faces-cells
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_intFaceToCells;

  /// pointer to the connectivity internal faces-nodes
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_intFaceToNodes;

  /// pointer to the connectivity faces-internal face idx and orientation or external face boundary condition type
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_facesToIntFaceIdxAndOrientOrBCType;

  /// handle to the State's storage
  Framework::DataHandle< Framework::State*,Framework::GLOBAL> m_states;

  /// handle to the Node's storage
  Framework::DataHandle< Framework::Node*,Framework::GLOBAL> m_nodes;

  /// GeometricEntity's pool for Cells ordered by GeoType
  std::map<CFuint,
           std::vector<Framework::GeometricEntity*>
          > m_poolCells;

  /// GeometricEntity's pool for Faces ordered by GeoType
  /// Each cell has multiple faces, need to create as many as the maximum number of faces to a cell.
  std::map<CFuint,
           std::vector<Framework::GeometricEntity*>
          > m_poolFaces;

  /// face-node local indexes in a cell
  std::map<CFuint,
           Common::Table< CFuint >*
          > m_cellFaceNodes;

  /// storage of geotype of cell-faces
  std::map<CFuint,
           std::vector<CFuint>
          > m_cellFaceGeoTypes;

  /// the built GeometricEntity of CFGeoEnt::Type CELL
  std::vector< Framework::GeometricEntity* > m_builtGeoCells;

  /// list of the built GeometricEntity's of CFGeoEnt::Type FACE
  std::vector<Framework::GeometricEntity*> m_builtGeoFaces;

  /// temporary cell for working
  Framework::GeometricEntity * m_tmpCell;

  /// booleans telling whether a face is boundary
  std::vector< bool > m_isFaceOnBoundary;

  /// side of the neigbouring cell with respect to the faces
  std::vector< CFuint > m_neighbrCellSide;

  /// side of the current cell with respect to the faces
  std::vector< CFuint > m_currentCellSide;

  /// orientation of the faces
  std::vector< CFuint > m_faceOrient;

  /// boundary condition index for the boundary faces
  std::vector< CFuint > m_faceBCIdx;

}; // end of class CellToFaceGEBuilder

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CellToFaceGEBuilder_hh
