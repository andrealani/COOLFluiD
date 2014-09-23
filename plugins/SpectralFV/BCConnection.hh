#ifndef COOLFluiD_Numerics_SpectralFV_BCConnection_hh
#define COOLFluiD_Numerics_SpectralFV_BCConnection_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/ReconstructStatesSpectralFV.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a connection boundary condition
 * @warning this will not work for implicit, since the `ghost values' are not perturbed
 *
 * @author Kris Van den Abeele
 */
class BCConnection : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCConnection(const std::string& name);

  /// Destructor
  ~BCConnection();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCConnection";
  }

  /// Set up private data and data
  void setup();

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                             std::vector< std::vector< RealVector* > >& ghostGrads,
                             const std::vector< RealVector >& normals,
                             const std::vector< RealVector >& coords);

private: // methods

  /**
   * set data related to the current boundary face
   */
  void setBndFaceData();

  /**
   * create the face to cell and face orientation data
   */
  void createFaceConnectionData();

  /**
   * @return the ID of the node that is connected to the given node
   */
  CFuint getMatchingNodeID(const CFuint givenNode);

protected: // data

  /// handle to the State's storage
  Framework::DataHandle< Framework::State*,Framework::GLOBAL> m_states;

  /// connectivity cells-states
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_cellStatesConn;

  /// Strategy that reconstructs the states in a given number of nodes
  Common::SafePtr< ReconstructStatesSpectralFV > m_statesReconstr;

  /// IDs of nodes that are connected
  std::vector< CFuint > m_connectedNodes;

  /// map of boundary face global ID's to `ghost' cell IDs and orientation
  std::map< CFuint , std::pair< CFuint , CFuint > > m_faceToGCellAndFaceOrient;

  /// flux point on SV faces reconstruction coefficients
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
      m_svFaceflxPntsRecCoefs;

  /// ghost cell states
  std::vector< Framework::State* > m_ghostCellStates;

  /// ghost cell gradients
  std::vector< std::vector< RealVector >* > m_ghostCellGrads;

  /// orientation of current face
  CFuint m_orient;

  /// name of file containing connected nodes
  std::string m_connNodesFileName;

}; // class BCConnection

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_BCConnection_hh
