#ifndef COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh
#define COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/CellTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace Numerics {
    
    namespace FiniteVolume {
      class CellCenterFVMData;
      
//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class implements a shock-fitting mesh algorithm to
 * redistribute mesh points close to shocks
 *
 * @author Francisco Moreira Huhn
 *
 */
template <typename MODEL>
class MeshFittingAlgorithm : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  MeshFittingAlgorithm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MeshFittingAlgorithm();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();

  /**
   * Execute on a set of dofs
   */
  virtual void execute();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:
  
  /**
   * Creates the "note to node" and "node to cell" connectivity
   */
  void createConnectivity();

  /**
   * Determine and register non-modifiable nodes
   */
  void determineNonModNodes();
  
  /**
   * Compute gradient of the monitor variable Phi
   */
  void computeGradientPhi();
  
  /**
   * Compute the spring constants for each edge
   */
  void computeSprings();
  
  /**
   * Fill the matrix and RHS of the system to be solved
   */
  void fillSystem();
  
  /**
   * Update nodes with new positions
   */
  void updateNodePositions();
  
  /**
   * Synchronize nodes and recompute mesh data (normals, volumes, etc.)
   */
  void syncNodesAndRecomputeMeshData();
  
  /**
   * Run the shock sensor on the monitor variable Phi
   */
  void shockSensor();
  
  /**
   * Compute Euclidean distance between two Nodes
   */
  CFreal getDistance(Framework::Node& a, Framework::Node& b);

  bool isNodeLocked(Framework::Node* node);

  
CFreal computeSpringConstant(const Framework::Node* const firstNode, const Framework::Node* const secondNode) ;

private: //data

  /// storage of shockSensor values
  Framework::DataSocketSource<CFreal> socket_shockSensor;
  
  /// storage of the gradient of the monitor variable Phi
  Framework::DataSocketSource<CFreal> socket_gradientPhi;
  
  /// storage of the pressure
  Framework::DataSocketSource<CFreal> socket_pressure;
  
  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of the face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the RHS
  Framework::DataSocketSink<CFreal> socket_rhs;
  
  /// storage of the flag for IDs of cells for which normals are outward
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;
  
  // pointer to the data of the cell centered FVM method
  Common::SafePtr<Numerics::FiniteVolume::CellCenterFVMData> m_fvmccData;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  
  /// Euler physical data
  RealVector m_physicalData;
  
  // accumulator for LSSMatrix
  // std::auto_ptr<Framework::BlockAccumulator> m_acc;
  
  /// minimum pressure ratio for shock sensor
  CFreal m_phiMinSS;
  
  /// maximum pressure ratio for shock sensor
  CFreal m_phiMaxSS;
  
  /// Minimum percentile limit for spring constants
  CFreal m_minPerc;

  /// Maximum percentile limit for spring constants
  CFreal m_maxPerc;
  
  /// Step for mesh adaptation
  CFreal m_alpha;
  
  /// Shock sensor rate
  CFuint m_ssRate;
  
  /// Mesh adaptation rate
  CFuint m_maRate;

  /// Monitor variable string
  std::string m_monitorVar;
  
  /// Monitor variable (nothing to do with Shock Sensor phi)
  CFuint m_phi;

  /// Degree from boundary for non-modifiable nodes
  CFuint m_boundaryDegreeNM;
  
  /// Stop adaptation but keep shock sensor
  CFuint m_stopAdaptationIter;
  
  /// Connectivity "node to node"
  /// connectedNodes[i] gives a set whose members
  /// are localIDs of nodes connected to node i
  std::vector< std::set<CFuint> > m_connectedNodes;

  /// Connectivity "node to cell"
  /// connectedCells[i] gives a set whose members
  /// are localIDs of states adjacent to node i
  std::vector< std::set<CFuint> > connectedCells;

  /// Connectivity "node to cell"
  /// connectedCells[i] gives a set whose members
  /// are localIDs of gstates adjacent to node i  
  std::vector< std::set<CFuint> > connectedGhosts;

  /// Connectivity "gstate to node"
  /// ghostNodes[i] gives a vector whose members
  /// are localIDs of nodes adjacent to gstate i
  std::vector< std::vector<CFuint> > ghostNodes;

  /// Set of Boundary Nodes
  std::set<CFuint> m_boundaryNodes;
    
  /// Set of Non-Modifiable Nodes
  std::set<CFuint> m_nonModifiableNodes;

  /// Visited array for boundary BFS
  std::vector<bool> visited;

  CFreal m_minLim;
  CFreal m_maxLim;
  CFreal m_avgKe;

}; // end of class MeshFittingAlgorithm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/MeshFittingAlgorithm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh
