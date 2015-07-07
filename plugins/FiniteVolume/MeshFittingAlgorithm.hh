#ifndef COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh
#define COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "SimpleEdgeGraph.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace Numerics {
    
    namespace FiniteVolume {
      class CellCenterFVMData;
      
//////////////////////////////////////////////////////////////////////////////

//Storage of truncation data for the springs
struct SpringTrucationData{
  CFreal minLimit;
  CFreal maxLimit;
  CFreal mean;
};

//////////////////////////////////////////////////////////////////////////////


/**
 *
 * This class implements a shock-fitting mesh algorithm to
 * redistribute mesh points close to shocks
 *
 * @author Pedro Santos
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
   * Executes the mesh fitting
   */
  virtual void execute();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure(Config::ConfigArgs& args);

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
   * Creates the nodal connectivity
   */
  void createNodalConnectivity();

  /**
   * Find boundary nodes
   */
  void findBoundaryNodes();

  /**
   * Compute the normals for the nodes moving in the boundary
   */
  void computeMovingInBoundaryNodeNormals(); 

  /**
   * create a mapping from the moving in boundary node ID node to faceID
   */
  std::vector<std::set<CFuint> > getMapMovingInBoundaryNodeID2FaceID();
 
  /**
   * Compute the spring max, min and average spring stiffness
   */
  void computeSpringTruncationData();
 
    /**
   * computes the spring constant between two nodes
   */
  CFreal computeSpringConstant(const Framework::Node* const firstNode, 
                               const Framework::Node* const secondNode) ;
  
  /**
   * truncates the spring constant
   */
  CFreal truncateSpringConstant(CFreal springConstant); 


  /**
   * Solve Linear System
   */
  void solveLinearSystem();
 
  /**
   * resize linear system matrix and rhs to solve for the states
   */
  void resizeSystemSolverToStateData();

  /**
   * resize linear system matrix and rhs to solve for the nodes 
   */
  void resizeSystemSolverToNodalData();

  /**
   * assemble the linear system
   */
  void assembleLinearSystem();
 
  /**
  * Queries if the node moves in the boundary
  */
  bool isNodeMovingInBoundary(Framework::Node* node);

  /**
   * assemble the system lines for a node moving in along the boundary
   */
  void assembleMovingInBoundaryNode(const Framework::Node* node);

  /**
  * Queries if the node is locked 
  */
  bool isNodeLocked(Framework::Node* node);
 
  /**
   * assemble the system lines for a locked Node
   */
  void assembleLockedNode(const Framework::Node* node);

  /**
   * assemble the system lines for an inner Node
   */
  void assembleInnerNode(const Framework::Node* node);
 
  /**
   * Update nodes with new positions
   */
  void updateNodePositions();
  
  /**
   * Send message to trigger re-computation of mesh information
   */
  void triggerRecomputeMeshData();
  
private: //data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the RHS
  Framework::DataSocketSink<CFreal> socket_rhs;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;
  
  // pointer to the data of the cell centered FVM method
  Common::SafePtr<Numerics::FiniteVolume::CellCenterFVMData> m_fvmccData;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceTRSBuilder;
  
  // accumulator for LSSMatrix
  // std::auto_ptr<Framework::BlockAccumulator> m_acc;
  
  /// Minimum percentile limit for spring constants
  CFreal m_minPercentile;

  /// Maximum percentile limit for spring constants
  CFreal m_maxPercentile;
  
  /// Step for mesh adaptation
  CFreal m_meshAcceleration;

  CFreal m_ratioBoundaryToInnerEquilibriumSpringLength;

  std::vector<std::string> m_unlockedBoundaryTRSs;
  
  /// Monitor variable string
  std::string m_monitorVar;
  
  /// Monitor variable ID)
  CFuint m_monitorVarID;

  ///Equilibrium spring length 
  CFreal m_equilibriumSpringLength;

  /// Set of Boundary Nodes
  std::set<Framework::Node*> m_boundaryNodes;

  /// Map node ID to normal vector
  std::map<CFuint, RealVector> m_mapNodeIDNormal;

  //Storage for the truncation data of the springs
  SpringTrucationData m_springTruncationData;
    
  //Simple edge graph 
  FiniteVolume::SimpleEdgeGraph m_edgeGraph;

}; // end of class MeshFittingAlgorithm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/MeshFittingAlgorithm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh
