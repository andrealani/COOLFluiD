#ifndef COOLFluiD_FluxReconstructionMethod_MeshFittingAlgorithmFRQ2_hh
#define COOLFluiD_FluxReconstructionMethod_MeshFittingAlgorithmFRQ2_hh


#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/TRSDistributeData.hh"

#include "SimpleEdgeGraphFRQ2.hh"

#include "Common/CFMultiMap.hh"
#include "Framework/GeometricEntityPool.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolver.hh" 
#include "FluxReconstructionMethod/FluxReconstruction.hh"
////////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }


    namespace FluxReconstructionMethod {
      class FluxReconstructionSolverData;   /// To be changed CellCenterFVMData
      
////////////////////////////////////////////////////////////////////////////////
      
      //Storage of truncation data for the springs
struct SpringTrucationData{
  CFreal minLimit;
  CFreal maxLimit;
  CFreal mean;
};

////////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class implements a r-refinement mesh algorithm to
 * redistribute mesh points close to shocks for FR
 *
 * @authors Firas Ben Ameur
 *          Joachim Balis
 *
 */
class MeshFittingAlgorithmFRQ2 : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  MeshFittingAlgorithmFRQ2(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MeshFittingAlgorithmFRQ2();

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
  void  createNodalConnectivity();


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
   * Computing the nodalState
   */
  
  CFreal computeState(const Framework::Node*  Node);
  
  /**
   * computes the spring constant between two nodes
   */
  CFreal computeSpringConstant(const Framework::Node* const firstNode, 
                               const Framework::Node* const secondNode) ;
  

  void computeNodeStates();
  /**
   * Compute the spring constant between two nodes using
   * the semi torsional Spring Analogy for 2D trianglar meshes
   */
  
  CFreal computetorsionConstant(const  Framework::Node* const  centralNode, 
				const  Framework::Node* const firstNode,
				const  Framework::Node* const  secondNode);
  
  /**
   * Compute the spring constant between two nodes using
   * the semi torsional Spring Analogy and/or the ortho-semi 
   * torisonal spring analogy for 3D tetrahedral meshes
   */
  
  CFreal computetorsionConstant3D(CFuint nodeID1, CFuint nodeID2, 
				  const  Framework::Node* const firstNode,
				  const  Framework::Node* const  secondNode);
  
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
   * Compute the spring constant between two nodes using
   * the semi torsional Spring Analogy for 2D quadrilateral meshes
   */  
  CFreal computeConstantquads(const  RealVector xyi , 
			      const  Framework::Node* const  firstNode,
			      const  Framework::Node* const  secondNode,
			      const  Framework::Node* const  thirdNode,
			      const  Framework::Node* const  fourthNode,
			      CFreal elementArea);
 
  
  /**
   * assemble the system lines for a node moving in along the boundary
   */
  void assembleMovingInBoundaryNode(const Framework::Node* node);

  /**
   * Queries if the node is locked 
   */
  bool isNodeLocked( Framework::Node* node);
  /**
   * Queries if the node is inside the region 
   */
  //bool insideRegion(Framework::Node* node);
  
  /**
   * Compute the area of a quadrilateral element 
   */
  CFreal computeElementArea2dQuads(const  Framework::Node* const  a,const  Framework::Node* const  c,
				   const  Framework::Node* const  b, const  Framework::Node* const  d);
  /**
   *  Compute the aspect ration of a quadrilateral element 
   */
  CFreal computeAspectRatio2dQuads(const  Framework::Node* const  a,const  Framework::Node* const  b,
				   const  Framework::Node* const  c, const  Framework::Node* const  d);
  
  /**
   * assemble the system lines for a locked Node
   */
  void assembleLockedNode(const Framework::Node* node);
  
  /**
   * assemble the system lines for a 2D triangular element
   */ 
  void assembleinRegionNode2DTriag(const  Framework::Node* node);
  
  /**
   * assemble the systme lines for 2D quadrilateral element
   */
  void assembleinRegionNode2DQuads(const  Framework::Node* node);
  
  /**
 * assemble system lines for 3D tetrahedral elemnt 
 */
  void assembleinRegionNode3DTet(const  Framework::Node* node);
  
  /**
   * Compute the intersection point between two quadrilateral diagonals
   */
  RealVector computeIntersection(const  Framework::Node* const  a, const  Framework::Node* const  c,
				 const  Framework::Node* const  b, const  Framework::Node* const  d);
  
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


  void createGeneralConnectivityFR();

  void findCenterNodes();
  bool isBoundaryNode( Framework::Node* node);
  bool isInSideCell( Framework::Node* node);
  void nbOfConnectedFacesToaNode();
  
private: //data
  /// the socket of the nodes connectivity for 2D quadrilatral mesh
  Common::CFMultiMap <CFuint , CFuint>  m_mapNodeNode1;
  
  /// the socket to the data handle of the nodal stiffness
  Framework::DataSocketSource < CFreal > socket_stiffness;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the RHS
  Framework::DataSocketSink<CFreal> socket_rhs;


  /// socket for the wallDistance storage
  //Framework::DataSocketSink<CFreal> socket_wallDistance;

  // pointer to the data of the cell centered FVM method
  Common::SafePtr<FluxReconstructionSolverData> m_frData;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;
  
  
  /// builder of geometric entities

  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_geoBuilder;
    /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;
  
  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  // physical data array
  RealVector m_pdata;
 
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder1;
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceTRSBuilder;

  // dummy state vector
  std::auto_ptr<Framework::State> m_state; 
  
  /// Minimum percentile limit for spring constants
  CFreal m_minPercentile;

  /// Maximum percentile limit for spring constants
  CFreal m_maxPercentile;
  
  /// Step for mesh adaptation
  CFreal m_meshAcceleration;

  /// Tolerance on the mesh mouvement for RSI

  CFreal m_ratioBoundaryToInnerEquilibriumSpringLength;

  std::vector<std::string> m_unlockedBoundaryTRSs;
  
  /// Monitor variable ID from State
  CFuint m_monitorVarID;
  
  /// Monitor variable ID from physical data
  CFuint m_monitorPhysVarID;

  /// MQIvalue to be defined 
  
  ///Equilibrium spring length 
  CFreal m_equilibriumSpringLength;

  /// Set of Boundary Nodes
  std::set< Framework::Node* > m_boundaryNodes;
  std::set< Framework::Node* > m_CenterNodes;

  /// Map node ID to normal vector
  std::map<CFuint, RealVector> m_mapNodeIDNormal;

  //Storage for the truncation data of the springs
  SpringTrucationData m_springTruncationData;
    

  /// number of solution points for current element type
  CFuint m_nbrSolPnts;

  /// coefs to extrapolate the states to the flx pnts

  std::vector< std::vector< CFreal > > m_solPolyValsAtNodes;

//  Common::SafePtr< std::vector< std::vector< CFreal > > >

  /// nb of nodes in the element
  CFuint m_nbrNodesElem;

  /// order of FR method
  CFuint m_order;

  //Simple edge graph 
  FluxReconstructionMethod::SimpleEdgeGraphFRQ2 m_edgeGraph; 

  std::vector< CFuint > m_nbOfNeighborCellsToaNode;

  std::vector< CFuint > nbOfConnectedFaces;

  std::vector< RealVector > m_vecNodeCoords;

  RealVector m_NodeCoords;
}; // end of class MeshFittingAlgorithm
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FR


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FR_MeshFittingAlgorithmFRQ2_hh
