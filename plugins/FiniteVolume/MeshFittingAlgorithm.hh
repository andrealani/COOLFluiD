#ifndef COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh
#define COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "SimpleEdgeGraph.hh"
#include "Framework/TRSDistributeData.hh"
#include "MeshTools/ComputeWallDistanceVector2CCMPI.hh"
#include "Common/CFMultiMap.hh"
#include "Framework/GeometricEntityPool.hh"


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


  CFreal computeSecondSpringConstant(const Framework::Node* const firstNode, 
			             const Framework::Node* const secondNode);
  
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
  
  CFreal computetorsionConstant3DSemiOnly(CFuint  nodeID1,CFuint nodeID2 ,
					  const  Framework::Node* const  firstNode,
					  const  Framework::Node* const  secondNode);
  
  /**
   * Compute the volume of a tetrahedral element
   */
  
  CFreal  ComputeTvolume(CFuint  n1, CFuint  n2,
			 CFuint n3, CFuint n4);
  
  /**
   * Compute the face area of a tetrahedral element
   */
  CFreal ComputeTFacesurface(Framework::Node* node1, 
			     Framework::Node* node2,
			     Framework::Node* node3);
  
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

  CFreal computeConstantquads3D(const  Framework::Node* const  firstNode,
				const  Framework::Node* const  secondNode,
				const Framework::Node* const  thirdNode,
				const Framework::Node* const  fourthNode);

  CFreal computeConstantquads3DElemArea(const  RealVector xyi , const  Framework::Node* const  firstNode,
				const  Framework::Node* const  secondNode,
				const Framework::Node* const  thirdNode,
				const Framework::Node* const  fourthNode, 
				CFreal elementArea);
  /** Compute the skewness of a quadrilateral element 
   *
   */
  CFreal computeSkewness2dQuads(const  Framework::Node* const  a1,const  Framework::Node* const  b1,
				const  Framework::Node* const  c1, const  Framework::Node* const  d1);
  
  /**
   * assemble the system lines for a node moving in along the boundary
   */
  void assembleMovingInBoundaryNode(const Framework::Node* node);
  void assembleMovingInBoundaryNode3DHexa(const Framework::Node* node);
  /**
   * Queries if the node is locked 
   */
  bool isNodeLocked( Framework::Node* node);
  /**
   * Queries if the node is inside the region 
   */
  bool insideRegion(Framework::Node* node);
  RealVector  insideRegionExtrapolate(Framework::Node* node);
 
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
  
  void assembleinRegionNode2DQuadsGeoBased(const  Framework::Node* node);

   /**
   * checks if 4 given points are coplanar
   */
   bool areCoplanar(const  Framework::Node* const  firstNode,
		    const  Framework::Node* const  secondNode,
		    const Framework::Node* const  thirdNode,
		    const Framework::Node* const  fourthNode);
  /**
   * assemble the systme lines for 3D hexahedral element
   */
  void assembleinRegionNode3DHexa(const  Framework::Node* node);

  /**
 * assemble system lines for 3D tetrahedral elemnt 
 */
  void assembleinRegionNode3DTet(const  Framework::Node* node);
  
  /**
   * Compute the intersection point between two quadrilateral diagonals
   */
  RealVector computeIntersection(const  Framework::Node* const  a, const  Framework::Node* const  c,
				 const  Framework::Node* const  b, const  Framework::Node* const  d);

  RealVector computeIntersectionQuad3D(const  Framework::Node* const  a,const  Framework::Node* const  c,
						      const  Framework::Node* const  b, const  Framework::Node* const  d);
  
  /**
   * assemble the system lines for an inner Node
   */
  void assembleInnerNode(const Framework::Node* node);
  
  /**
   * Update nodes with new positions
   */
  void updateNodePositions();
  
  void smooth(const  Framework::Node* node);

  void computeWallDistanceExtrapolate();

  /**
   * Send message to trigger re-computation of mesh information
   */
  void triggerRecomputeMeshData();
  

  void stateInterpolator();

  void interpolateNewMeshProperties();

  void saveOldMeshProperties();


//// Multiple springs
 void assembleInnerNodeSecond(const Framework::Node* node);

private: //data
  /// the socket of the nodes connectivity for 2D quadrilatral mesh
  Common::CFMultiMap <CFuint , CFuint>  m_mapNodeNode1;
  
  /// the socket of the 3D terahedral mesh connectivity 
  Common::CFMultiMap <CFuint , CFuint>  m_mapCellNode1;

  /// the socket of the 3D hexahedral mesh connectivity 
  Common::CFMultiMap <CFuint , CFuint>  m_mapNodeCell1;
  
  /// the socket of cell internal radius at setup 
  Common::CFMultiMap <CFuint , CFreal>  m_mapNodeRadius1;
  
  /// the socket of states at setup 
  Common::CFMultiMap <CFuint , CFreal>  m_mapNodeNState1;
  
  /// the socket to the data handle of the nodal stiffness
  Framework::DataSocketSource < CFreal > socket_stiffness;
  
  /// the socket to the data handle of the relative error 
  Framework::DataSocketSource < CFreal > socket_relativeError;
  
  /// the socket to the data handle of the internal radius circle
  Framework::DataSocketSource < CFreal > socket_iradius;
  
  /// the socket to the data handle of the skewness
  Framework::DataSocketSource < CFreal > socket_skewness;

  /// the socket to the data handle of the aspect ratio
  Framework::DataSocketSource < CFreal > socket_AR;

  /// the socket to the data handle of the internal radius of sphere
  Framework::DataSocketSource < CFreal > socket_isphere;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  Framework::DataSocketSink<std::vector< Framework::State*> > socket_stencil;



  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the RHS
  Framework::DataSocketSink<CFreal> socket_rhs;


  /// socket for the wallDistance storage
  Framework::DataSocketSink<CFreal> socket_wallDistance;

  /// socket for the node inside Region 
  // True if node is inside region 
  // False if the node is outside the region 
  Framework::DataSocketSink<bool> socket_nodeisAD; 
 
  Framework::DataSocketSink<CFreal> socket_nodeDistance;

 
  /// data handle cashed for efficiency  reasons
  Framework::DataHandle<CFreal> m_wallDistance;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;
  
  // pointer to the data of the cell centered FVM method
  Common::SafePtr<Numerics::FiniteVolume::CellCenterFVMData> m_fvmccData;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder1;
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceTRSBuilder;

  
  // physical data array
  RealVector m_pdata;
 
  // dummy state vector
  std::auto_ptr<Framework::State> m_state; 
  
  /// Minimum percentile limit for spring constants
  CFreal m_minPercentile;

  /// Maximum percentile limit for spring constants
  CFreal m_maxPercentile;
  
  /// Step for mesh adaptation
  CFreal m_meshAcceleration;

  /// Tolerance on the mesh mouvement for RSI
  CFreal m_tolerance;

  CFreal m_ratioBoundaryToInnerEquilibriumSpringLength;

  std::vector<std::string> m_unlockedBoundaryTRSs;
  
  /// Monitor variable ID from State
  CFuint m_monitorVarID;
  
  /// Monitor variable ID from physical data
  CFuint m_monitorPhysVarID;

  /// MQIvalue to be defined 
  CFuint m_MQIvalue;
  
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

  FiniteVolume::SimpleEdgeGraph m_edgeGraphN;

  // Storage of element skewness 
  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeSkew1;

  // Storage of element aspect ratio  
  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeAR1;

  // Storage of internal radius of the sphere
  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeTS1;

  //RealVector Extrapolated;
  //RealVector m_vectorNodeElementVolume; 
  //RealVector m_vectorNodeElementSkewness;
  //RealVector m_vectorNodeElementSkewnessOut;

  // Hybrid spring analogy
  CFreal m_acceptableDistance;

  // ST based on the middle angle for quad meshes
  bool  m_thetaMid;
  bool m_interpolateState;
  bool m_smoothSpringNetwork;
  bool m_smoothNodalDisp;
  RealVector oldStates;

  RealVector oldCoordinates;

  // Global RSI 
  CFreal m_erreurG;

  CFuint movingIter;

}; // end of class MeshFittingAlgorithm
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MeshFittingAlgorithm_hh
