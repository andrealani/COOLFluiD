#ifndef COOLFluiD_FluxReconstructionMethod_MLPLimiter_hh
#define COOLFluiD_FluxReconstructionMethod_MLPLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that applies an elementwise MLP limiter to the solution
 *
 * @author Ray Vandenhoeck
 *
 */
class MLPLimiter : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit MLPLimiter(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~MLPLimiter();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
      providesSockets();
      
   /**
   * apply the limiter to perturbed states
   */
   void limitPerturbedStates(CFuint elemIdx, CFuint pertSol, CFuint pertVar);
   
   /**
   * restore the limited perturbed states
   */
   void restoreStates(CFuint elemIdx);

protected: // functions

  /**
   * reset the node min and max cell averaged states
   */
  void resetNodeNghbrCellAvgStates();

  /**
   * set the reconstruction data for averaged solution in current element type
   */
  void setAvgReconstructionData(CFuint iElemType);

  /**
   * set all the reconstruction data for current element type
   */
  void setAllReconstructionData(CFuint iElemType);

  /**
   * compute the cell averaged state
   */
  virtual void reconstructCellAveragedState();

  /**
   * compute a cell averaged variable
   */
  virtual void reconstructCellAveragedVariable(const CFuint iEq);

  /**
   * compute derivative in the cell center of a variable
   */
  virtual void computeCellCenterDerivVariable(const CFuint iEq);

  /**
   * set minimum and maximum node neighbouring cell average solutions
   * @pre reconstructCellAveragedState()
   */
  void setMinMaxNodeNghbrCellAvgStates();

  /**
   * set minimum and maximum neighbouring cell average solutions
   */
  void setMinMaxNghbrCellAvgStates();
  
  /**
   * Compute the node states
   */
  void computeNodeStates(std::vector< RealVector > states, std::vector< RealVector >& statesNodes);
  
  /**
   * Compute states projected to a lower order
   */
  void computeProjStates(std::vector< RealVector >& projStates, CFuint order);
  
  /**
   * Apply additional checks to the state
   */
  virtual void applyChecks(CFreal phi) {};
  
  /**
   * Check if the states are physical
   */
  virtual bool checkPhysicality() {return true;};
  
  /**
   * Compute the flx pnt states
   */
  void computeFlxPntStates(std::vector< RealVector > states, std::vector< RealVector >& statesFlxPnt);
  
  /**
   * Apply the previously computed limiter instead of recomputing the limiter
   */
  void applyPrevLimiter();
  
  /**
   * Execute the slope limiter on the P1 states
   */
  void executeSlopeLimiter(const CFuint elemIdx, const bool useMin);
  
  /**
   * Update the limiting output socket for current cell
   */
  void updateLimitingOutput(CFreal limitingType);
  
  /**
   * check for special physics-dependent conditions if limiting is necessary
   */
  virtual bool checkSpecialLimConditions() {return true;};

protected: // data

  /// socket for minimum nodal state's
  Framework::DataSocketSource< RealVector > socket_nodeNghbCellMinAvgStates;

  /// socket for maximum nodal state's
  Framework::DataSocketSource< RealVector > socket_nodeNghbCellMaxAvgStates;
  
  /// socket for output of the MLP limiting information
  /// Values: 0.0 = no limiting, 1.0 = order reduction, 2.0 = slope limiting
  Framework::DataSocketSource< CFreal > socket_outputLimiting;
  
  /// socket for limiter values
  std::vector< RealVector > m_limiterValues;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;
  
  /// vector containing backup states
  std::vector< RealVector > m_cellStatesBackup;

  /// vector containing pointers to the nodes in a cell
  std::vector< Framework::Node*  >* m_cellNodes;

  /// coefficients for the computation of the cell averaged solution
  Common::SafePtr< RealVector > m_cellAvgSolCoefs;

  /// coefficients for the computation of the solution derivative
  Common::SafePtr< std::vector< RealVector > > m_cellCenterDerivCoefs;

  /// reconstruction coefficients for the flux points
  Common::SafePtr< RealMatrix > m_flxPntsRecCoefs;

  /// indexes of all the flux points
  Common::SafePtr< std::vector< CFuint > > m_allFlxPntIdxs;

  /// mapped coordinates of the solution points
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// cell averaged state
  RealVector m_cellAvgState;

  /// derivative in cell center of a physical variable
  RealVector m_cellCenterDerivVar;

  /// minimum neighbouring cell averaged state
  RealVector m_minAvgState;

  /// maximum neighbouring cell averaged state
  RealVector m_maxAvgState;

  /// minimum cell averaged state on the mesh
  RealVector m_minAvgStateAll;

  /// maximum cell averaged state on the mesh
  RealVector m_maxAvgStateAll;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// dimensionality of the physical model
  CFuint m_dim;

  /// number of flux points for current element type
  CFuint m_nbrFlxPnts;

  /// number of solution points for current element type
  CFuint m_nbrSolPnts;

  /// number of corner nodes for current element type
  CFuint m_nbrCornerNodes;

  /// boolean telling whether to apply the limiter
  std::vector< bool > m_applyLimiter;

  /// boolean telling whether to apply the limiter to a physical variable
  std::vector< bool > m_applyLimiterToPhysVar;

  /// factor used to determine whether to limit the solution
  CFreal m_tvbLimitFactor;
  
  /// residual after which the limiter is frozen
  CFreal m_freezeLimiterRes;
  
  /// iteration after which the limiter is frozen
  CFuint m_freezeLimiterIter;

  /// exponent for the computation of the length scale
  CFreal m_lengthScaleExp;
  
  /// extrapolated states in the nodes of the cell
  std::vector< RealVector > m_cellStatesNodes;
  
  /// extrapolated P1 states in the nodes of the cell
  std::vector< RealVector > m_cellStatesNodesP1;
  
  /// maximum number of flux points in a cell
  CFuint m_maxNbrFlxPnts;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtNodes;
  
  /// vector of transformation matrices to lower orders
  std::vector< RealMatrix > m_transformationMatrices;
  
  /// states projected on P1
  std::vector< RealVector > m_statesP1;
  
  /// states projected on P1
  std::vector< RealVector > m_states2;
  
  /// nb of nodes in the element
  CFuint m_nbrNodesElem;
  
  /// coefs to compute the derivative of the states in the nodes
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtNodes;
  
  /// order of FR method
  CFuint m_order;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  /// extrapolated states in the flux points of the cell
  std::vector< RealVector > m_cellStatesFlxPnt;
  
  /// limited states of the previous iteration
  std::vector< RealVector > m_prevStates;

  /// showrate of limiter info
  CFuint m_showrate;

}; // class MLPLimiter

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_MLPLimiter_hh
