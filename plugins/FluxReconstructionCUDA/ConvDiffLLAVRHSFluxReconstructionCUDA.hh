#ifndef COOLFluiD_FluxReconstructionCUDA_ConvDiffLLAVRHSFluxReconstructionCUDA_hh
#define COOLFluiD_FluxReconstructionCUDA_ConvDiffLLAVRHSFluxReconstructionCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/KernelData.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#include "Framework/MathTypes.hh"
#include "Framework/VarSetTransformerT.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class CellConn;
  }

    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the convective-diffusive RHS using
 * FR with CUDA bindings
 *
 * @author Ray Vandenhoeck
 *
 */
template <typename SCHEME, typename PHYSICS, typename PHYSICSNS, CFuint ORDER, CFuint NB_BLOCK_THREADS>
class ConvDiffLLAVRHSFluxReconstructionCUDA : public ConvRHSFluxReconstruction {
public:

  /**
   * Constructor
   */
  explicit ConvDiffLLAVRHSFluxReconstructionCUDA(const std::string& name);
  
  /**
   * Destructor
   */
  virtual ~ConvDiffLLAVRHSFluxReconstructionCUDA();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Un Setup private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  virtual void execute();
      
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:
  
  /// Initialize the computation of RHS
  virtual void initializeComputationRHS();
  
  /// Store the stencil data
  virtual void storeStencilData();
  
  /// Copy the local connectivity data to GPU
  void copyLocalCellConnectivity();
  
protected:
  
  /// storage for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage for the gradients
  Framework::DataSocketSink< CFreal >  socket_gradientsCUDA;
  
  /// storage for the gradients
  Framework::DataSocketSink< CFreal >  socket_gradientsAVCUDA;
  
  /// storage for solution point geometric jacobians
  Framework::DataSocketSink< CFreal >  socket_volumes;
  
  /// storage for the cell volumes
  Framework::DataSocketSink< CFreal >  socket_cellVolumes;
  
  /// number of cell flux points
  CFuint m_nbrFlxPnts;
  
  /// storage for the solution point normals
  Framework::DataSocketSink< CFreal > socket_solPntNormals;
  
  /// storage for the flux point normals
  Framework::DataSocketSink< CFreal > socket_flxPntNormals;

  /// storage for the face directions
  Framework::DataSocketSink< CFint > socket_faceDir;
  
  /// cell-face connectivity
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_cellFaces;
  
  /// cell-nodes connectivity
  Common::SafePtr< Common::ConnectivityTable<CFuint> > m_cellNodes;
  
  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder2;
  
  /// pointer to booleans telling whether a face is on the boundary
  Common::SafePtr< std::vector< bool > > m_isFaceOnBoundaryCell;
  
  // cell connectivity
  Framework::LocalArray<Framework::CellConn>::MALLOC_TYPE m_cellConn;
  
  // cell-state connectivity
  Framework::LocalArray<Framework::CellConn>::MALLOC_TYPE m_cellStateConn;
  
  /// number of blocks in x 
  CFuint m_nbBlocksPerGridX;
  
  /// number of blocks in y
  CFuint m_nbBlocksPerGridY;
  
  /// number of cells per block
  CFuint m_nbCellsPerBlock;

  /// number of OpenMP threads
  CFuint m_nbThreadsOMP;

  /// controlling parameter kappa
  CFreal m_kappa;
  
  /// peclet number
  CFreal m_peclet;
  
  /// number of neighbors for each node
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_nbNodeNeighbors;
  
  /// number of corner nodes for current element type
  CFuint m_nbrCornerNodes;
  
  /// vector containing pointers to the nodes in a face
  std::vector< Framework::Node*  >* m_faceNodes;
  
  /// flag telling whether to compute the number of node neighbors
  bool m_flagComputeNbNghb;
  
  /// polynomial coefficients for reconstruction of the artificial viscosity at the flx pnts
  std::vector< std::vector< CFreal > > m_nodePolyValsAtFlxPnts;
  
  /// polynomial coefficients for reconstruction of the artificial viscosity at the sol pnts
  std::vector< std::vector< CFreal > > m_nodePolyValsAtSolPnts;
  
//  /// cell node connectivity table
//  Common::SafePtr< Framework::MeshData::ConnTable > m_cellNodesConn;
  
   /// iteration after which the limiter is frozen
  CFuint m_freezeLimiterIter;
  
  /// boolean telling whether to use max artificial viscosity wrt previous iteration
  bool m_useMax;
  
  /// transformation matrices to order P-1
  RealMatrix m_transformationMatrix;
  
  /// index of the monitored variable for LLAV
  CFuint m_monitoredVar;
  
  /// subcell resolution
  CFreal m_subcellRes;
  
  /// index of the monitored physical variable for LLAV
  CFuint m_monitoredPhysVar;
  
  /// reference smoothness
  CFreal m_s0;
  
  /// bool whether to impose zero LLAV BC
  bool m_LLAVBCZero;
  
  /// IDs of the neighbor node IDs per cellID
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_neighbNodeIDs;
  
  /// node eps
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_nodeEpsilons;
  
  /// cell eps
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_cellEpsilons;
  
  /// node IDs of face nodes
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_faceNeighbNodeIDs;
  
  /// number of face nodes
  CFuint m_nbFaceNodes;
  
  /// IDs of the solution points per cellID
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_stateIDs;
  
  /// IDs of the neighbor cell IDs per cellID
  Framework::LocalArray<CFint>::MALLOC_TYPE m_neighbCellIDs;

  /// bools telling whether the inner cell is LEFT or RIGHT 
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_innerCellIsLeft;
  
  /// IDs of the neighbor face IDs per cellID
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_neighbFaceIDs;
  
  /// sol sol dep in different format
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_solSolDep2;
  
  /// sol flx dep in different format
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_solFlxDep2;
  
  /// flx sol dep in different format
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_flxSolDep2;
  
  /// derivatives of base polynomials in different format
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_solPolyDerivAtSolPnts2;
  
  /// values of base polynomials in different format
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_solPolyValsAtFlxPnts2;
  
  /// flx pnt normal directions in different format
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_flxPntFlxDim2;
  
  /// divergence of correction functions in different format
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_corrFctDiv2;
  
  /// connectivity between faces and flx pnt indexes in different format
  Framework::LocalArray<CFuint>::MALLOC_TYPE m_faceFlxPntConn2;

  /// face integration coefficients
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_faceIntegrationCoefs2;
  
  /// transformation matrices to order P-1
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_transformationMatrix2;
  
  /// polynomial coefficients for reconstruction of the artificial viscosity at the flx pnts
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_nodePolyValsAtFlxPnts2;
  
  /// polynomial coefficients for reconstruction of the artificial viscosity at the sol pnts
  Framework::LocalArray<CFreal>::MALLOC_TYPE m_nodePolyValsAtSolPnts2;
   
  /// flag telling to solve on GPU
  bool m_onGPU;
  
  /// ratio between conv and diff CFL
  CFreal m_cflConvDiffRatio;
  
  /// boolean telling whether to add the contribution of the artificial flux to the update coefficients
  bool m_addUpdCoeff;
  
}; // class FVMCC_ComputeRHSCell

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ConvDiffLLAVRHSFluxReconstructionCUDA.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionCUDA_ConvDiffLLAVRHSFluxReconstructionCUDA_hh
