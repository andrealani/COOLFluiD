#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCell_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCell_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeCUDA/FVMCC_ComputeSourceRHSCell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS with a source term and system jacobian using
 * standard cell center FVM schemes with CUDA bindings
 *
 * @author Andrea Lani
 * @author Isaac Alonso
 *
 */
template <typename SCHEME, typename PHYSICS, typename SOURCE, typename POLYREC, typename LIMITER, CFuint NB_BLOCK_THREADS>
class FVMCC_ComputeSourceRhsJacobCell : public FVMCC_ComputeSourceRHSCell<SCHEME, PHYSICS, SOURCE, POLYREC, LIMITER, NB_BLOCK_THREADS> {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeSourceRhsJacobCell(const std::string& name);
  
  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeSourceRhsJacobCell();
  
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
   * Execute Processing actions
   */
  virtual void execute();
      
protected:
  
  /// Initialize the computation of RHS
  virtual void initializeComputationRHS();
  
  /// Update the system matrix from given kernel 
  virtual void updateSystemMatrix(CFuint kernelID);
  
  /// Finalize the computation of RHS
  virtual void finalizeComputationRHS();
  
  /// Compute the jacobian contribution of the current (boundary) face
  virtual void computeBoundaryJacobianTerm();
  
  /// execute on BCs
  virtual void executeBC();
  
protected:
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;
  
  /// storage of the block jacobian matrices
  CudaEnv::CudaVector<CFuint, CudaEnv::MallocHostAlloc<CFuint> > m_blockStart;
  
  /// number of cells per set (=block of cells to compute on device at once)  
  std::vector<CFuint> m_nbCellsInKernel;
  
  /// arrat storing the start IDs for the matrix contributions coming from each kernel
  std::vector<CFuint> m_blockStartKernel;
  
  /// arrat storing the start cell IDs for matrix contributions coming from each kernel
  std::vector<CFuint> m_blockStartKernelCellID;
  
  /// storage of the block jacobian matrices
  CudaEnv::CudaVector<CFreal> m_blockJacobians;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::NumericalJacobian> _numericalJacob;
  
  /// vector for the temporary flux corresponding to a perturbed state
  RealVector _pertFlux;
  
  /// vector for the temporary flux finite difference
  RealVector _fluxDiff;

  /// vector for the original state values
  RealVector _origState;
  
  /// accumulator for LSSMatrix for boundary contribution
  std::auto_ptr<Framework::BlockAccumulator> _bAcc;
  
  /// number of kernel blocks
  CFuint m_nbKernelBlocks;
  
}; // class FVMCC_ComputeSourceRhsJacobCell

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeSourceRhsJacobCell.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCell_hh
