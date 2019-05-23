#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCellExec_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCellExec_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeCUDA/FVMCC_ComputeSourceRhsJacobCell.hh"

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
class FVMCC_ComputeSourceRhsJacobCellExec : public FVMCC_ComputeSourceRhsJacobCell<SCHEME, PHYSICS, SOURCE, POLYREC, LIMITER, NB_BLOCK_THREADS> {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeSourceRhsJacobCellExec(const std::string& name) :
    FVMCC_ComputeSourceRhsJacobCell<SCHEME, PHYSICS, SOURCE, POLYREC, LIMITER, NB_BLOCK_THREADS>(name) {}
  
  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeSourceRhsJacobCellExec() {}
  
  /**
   * Execute Processing actions
   */
  virtual void execute();
  
}; // class FVMCC_ComputeSourceRhsJacobCellExec

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCellExec_hh
