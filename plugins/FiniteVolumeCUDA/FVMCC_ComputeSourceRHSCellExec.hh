#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRHSCellExec_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRHSCellExec_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeCUDA/FVMCC_ComputeSourceRHSCell.hh"
#include "FiniteVolume/KernelData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class CellConn;
  }
  
  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes with CUDA bindings
 *
 * @author Andrea Lani
 *
 */

//This template include one more parameter, the source term
template <typename SCHEME, typename PHYSICS, typename SOURCE, typename POLYREC, typename LIMITER, CFuint NB_BLOCK_THREADS>
class FVMCC_ComputeSourceRHSCellExec: public FVMCC_ComputeSourceRHSCell<SCHEME, PHYSICS, SOURCE, POLYREC, LIMITER, NB_BLOCK_THREADS> {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeSourceRHSCellExec(const std::string& name) :
    FVMCC_ComputeSourceRHSCell<SCHEME, PHYSICS, SOURCE, POLYREC, LIMITER, NB_BLOCK_THREADS>(name) {}
  
  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeSourceRHSCellExec() {}
  
  /**
   * Execute Processing actions
   */
  virtual void execute();
    
}; // class FVMCC_ComputeSourceRHSCellExec

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRHSCellExec_hh
