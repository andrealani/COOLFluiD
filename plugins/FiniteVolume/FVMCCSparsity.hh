#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCCSparsity_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCCSparsity_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GlobalJacobianSparsity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * Computes the sparticity of the global jacobian matrix for a cell-centered
 * space discretization methods.
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 */
class FVMCCSparsity : public Framework::GlobalJacobianSparsity {

public: // functions

  /**
   * Constructor.
   */
  FVMCCSparsity();

  /**
   * Destructor.
   */
  virtual ~FVMCCSparsity();

  /**
   * Computes the non zero entries in the global jacobian matrix
   */
  virtual void computeNNz(std::valarray<CFint>& nnz,
                          std::valarray<CFint>& ghostNnz);

  /// Computes the non zero entries in the global jacobian matrix
  /// using nodes instead of states
  virtual void computeNNzNodeBased(std::valarray<CFint>& nnz,
				   std::valarray<CFint>& ghostNnz);
  
  /**
   * Computes the non zero entries in the global jacobian matrix.
   * Also computes the vertex-vertex connectivity.
   */
  virtual void computeMatrixPattern(
    std::valarray<CFint>& nnz,
    std::valarray<CFint>& ghostNnz,
    std::vector<std::vector<CFuint> >& matrixPattern);
 
  /// Computes the matrix patern and stores it in a connectivity table
  virtual void computeMatrixPattern
  (Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
   Common::ConnectivityTable<CFuint>& matrixPattern);
  
}; // class FVMCCSparsity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCCSparsity_hh
