#ifndef COOLFluiD_Numerics_FluctSplit_DiffusionRhsJacob_hh
#define COOLFluiD_Numerics_FluctSplit_DiffusionRhsJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeRHS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework
  {
    class State;
    class BlockAccumulator;
  }

    namespace FluctSplit {

      class ComputeJacobStrategy;

//////////////////////////////////////////////////////////////////////////////

/// A command that computes the RHS and the Jacobian for discretizing
/// diffusion equations
/// @author Andrea Lani
class FluctSplit_API DiffusionRhsJacob : public ComputeRHS {
public:

  /// Constructor.
  explicit DiffusionRhsJacob(const std::string& name);

  /// Destructor.
  virtual ~DiffusionRhsJacob();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

protected: // functions

  /// Execute the command on the current TRS
  virtual void executeOnTrs();

protected:

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> _lss;
  
  /// vector of LSSMatrix accumulators (one for each cell type)
  std::vector<Framework::BlockAccumulator*> _acc;
  
  /// strategy object to compute the jacobian contributions
  Common::SafePtr<ComputeJacobStrategy> _jacobStrategy;
  
}; // class DiffusionRhsJacob

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_DiffusionRhsJacob_hh
