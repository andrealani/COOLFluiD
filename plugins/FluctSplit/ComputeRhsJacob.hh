#ifndef COOLFluiD_Numerics_FluctSplit_ComputeRhsJacob_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeRhsJacob_hh

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

/// A command that computes the RHS and the Jacobian fro RDS
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API ComputeRhsJacob : public ComputeRHS {
public:

  /// Constructor.
  explicit ComputeRhsJacob(const std::string& name);

  /// Destructor.
  virtual ~ComputeRhsJacob();

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
  
}; // class ComputeRhsJacob

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeRhsJacob_hh
