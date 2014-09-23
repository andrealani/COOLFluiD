#ifndef COOLFluiD_Numerics_FluctSplit_ComputeRhsJacobCoupling_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeRhsJacobCoupling_hh

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
class FluctSplit_API ComputeRhsJacobCoupling : public ComputeRHS {
public:

  /// Constructor.
  explicit ComputeRhsJacobCoupling(const std::string& name);

  /// Destructor.
  virtual ~ComputeRhsJacobCoupling();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

protected: // functions

  /// Execute the command on the current TRS
  virtual void executeOnTrs();
   
protected:

  /// acquaintance of all the linear system solvers
  std::vector<Common::SafePtr<Framework::LinearSystemSolver> > _lss;
   
  // accumulator for all the LSSMatrix's
  std::vector<std::vector<Framework::BlockAccumulator*> > _acc;
  
  // array of the equation IDs for each LSS
  std::vector<Common::SafePtr<std::vector<CFuint> > > _equations;
  
  /// strategy object to compute the jacobian contributions
  Common::SafePtr<ComputeJacobStrategy> _jacobStrategy;
 
}; // class ComputeRhsJacobCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeRhsJacobCoupling_hh
