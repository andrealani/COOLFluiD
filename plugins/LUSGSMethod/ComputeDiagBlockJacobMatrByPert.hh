#ifndef COOLFluiD_Numerics_LUSGSMethod_ComputeDiagBlockJacobMatrByPert_hh
#define COOLFluiD_Numerics_LUSGSMethod_ComputeDiagBlockJacobMatrByPert_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "Framework/DataSocketSink.hh"
#include "LUSGSMethod/LUSGSIteratorData.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class computes the block diagonal Jacobian matrices by perturbation of the states.
   * @author Matteo Parsani
   * @author Kris Van den Abeele
   */
class ComputeDiagBlockJacobMatrByPert : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit ComputeDiagBlockJacobMatrByPert(std::string name);

  /**
   * Destructor.
   */
  ~ComputeDiagBlockJacobMatrByPert() {}

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks.
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: // functions

  void factorizeMatrix(const CFuint diag, RealMatrix& matrix);

protected: //Data

  /// socket for diagonal block Jacobian matrices
  Framework::DataSocketSink < RealMatrix > socket_diagBlockJacobMatr;

  /// handle to list of booleans telling whether a states set is parallel updatable
  Framework::DataSocketSink< bool > socket_isStatesSetParUpdatable;

  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// socket for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to the IDs of the states in each set of states
  Framework::DataSocketSink< std::vector< CFuint > > socket_statesSetStateIDs;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// pointer to the auxiliary rhs variable
  Common::SafePtr< RealVector > m_unpertRhs;

  /// number of equations
  CFuint m_nbrEqs;

  /// iterator for the states in a state set
  CFuint m_stateIdx;

  /// iterator for the physical variables in a state
  CFuint m_eqIdx;

}; // class ComputeDiagBlockJacobMatrByPert

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_ComputeDiagBlockJacobMatrByPert_hh

