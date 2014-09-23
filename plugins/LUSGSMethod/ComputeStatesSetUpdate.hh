#ifndef COOLFluiD_Numerics_LUSGSMethod_ComputeStatesSetUpdate_hh
#define COOLFluiD_Numerics_LUSGSMethod_ComputeStatesSetUpdate_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "LUSGSMethod/LUSGSIteratorData.hh"
#include "Framework/DataSocketSink.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class updates the index of the states set which will be updated.
   * @author Matteo Parsani
   * @author Kris Van den Abeele
   */
class ComputeStatesSetUpdate : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit ComputeStatesSetUpdate(std::string name);

  /**
   * Destructor.
   */
  ~ComputeStatesSetUpdate() {}

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

  void solveTriangularSystems(const RealMatrix& lhsMatrix, RealVector& rhs);

protected: // data

  /// socket for diagonal block Jacobian matrices
  Framework::DataSocketSink < RealMatrix > socket_diagBlockJacobMatr;

  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// handle to list of booleans telling whether a states set is parallel updatable
  Framework::DataSocketSink< bool > socket_isStatesSetParUpdatable;

  /// pointer to the auxiliary rhs variable
  Common::SafePtr< RealVector > m_resAux;

}; // class ComputeStatesSetUpdate

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_ComputeStatesSetUpdate_hh

