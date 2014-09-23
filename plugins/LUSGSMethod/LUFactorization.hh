#ifndef COOLFluiD_Numerics_LUSGSMethod_LUFactorization_hh
#define COOLFluiD_Numerics_LUSGSMethod_LUFactorization_hh

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
   * This class computes the LU factorization of the block diagonal Jacobian matrices.
   * @author Matteo Parsani
   * @author Kris Van den Abeele
   */
class LUFactorization : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit LUFactorization(std::string name);

  /**
   * Destructor.
   */
  ~LUFactorization() {}

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

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// number of equations
  CFuint m_nbrEqs;

}; // class LUFactorization

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_LUFactorization_hh

