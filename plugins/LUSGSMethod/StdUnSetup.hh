#ifndef COOLFluiD_Numerics_LUSGSMethod_StdUnSetup_hh
#define COOLFluiD_Numerics_LUSGSMethod_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////
#include <vector>

#include "Framework/DataSocketSink.hh"
#include "LUSGSMethod/LUSGSIteratorData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * LUSGS method.
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class StdUnSetup : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit StdUnSetup(std::string name);

  /**
   * Destructor.
   */
  ~StdUnSetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected:

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink< CFreal > socket_updateCoeff;

  /// socket for the past states
  Framework::DataSocketSink< Framework::State* > socket_pastStates;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// socket for the pivot element of the LU factorization
  Framework::DataSocketSink< std::vector< CFuint > > socket_pivotLUFactorization;


}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_StdUnSetup_hh
