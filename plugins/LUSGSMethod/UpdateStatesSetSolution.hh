#ifndef COOLFluiD_Numerics_LUSGSMethod_UpdateStatesSetSolution_hh
#define COOLFluiD_Numerics_LUSGSMethod_UpdateStatesSetSolution_hh

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
class UpdateStatesSetSolution : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit UpdateStatesSetSolution(std::string name);

  /**
   * Destructor.
   */
  ~UpdateStatesSetSolution() {}

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

protected: // data

  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// socket for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to the IDs of the states in each set of states
  Framework::DataSocketSink< std::vector< CFuint > > socket_statesSetStateIDs;

  /// handle to list of booleans telling whether a states set is parallel updatable
  Framework::DataSocketSink< bool > socket_isStatesSetParUpdatable;

  /// number of equations
  CFuint m_nbrEqs;

}; // class UpdateStatesSetSolution

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_UpdateStatesSetSolution_hh

