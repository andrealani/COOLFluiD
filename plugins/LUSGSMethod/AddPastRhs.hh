#ifndef COOLFluiD_Numerics_LUSGSMethod_AddPastRhs_hh
#define COOLFluiD_Numerics_LUSGSMethod_AddPastRhs_hh

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
   * This class stores the rhs of the past time step.
   * @author Kris Van den Abeele
   */
class AddPastRhs : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit AddPastRhs(std::string name);

  /**
   * Destructor.
   */
  ~AddPastRhs() {}

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks.
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// handle to the IDs of the states in each set of states
  Framework::DataSocketSink< std::vector< CFuint > > socket_statesSetStateIDs;

  /// handle to list of booleans telling whether a states set is parallel updatable
  Framework::DataSocketSink< bool > socket_isStatesSetParUpdatable;

  /// socket for rhs of past time step
  Framework::DataSocketSink< CFreal > socket_pastRhs;

  /// number of equations
  CFuint m_nbrEqs;

}; // class AddPastRhs

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_AddPastRhs_hh

