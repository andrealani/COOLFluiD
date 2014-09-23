#ifndef COOLFluiD_Numerics_LUSGSMethod_StdPrepare_hh
#define COOLFluiD_Numerics_LUSGSMethod_StdPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "LUSGSIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class is a prepare step and it is used to backup the old solution.
   * @author Matteo Parsani 
   * @author Kris Van den Abeele
   */
class StdPrepare : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit StdPrepare(std::string name);

  /**
   * Destructor.
   */
  ~StdPrepare() {}

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected:

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to past states
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// handle to the IDs of the states in each set of states
  Framework::DataSocketSink< std::vector< CFuint > > socket_statesSetStateIDs;

  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// pointer to the auxiliary rhs variable
  Common::SafePtr< RealVector > m_resAux;


}; // class StdPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_StdPrepare_hh

