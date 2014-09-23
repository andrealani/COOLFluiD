#ifndef COOLFluiD_Numerics_MatrixStabilityMethodWriter_SetStates_hh
#define COOLFluiD_Numerics_MatrixStabilityMethodWriter_SetStates_hh

//////////////////////////////////////////////////////////////////////////////

#include "MatrixStabilityMethodWriter/MatrixStabilityMethodData.hh"
#include "MathTools/RealVector.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This commands sets the states to the appropriate values (all zero but one)
   * @author Kris Van den Abeele
   */
class SetStates : public MatrixStabilityMethodCom {
public:

  /**
   * Constructor.
   */
  explicit SetStates(const std::string& name);

  /**
   * Destructor.
   */
  ~SetStates()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class SetStates

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MatrixStabilityMethodWriter_SetStates_hh
