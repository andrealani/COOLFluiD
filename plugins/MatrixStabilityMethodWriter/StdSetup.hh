#ifndef COOLFluiD_Numerics_StdSetup_hh
#define COOLFluiD_Numerics_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "MatrixStabilityMethodData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/State.hh"
#include "Framework/Storage.hh"
#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a sets up the MatrixStabilityMethodWriter method.
 *
 * @author Kris Van den Abeele
 */
class StdSetup : public MatrixStabilityMethodCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdSetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

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

  /// socket for rhs
  Framework::DataSocketSource<CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSource<CFreal> socket_updateCoeff;

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_StdSetup_hh

