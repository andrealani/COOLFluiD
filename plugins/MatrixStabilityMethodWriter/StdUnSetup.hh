#ifndef COOLFluiD_Numerics_StdUnSetup_hh
#define COOLFluiD_Numerics_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "MatrixStabilityMethodWriter/MatrixStabilityMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a unsets the MatrixStabilityMethodWriter method.
 *
 * @author Kris Van den Abeele
 */
class StdUnSetup : public MatrixStabilityMethodCom {
public:

  /**
   * Constructor.
   */
  explicit StdUnSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdUnSetup()
  {
  }

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected:

  /// socket for rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_StdUnSetup_hh
