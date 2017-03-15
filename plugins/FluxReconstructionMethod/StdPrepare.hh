#ifndef COOLFluiD_FluxReconstructionMethod_StdPrepare_hh
#define COOLFluiD_FluxReconstructionMethod_StdPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to prepare
 * the computation by clearing the rhs and the update coefficients
 * @author Tiago Quintino
 * @author Kris Van den Abeele
 */
class StdPrepare : public FluxReconstructionSolverCom {

public: // functions

  /**
   * Constructor.
   */
  explicit StdPrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~StdPrepare();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // data

  /// socket for Rhs
  Framework::DataSocketSink<
                              CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<
                              CFreal> socket_updateCoeff;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

}; // class Prepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdPrepare_hh
