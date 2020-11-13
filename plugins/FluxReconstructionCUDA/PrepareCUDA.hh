#ifndef COOLFluiD_FluxReconstructionMethod_PrepareCUDA_hh
#define COOLFluiD_FluxReconstructionMethod_PrepareCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/StdPrepare.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to prepare
 * the computation on GPU by clearing the rhs and the update coefficients
 * @author Ray Vandenhoeck
 */
class PrepareCUDA : public StdPrepare {

public: // functions

  /**
   * Constructor.
   */
  explicit PrepareCUDA(const std::string& name);

  /**
   * Destructor.
   */
  ~PrepareCUDA();

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

  /// socket for gradients
  Framework::DataSocketSink< CFreal > socket_gradientsCUDA;
  
  /// socket for AV gradients
  Framework::DataSocketSink< CFreal > socket_gradientsAVCUDA;
  

}; // class Prepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PrepareCUDA_hh
