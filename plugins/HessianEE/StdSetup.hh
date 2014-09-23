#ifndef COOLFluiD_Numerics_HessianEE_StdSetup_hh
#define COOLFluiD_Numerics_HessianEE_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "HessEEData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * executed in order to setup the HessianEE
 *
 * @author Jurek Majevski
 */
class StdSetup : public HessEECom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(std::string name);

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
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected:

  /// @todo missing documentation
  Framework::DataSocketSource<CFreal> socket_adapt_func;

  /// storage for the gradient of solution
  Framework::DataSocketSource<RealVector> socket_grad;

  /// storage for the hessian at each state
  Framework::DataSocketSource<RealMatrix> socket_hessian;

  /// storage for the global metric at each state
  Framework::DataSocketSource<RealMatrix> socket_glob_metric;

  /// storage for the metric for remeshing
  Framework::DataSocketSource<RealMatrix> socket_metric;

  /// storage for auxiliary weight
  Framework::DataSocketSource<CFreal> socket_adapt_wght;

  /// storage for the auxiliary matrix at each state
  Framework::DataSocketSource<RealMatrix> socket_adapt_matrix;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_StdSetup_hh

