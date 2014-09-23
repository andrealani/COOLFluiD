#ifndef COOLFluiD_Numerics_HessianEE_StdMetricCalc_hh
#define COOLFluiD_Numerics_HessianEE_StdMetricCalc_hh

//////////////////////////////////////////////////////////////////////////////

#include "HessEEData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to ComputeSpaceResidual the MeshData.
 * @author Jurek Majewski
*/
class StdMetricCalc : public HessEECom
{
public:

  /**
   * Constructor.
   */
  explicit StdMetricCalc(const std::string& name);

  /**
   * Destructor.
   */
  ~StdMetricCalc();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Execute Processing actions
   */
  void execute();

private: // helper methods

  /// storage for auxiliary weight
  Framework::DataSocketSink<RealMatrix> socket_metric;

  /// storage for the hessian at each state
  Framework::DataSocketSink<RealMatrix> socket_hessian;

  /// storage for the gradient of solution
  Framework::DataSocketSink<RealVector> socket_grad;

  /// @todo missing documentation
  Framework::DataSocketSink<CFreal> socket_adapt_func;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  void    calcMetric();

  CFreal  calcConst();

  void    limitMetric( const CFreal& c);

  void  writeTEC( const std::string& fname);

}; // class StdMetricCalc

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_StdMetricCalc_hh

