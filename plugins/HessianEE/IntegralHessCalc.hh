#ifndef COOLFluiD_Numerics_HessianEE_IntegralHessCalc_hh
#define COOLFluiD_Numerics_HessianEE_IntegralHessCalc_hh

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
class IntegralHessCalc : public HessEECom
{
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit IntegralHessCalc(const std::string& name);

  /**
   * Destructor.
   */
  ~IntegralHessCalc();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Execute Processing actions
   */
  void execute();
  //void executeOnTrs();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


private: // helper methods
  void  writeTEC( const std::string& fname);
  void  calcGrad();
  void  calcHess();

private: // data

  /// storage for the gradient of solution
  Framework::DataSocketSink<RealVector> socket_grad;

  /// @todo missing documentation
  Framework::DataSocketSink<CFreal> socket_adapt_func;

  /// storage for the hessian at each state
  Framework::DataSocketSink<RealMatrix> socket_adapt_hess;

  /// storage for auxiliary weight
  Framework::DataSocketSink<CFreal> socket_adapt_wght;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// reference lenght in mesh
  CFreal _refH;

}; // class IntegralHessCalc

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_IntegralHessCalc_hh

