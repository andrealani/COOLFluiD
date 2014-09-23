#ifndef COOLFluiD_Numerics_HessianEE_FunctionCalc_hh
#define COOLFluiD_Numerics_HessianEE_FunctionCalc_hh

//////////////////////////////////////////////////////////////////////////////

#include "HessEEData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"

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
class FunctionCalc : public HessEECom
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
  explicit FunctionCalc(const std::string& name);

  /**
   * Destructor.
   */
  ~FunctionCalc();

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

private: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// @todo missing documentation
  Framework::DataSocketSink<CFreal> socket_adapt_func;

  /// the socket to the data handle of the nStatesProxy
  Framework::DataSocketSink<Framework::ProxyDofIterator<RealVector>* > socket_nstatesProxy;

  /// reference lenght in mesh
  CFreal _refH;

}; // class FunctionCalc

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_FunctionCalc_hh

