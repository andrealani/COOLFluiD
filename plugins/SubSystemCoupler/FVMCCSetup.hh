#ifndef COOLFluiD_Numerics_SubSystemCoupler_FVMCCSetup_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FVMCCSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "SubSystemCoupler/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to setup the SubSystemCoupler Method
 *
 * @author Thomas Wuilbaut
 *
 */
class FVMCCSetup : public StdSetup {
public:

  /**
   * Constructor.
   */
  explicit FVMCCSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCCSetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected: // data

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// socket for Ghost State's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

}; // class FVMCCSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FVMCCSetup_hh

