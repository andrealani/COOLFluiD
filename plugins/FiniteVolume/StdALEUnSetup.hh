#ifndef COOLFluiD_Numerics_FiniteVolume_StdALEUnSetup_hh
#define COOLFluiD_Numerics_FiniteVolume_StdALEUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Thomas Wuilbaut
 */
class StdALEUnSetup : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit StdALEUnSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdALEUnSetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  virtual void execute();

protected:

  /// Socket for the past coordinates of the Node's
  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;

  /// Socket for the future coordinates of the Node's
  Framework::DataSocketSink<Framework::Node*> socket_futureNodes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StdALEUnSetup_hh

