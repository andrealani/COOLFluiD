#ifndef COOLFluiD_Numerics_FiniteVolume_BDF2ALEUnSetup_hh
#define COOLFluiD_Numerics_FiniteVolume_BDF2ALEUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "CellCenterFVMData.hh"
#include "StdALEUnSetup.hh"

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
class BDF2ALEUnSetup : public StdALEUnSetup {
public:

  /**
   * Constructor.
   */
  explicit BDF2ALEUnSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~BDF2ALEUnSetup();

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

  /// Socket for the past coordinates of the Node's (at time n-1)
  Framework::DataSocketSink<Framework::Node*> socket_pastPastNodes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_BDF2ALEUnSetup_hh

