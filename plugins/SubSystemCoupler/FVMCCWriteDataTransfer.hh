#ifndef COOLFluiD_Numerics_SubSystemCoupler_FVMCCWriteDataTransfer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FVMCCWriteDataTransfer_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdWriteDataTransfer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a FVMCCWriteDataTransfer Coupler command
 *
 * @author Thomas Wuilbaut
 *
 */
class FVMCCWriteDataTransfer : public StdWriteDataTransfer {
public:

  /**
   * Constructor.
   */
  explicit FVMCCWriteDataTransfer(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCCWriteDataTransfer();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute Processing actions
   */
  void executeWrite(const CFuint iProc);

private: // data

  // the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;

}; // class FVMCCWriteDataTransfer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FVMCCWriteDataTransfer_hh

