#ifndef COOLFluiD_Numerics_SubSystemCoupler_FVMCCReadDataTransfer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FVMCCReadDataTransfer_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdReadDataTransfer.hh"
#include "Framework/FaceTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a FVMCCReadDataTransfer Coupler command
 *
 * @author Thomas Wuilbaut
 *
 */
class FVMCCReadDataTransfer : public StdReadDataTransfer {
public:

  /**
   * Constructor.
   */
  explicit FVMCCReadDataTransfer(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCCReadDataTransfer();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /**
   * Use the PostVariableTransformer to transform
   * the values received through the datahandle
   */
  virtual void transformReceivedGhostData();

  /**
   * Use the PostVariableTransformer to transform
   * the values received through the datahandle
   */
  virtual void transformReceivedNodalData();

private:

  /// socket for nodal States
  Framework::DataSocketSink< RealVector> socket_nstates;

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;


}; // class FVMCCReadDataTransfer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FVMCCReadDataTransfer_hh

