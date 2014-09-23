#ifndef COOLFluiD_Numerics_HessianEE_StdUpdateCom_hh
#define COOLFluiD_Numerics_HessianEE_StdUpdateCom_hh

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
class StdUpdateCom : public HessEECom
{
public:

  /**
   * Constructor.
   */
  explicit StdUpdateCom(const std::string& name);

  /**
   * Destructor.
   */
  ~StdUpdateCom();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // helper methods
  void  writeTEC( const std::string& fname);

private: // data

  // the socket to the data handle of the node's
  Framework::DataSocketSink<RealMatrix> socket_glob_metric;

  // the socket to the data handle of the node's
  Framework::DataSocketSink<RealMatrix> socket_metric;

  // the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// reference lenght in mesh
  CFreal _refH;

}; // class StdUpdateCom

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_StdUpdateCom_hh

