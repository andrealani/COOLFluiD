#ifndef COOLFluiD_Numerics_HessianEE_StdSmoothCom_hh
#define COOLFluiD_Numerics_HessianEE_StdSmoothCom_hh

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
class StdSmoothCom : public HessEECom
{
public:

  /**
   * Constructor.
   */
  explicit StdSmoothCom(const std::string& name);

  /**
   * Destructor.
   */
  ~StdSmoothCom();

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

  /// storage for the matrix
  Framework::DataSocketSink<RealMatrix> socket_adapt_matrix;

  /// storage for auxiliary weight
  Framework::DataSocketSink<CFreal> socket_adapt_wght;

  /// storage for the matrix
  Framework::DataSocketSink<RealMatrix> socket_hessian;

  /// storage for the matrix
  Framework::DataSocketSink<RealMatrix> socket_metric;

}; // class StdSmoothCom

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_StdSmoothCom_hh

