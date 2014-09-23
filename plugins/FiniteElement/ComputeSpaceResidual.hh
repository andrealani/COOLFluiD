#ifndef COOLFluiD_Numerics_FiniteElement_ComputeSpaceResidual_hh
#define COOLFluiD_Numerics_FiniteElement_ComputeSpaceResidual_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Common/CFMap.hh"
#include "FElemTypeData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to compute the space residual
 */
class ComputeSpaceResidual : public FiniteElementMethodCom {
public:

  /**
   * Constructor.
   */
  explicit ComputeSpaceResidual(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ComputeSpaceResidual();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup the private data and data of the aggregated classes
   * in this command after the  processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  /// map of LSSMatrix accumulators, one for each cell type
  Common::CFMap<CFuint,FElemTypeData> m_map_femdata;

  /// socket for Rhs
  Framework::DataSocketSink< CFreal> socket_rhs;

}; // class ComputeSpaceResidual

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ComputeSpaceResidual_hh

