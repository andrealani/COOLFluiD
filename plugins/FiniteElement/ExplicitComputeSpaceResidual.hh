#ifndef COOLFluiD_Numerics_FiniteElement_ExplicitComputeSpaceResidual_hh
#define COOLFluiD_Numerics_FiniteElement_ExplicitComputeSpaceResidual_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSpaceResidual.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to ComputeSpaceResidual in an explicit manner.
 */
class ExplicitComputeSpaceResidual : public ComputeSpaceResidual {
public:

  /**
   * Constructor.
   */
  explicit ExplicitComputeSpaceResidual(const std::string& name);

  /**
   * Destructor.
   */
  ~ExplicitComputeSpaceResidual();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Execute Processing actions
   */
  void executeOnTrs();

private: // data

}; // class ExplicitComputeSpaceResidual

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ExplicitComputeSpaceResidual_hh

