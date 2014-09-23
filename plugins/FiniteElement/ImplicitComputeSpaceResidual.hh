#ifndef COOLFluiD_Numerics_FiniteElement_ImplicitComputeSpaceResidual_hh
#define COOLFluiD_Numerics_FiniteElement_ImplicitComputeSpaceResidual_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSpaceResidual.hh"
#include "ComputeJacobStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to ComputeSpaceResidual in an implicit manner.
 */
class ImplicitComputeSpaceResidual : public ComputeSpaceResidual {
public:

  /**
   * Constructor.
   */
  explicit ImplicitComputeSpaceResidual(const std::string& name);

  /**
   * Destructor.
   */
  ~ImplicitComputeSpaceResidual();

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

}; // class ImplicitComputeSpaceResidual

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ImplicitComputeSpaceResidual_hh

