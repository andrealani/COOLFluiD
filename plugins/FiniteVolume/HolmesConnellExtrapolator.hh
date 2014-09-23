#ifndef COOLFluiD_Numerics_FiniteVolume_HolmesConnellExtrapolator_hh
#define COOLFluiD_Numerics_FiniteVolume_HolmesConnellExtrapolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NodalStatesExtrapolator.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object
 *
 * @author Andrea Lani
 *
 */
class HolmesConnellExtrapolator : public Framework::NodalStatesExtrapolator<CellCenterFVMData> {
public:

  /**
   * Constructor
   */
  HolmesConnellExtrapolator(const std::string& name);

  /**
   * Default destructor
   */
  ~HolmesConnellExtrapolator();

  /**
   * Set up private data needed by the computation
   */
  void setup();

  /**
   * Extrapolate the solution in all mesh nodes
   */
  void extrapolateInAllNodes();

  /**
   * Extrapolate the solution in the given nodes
   */
  void extrapolateInNodes(const std::vector<Framework::Node*>& nodes);

private: // helper function

  /**
   * Compute all the reconstruction coefficients
   */
  void computeCoeffs();

private:

  /// storage of the reconstruction coefficients
  RealVector _coeffs;

}; // end of class HolmesConnellExtrapolator

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HolmesConnellExtrapolator_hh
