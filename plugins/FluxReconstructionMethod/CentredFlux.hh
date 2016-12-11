#ifndef COOLFluiD_FluxReconstructionMethod_VCJH_hh
#define COOLFluiD_FluxReconstructionMethod_VCJH_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Vincent-Castonguay-Jameson-Huyn Correction Function
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 */
class CentredFlux : public RiemannFlux {

public:  // methods

  /// Constructor
  CentredFlux(const std::string& name);

  /// Destructor
  ~CentredFlux();

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal);
                                  
  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState,
                                  Framework::State& rState,
                                  const RealVector& normal);


  /// Gets the Class name
  static std::string getClassName()
  {
    return "CentredFlux";
  }

  /// Set up private data and data
  virtual void setup();
  
  /// Unset up private data and data
  virtual void unsetup();

private: // data

  /// array storing the sum of the right and left flux
  RealVector  m_sumFlux;

}; // class CentredFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_CentredFlux_hh
