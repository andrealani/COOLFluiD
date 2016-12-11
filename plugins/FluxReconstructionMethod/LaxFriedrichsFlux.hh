#ifndef COOLFluiD_FluxReconstructionMethod_LaxFriedrichsFlux_hh
#define COOLFluiD_FluxReconstructionMethod_LaxFriedrichsFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a Lax-Friedrichs/Rusanov flux
 *
 */
class LaxFriedrichsFlux : public RiemannFlux {

public:  // methods

  /// Constructor
  LaxFriedrichsFlux(const std::string& name);

  /// Destructor
  ~LaxFriedrichsFlux();
 
  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState,
				  Framework::State& rState,
				  const RealVector& normal);

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal);


  /// Gets the Class name
  static std::string getClassName()
  {
    return "LaxFriedrichsFlux";
  }

  /// Set up private data and data
  virtual void setup();
  
  /// Unset up private data and data
  virtual void unsetup();

private: // data

  /// array storing the sum of the right and left flux
  RealVector  m_sumFlux;


}; // class LaxFriedrichsFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_LaxFriedrichsFlux_hh
