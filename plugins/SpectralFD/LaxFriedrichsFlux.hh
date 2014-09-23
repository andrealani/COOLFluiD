#ifndef COOLFluiD_Numerics_SpectralFD_LaxFriedrichsFlux_hh
#define COOLFluiD_Numerics_SpectralFD_LaxFriedrichsFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/SpectralFDMethodData.hh"
#include "SpectralFD/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

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

private: // data

  /// array storing the sum of the right and left flux
  RealVector  m_sumFlux;


}; // class LaxFriedrichsFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_LaxFriedrichsFlux_hh
