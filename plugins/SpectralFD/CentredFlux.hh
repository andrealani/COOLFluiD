#ifndef COOLFluiD_Numerics_SpectralFD_CentredFlux_hh
#define COOLFluiD_Numerics_SpectralFD_CentredFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/SpectralFDMethodData.hh"
#include "SpectralFD/RiemannFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a centred flux
 *
 * @author Kris Van den Abeele
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

private: // data

  /// array storing the sum of the right and left flux
  RealVector  m_sumFlux;

}; // class CentredFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_CentredFlux_hh
