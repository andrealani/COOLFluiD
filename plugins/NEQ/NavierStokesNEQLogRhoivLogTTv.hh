#ifndef COOLFluiD_Physics_NEQ_NavierStokesNEQLogRhoivLogTTv_hh
#define COOLFluiD_Physics_NEQ_NavierStokesNEQLogRhoivLogTTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQRhoivt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * Diffusive variable set for NEQ with (ln rho_i, u, v, ln T, ln Tv) states:
 * logarithms of the partial densities rho_i, of T and of every Tv.
 *
 * The gradient variables are (ln y_i, u, v, ln T, ln Tv), with y_i the mass
 * fractions, so their traces are exp of polynomials like the stored
 * variables. Every consumer of the gradients gets the physical ones back:
 * grad y_i = y_i grad(ln y_i) and grad T = T grad(ln T), with y_i and T at
 * the evaluation point. The diffusion driving forces then vanish with y_i, so
 * the diffusion velocity J_i/rho_i stays bounded where a species goes to zero.
 *
 * @author Rayan Dhib
 */
template <class BASE>
class NavierStokesNEQLogRhoivLogTTv : public NavierStokesNEQRhoivt<BASE> {

public:

  /// Constructor
  NavierStokesNEQLogRhoivLogTTv(const std::string& name,
                                Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  virtual ~NavierStokesNEQLogRhoivLogTTv();

  /// Set up the private data
  virtual void setup();

  /// Set the composition in the library from the state
  virtual void setComposition(const RealVector& state,
                              const bool isPerturb,
                              const CFuint iVar);

  /// Gradient variables (ln y_i, u, v, ln T, ln Tv) of the states, one column per state
  virtual void setGradientVars(const std::vector<RealVector*>& states,
                               RealMatrix& values,
                               const CFuint stateSize);

  /// Mixture density, sum of exp(ln rho_i)
  virtual CFreal getDensity(const RealVector& state);

  /// Dynamic viscosity at T, Tv = exp(ln T), exp(ln Tv)
  virtual CFreal getDynViscosity(const RealVector& state,
                                 const std::vector<RealVector*>& gradients);

  /**
   * Transport properties. The base reads Tv straight from the state and
   * dy_i/dn from the gradients, so it gets a copy with physical temperatures
   * and the mass fraction gradients.
   */
  virtual void computeTransportProperties(const RealVector& state,
                                          const std::vector<RealVector*>& gradients,
                                          const RealVector& normal);

  using NavierStokesNEQRhoivt<BASE>::getFlux;

  /// Diffusive flux, with the temperature gradients rebuilt first
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);

  /// Heat flux, with the mass fraction and temperature gradients rebuilt first
  virtual CFreal getHeatFlux(const RealVector& state,
                             const std::vector<RealVector*>& gradients,
                             const RealVector& normal);

  /**
   * Source of the axisymmetric equations (planar form times r), with the mass fraction and
   * temperature gradients rebuilt first. With one Tv and no Te equation, the diffusion of the
   * Tv-mode energy of atoms and electrons (sum of hsEl J_i over the species) is added to the
   * Tv equation, as in getFlux; the base source has only the molecules' hsVib J_i.
   */
  virtual void getAxiSourceTerm(const RealVector& physicalData,
                                const RealVector& state,
                                const std::vector<RealVector*>& gradients,
                                const CFreal& radius,
                                RealVector& source);

protected:

  /// Set the gradient state (y_i, u, v, T, Tv), the library state and the pressure
  virtual void setGradientState(const RealVector& state);

private:

  /// fills m_rhoi with exp(ln rho_i) of state and returns the mixture density
  CFreal setPartialDensities(const RealVector& state);

  /// gradient list with the species entries replaced by grad y_i = y_i grad(ln y_i), y_i = _gradState[i]
  const std::vector<RealVector*>& setMassFractionGradients(const std::vector<RealVector*>& gradients);

  /// gradient list with the T and Tv entries replaced by grad T = T grad(ln T), T = _gradState[i]
  const std::vector<RealVector*>& setTemperatureGradients(const std::vector<RealVector*>& gradients);

private:

  /// copy of the state with physical temperatures
  RealVector m_physTState;

  /// grad y_i of the species
  std::vector<RealVector> m_gradY;

  /// gradient list with grad y_i
  std::vector<RealVector*> m_gradYPtrs;

  /// grad T and grad Tv
  std::vector<RealVector> m_gradTemps;

  /// gradient list with grad T and grad Tv
  std::vector<RealVector*> m_gradTempPtrs;

  /// partial densities of the current state
  RealVector m_rhoi;

  /// (T, Tv) of the current state, passed to the library with m_rhoi
  RealVector m_temps;

}; // end of class NavierStokesNEQLogRhoivLogTTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQLogRhoivLogTTv.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_NavierStokesNEQLogRhoivLogTTv_hh
