#ifndef COOLFluiD_Physics_EntropyMHD_EntropyMHDTerm_hh
#define COOLFluiD_Physics_EntropyMHD_EntropyMHDTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * Physical term for the EntropyMHD 3D MHD model with GLM projection.
 *
 * Entropy-variable MHD: eq 7 evolves sigma = p*rho^(1-gamma) (not total energy E).
 * State variables (primitive): (rho, u, v, w, Bx, By, Bz, sigma, phi)
 * Flux[7] = sigma*Vn, Source[7] = 0 (pure conservation law).
 * Pressure recovered as p = sigma*rho^(gamma-1). Avoids catastrophic cancellation at low beta.
 *
 * @author Rayan Dhib
 */
class EntropyMHDTerm : public Framework::BaseTerm {
public:

  /**
   * Physical-data slot layout. T holds sigma = p*rho^(1-gamma); P is the recovered
   * thermodynamic pressure.
   */
  enum {RHO=0, P=1, BX=2, BY=3, BZ=4, VX=5, VY=6, VZ=7, GAMMA=8,
        XP=9, YP=10, ZP=11, PHI=12, T=13};

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  EntropyMHDTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~EntropyMHDTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 14;
  }

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "EntropyMHDTerm";
  }

  /// @name Accessors
  //@{
  CFreal getGamma() const { return _gamma; }
  CFreal getRefSpeed() const { return _refSpeed; }
  CFreal getMass() const { return _mass; }
  CFreal getNRef() const { return _nRef; }
  CFreal getBRef() const { return _BRef; }
  CFreal getLRef() const { return _lRef; }
  CFreal getTRef() const { return _TRef; }

  /// Thermodynamic parameters
  CFreal getSpitzerKappa0() const { return _kappa0; }
  CFreal getHeatingH0() const { return _heatingH0; }
  CFreal getHeatingLambda() const { return _heatingLambda; }
  CFreal getTRBroadeningTcut() const { return _TRBroadeningTcut; }

  /// Non-dimensional gravity factor: g_nondim = -gravFactor/r^2 * r_hat
  CFreal getGravFactorNonDim() const { return _gravFactorNonDim; }

  /// Equilibrium pressure at non-dimensional radius r (for initial conditions)
  CFreal getEquilibriumPressure(const CFreal r, const CFreal rho0,
                                const CFreal p0) const
  {
    cf_assert(r > 0.0);
    const CFreal T0 = p0 / rho0;
    const CFreal alpha = (_gamma - 1.0) / _gamma * _gravFactorNonDim / T0;
    const CFreal Tratio = 1.0 - alpha * (1.0 - 1.0/r);
    return p0 * std::pow(Tratio, _gamma / (_gamma - 1.0));
  }

  /// Equilibrium density at non-dimensional radius r (for initial conditions)
  CFreal getEquilibriumDensity(const CFreal r, const CFreal rho0,
                               const CFreal p0) const
  {
    cf_assert(r > 0.0);
    const CFreal T0 = p0 / rho0;
    const CFreal alpha = (_gamma - 1.0) / _gamma * _gravFactorNonDim / T0;
    const CFreal Tratio = 1.0 - alpha * (1.0 - 1.0/r);
    return rho0 * std::pow(Tratio, 1.0 / (_gamma - 1.0));
  }
  //@}

protected:

  /// gamma (5/3 for thermodynamic, 1.05 for polytropic fallback)
  CFreal _gamma;

  /// reference speed for the projection scheme
  CFreal _refSpeed;

  /// mass of the external object (default: solar mass)
  CFreal _mass;

  /// reference quantities for non-dimensionalization
  CFreal _nRef;   ///< reference proton density [m^-3]
  CFreal _BRef;   ///< reference magnetic field [T]
  CFreal _lRef;   ///< reference length (R_sun) [m]
  CFreal _TRef;   ///< reference temperature [K]

  /// Thermodynamic parameters (for future full-energy model)
  CFreal _kappa0;              ///< Spitzer conductivity [W m^-1 K^-7/2]
  CFreal _heatingH0;           ///< Heating amplitude [SI]
  CFreal _heatingLambda;       ///< Heating scale height [R_sun]
  CFreal _TRBroadeningTcut;    ///< TR broadening cutoff temperature [K]

  /// Non-dimensional gravity factor (set in setupPhysicalData)
  CFreal _gravFactorNonDim;    ///< GM / (lRef * V0^2)

}; // end of class EntropyMHDTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_EntropyMHD_EntropyMHDTerm_hh
