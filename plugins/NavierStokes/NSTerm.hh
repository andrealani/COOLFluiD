#ifndef COOLFluiD_Physics_NavierStokes_NSTerm_hh
#define COOLFluiD_Physics_NavierStokes_NSTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"
#include "NavierStokes/DynamicViscosityLaw.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokesModel.
 *
 * @author Andrea Lani
 *
 */
class NSTerm : public Framework::BaseTerm {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {MU=0, LAMBDA=1, RE=2};

  /**
   * Constructor without arguments
   */
  NSTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NSTerm();

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 3;
  }
  
  /**
   * Compute the dynamic viscosity starting from the
   * DIMENSIONAL pressure and temperature
   */
  CFreal getDynViscosityDim(const CFreal& pdim,
			    const CFreal& Tdim)
  {
    cf_assert(_dynViscosity.isNotNull());
    return _dynViscosity->compute(pdim, Tdim);
  }

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();
  
  /**
   * Get constant thermal conductivity for incompressible flow
   */
  CFreal getThermConductivity() const
  {
    return _thermConductivity;
  }
  
  /**
   * Get the Prandtl number
   */
  CFreal getPrandtl() const
  {
    return _Prandtl;
  }

  /**
   * Get the Reynolds number
   */
  CFreal getReynolds() const
  {
    return _ReynoldsRef;
  }

  /**
   * Set the coefficient for the scaling of the viscous stresses term
   */
  void setCoeffTau(const CFreal coeffTau)
  {
    _coeffTau = coeffTau;
  }

  /**
   * Set the coefficient for the scaling of the heat term
   */
  void setCoeffQ(const CFreal coeffQ)
  {
    _coeffQ = coeffQ;
  }

  /**
   * Set the cp over Prandtl coefficient
   */
  void setCpOverPrandtl(const CFreal cpOverPr)
  { 
    cf_assert(cpOverPr > 0.0);
    _cpOverPrandtl = cpOverPr;
    cf_assert(_cpOverPrandtl > 0.0);
  }
  
  void setPrandtlTurb(const CFreal Pr)
  {
    _PrandtlTurb = Pr;
  }

  /**
   * Get the coefficient for the scaling of the viscous stresses term
   */
  CFreal getCoeffTau() const
  {
    cf_assert(_coeffTau > 0.0);
    return _coeffTau;
  }

  /**
   * Get the coefficient for the scaling of the heat term
   */
  CFreal getCoeffQ() const
  {
    cf_assert(_coeffQ > 0.0);
    return _coeffQ;
  }

  /**
   * Get Cp over Prandtl
   */
  CFreal getCpOverPrandtl() const
  {
    cf_assert(_cpOverPrandtl > 0.0);
    return _cpOverPrandtl;
  }

  /**
   * Get ratio of Prandtl numbers (laminar and turbulent)
   */
  CFreal getPrOverPrT() const
  {
    cf_assert(_Prandtl > 0.0);
    cf_assert(_PrandtlTurb > 0.0);
    return _Prandtl / _PrandtlTurb;
  }


  CFreal getPrandtlTurb() const
  {
//     cf_assert(_PrandtlTurb > 0.0);
    return _PrandtlTurb;
  }
  
  /**
   * Get the artificial viscosity coefficient
   */
  CFreal getArtDiffCoeff() const
  {
    cf_assert(_muDiff > 0.0);
    return _muDiff;
  }
  
  /**
   * Get the name
   */
  static std::string getName()
  {
    return "NSTerm";
  }

private:
  
  /// dynamic viscosity law
  Common::SelfRegistPtr<DynamicViscosityLaw> _dynViscosity;
  
  /// adimensional coefficient
  CFreal _coeffTau;

  /// adimensional coefficient
  CFreal _coeffQ;

  /// prandtl over Cp
  CFreal _cpOverPrandtl;
  
  // turbulent prandtl 
  CFreal _PrandtlTurb;
  
  /// reference Reynolds number
  CFreal _ReynoldsRef;
  
  /// Prandtl number
  CFreal _Prandtl;
  
  /// name of the dynamic viscosity law
  std::string _dynViscosityStr;

  /// constant thermal conductivity for incompressible flow
  CFreal _thermConductivity;

  /// coefficient for increasing artificially the dissipation
  CFreal _muDiff;
  
}; // end of class NSTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NSTerm_hh
