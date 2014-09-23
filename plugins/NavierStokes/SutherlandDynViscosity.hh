#ifndef COOLFluiD_Physics_NavierStokes_SutherlandDynViscosity_hh
#define COOLFluiD_Physics_NavierStokes_SutherlandDynViscosity_hh

//////////////////////////////////////////////////////////////////////////////

#include "DynamicViscosityLaw.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the dynamic viscosity with the Sutherland law
 *
 * @author Tiago Quintino
 *
 */
class SutherlandDynViscosity : public DynamicViscosityLaw
{
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SutherlandDynViscosity(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SutherlandDynViscosity();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Compute the dynamic viscosity starting from the dimensional
   * pressure and temperature
   */
  CFreal compute(const CFreal& pdim, const CFreal& Tdim);

protected: // data

  /// Compute the viscosity for this gas
  std::string m_gas;

  /// Sutherland constant for the fluid
  CFreal m_SuthConst;

  /// Reference viscosity
  CFreal m_ViscRef;

  /// Reference temperature at which the reference viscosity is measured
  CFreal m_TRef;

}; // end of class SutherlandDynViscosity

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_SutherlandDynViscosity_hh
