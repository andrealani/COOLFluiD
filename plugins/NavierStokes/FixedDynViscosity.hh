#ifndef COOLFluiD_Physics_NavierStokes_FixedDynViscosity_hh
#define COOLFluiD_Physics_NavierStokes_FixedDynViscosity_hh

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
 * @author Andrea Lani
 *
 */
class FixedDynViscosity : public DynamicViscosityLaw
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
  FixedDynViscosity(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FixedDynViscosity();

  /**
   * Compute the dynamic viscosity starting from the dimensional
   * pressure and temperature
   */
  CFreal compute(const CFreal& pdim, const CFreal& Tdim);

protected: // data

  /// Fixed viscosity
  CFreal m_FixedVisc;

}; // end of class FixedDynViscosity

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_SutherlandDynViscosity_hh
