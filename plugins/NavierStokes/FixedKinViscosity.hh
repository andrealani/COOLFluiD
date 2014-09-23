#ifndef COOLFluiD_Physics_NavierStokes_FixedKinViscosity_hh
#define COOLFluiD_Physics_NavierStokes_FixedKinViscosity_hh

//////////////////////////////////////////////////////////////////////////////

#include "DynamicViscosityLaw.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

class Euler2DVarSet;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the dynamic viscosity starting from a constant kinematic viscosity
 *
 * @author Kris Van den Abeele
 *
 */
class FixedKinViscosity : public DynamicViscosityLaw
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
  FixedKinViscosity(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FixedKinViscosity();

  /**
   * Compute the dynamic viscosity starting from the dimensional
   * pressure and temperature
   */
  CFreal compute(const CFreal& pdim, const CFreal& Tdim);

protected: // data

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_eulerVarSet;

  /// Fixed kinematic viscosity
  CFreal m_fixedKinVisc;

}; // end of class FixedKinViscosity

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_SutherlandDynViscosity_hh
