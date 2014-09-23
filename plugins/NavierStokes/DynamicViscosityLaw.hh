#ifndef COOLFluiD_Physics_NavierStokes_DynamicViscosityLaw_hh
#define COOLFluiD_Physics_NavierStokes_DynamicViscosityLaw_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ConcreteProvider.hh"
#include "Common/OwnedObject.hh"
#include "Common/NonCopyable.hh"
#include "Config/ConfigObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class offers a basic interface to compute the dynamic viscosity
 * with a certain law
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class DynamicViscosityLaw : public Common::NonCopyable<DynamicViscosityLaw>,
                            public Common::OwnedObject,
                            public Config::ConfigObject
{
public: // functions

  /// the provider of this type of classes
  typedef Environment::ConcreteProvider<DynamicViscosityLaw,1> PROVIDER;

  /// the first argument in the creation should be the name
  typedef const std::string& ARG1;

  /**
   * Constructor
   */
  DynamicViscosityLaw(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DynamicViscosityLaw();

  /**
   * Compute the dynamic viscosity starting from the dimensional
   * pressure and temperature
   */
  virtual CFreal compute(const CFreal& pdim, const CFreal& Tdim) = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "DynamicViscosityLaw";
  }

}; // end of class DynamicViscosityLaw

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_DynamicViscosityLaw_hh
