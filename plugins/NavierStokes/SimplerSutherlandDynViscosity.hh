#ifndef COOLFluiD_Physics_NavierStokes_SimplerSutherlandDynViscosity_hh
#define COOLFluiD_Physics_NavierStokes_SimplerSutherlandDynViscosity_hh

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
 * @author Tiago Quintino
 *
 */
class SimplerSutherlandDynViscosity : public DynamicViscosityLaw
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
  SimplerSutherlandDynViscosity(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SimplerSutherlandDynViscosity();

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

  /// Reference viscosity
  CFreal m_ViscRef;

  /// Reference temperature
  CFreal m_TRef;

}; // end of class SimplerSutherlandDynViscosity

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_SimplerSutherlandDynViscosity_hh
