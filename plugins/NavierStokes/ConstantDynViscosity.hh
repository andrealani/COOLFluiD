#ifndef COOLFluiD_Physics_NavierStokes_ConstantDynViscosity_hh
#define COOLFluiD_Physics_NavierStokes_ConstantDynViscosity_hh

//////////////////////////////////////////////////////////////////////////////

#include "DynamicViscosityLaw.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes a constant dynamic viscosity
 *
 * @author Andrea Lani
 *
 */
class ConstantDynViscosity : public DynamicViscosityLaw
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
  ConstantDynViscosity(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ConstantDynViscosity();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Compute the dynamic viscosity starting from the dimensional
   * pressure and temperature
   */
  CFreal compute(const CFreal& pdim, const CFreal& Tdim);
  
protected: // data

  /// value
  CFreal m_value;
  
}; // end of class ConstantDynViscosity

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_ConstantDynViscosity_hh
