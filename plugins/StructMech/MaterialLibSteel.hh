#ifndef COOLFluiD_Physics_StructMech_MaterialLibSteel_hh
#define COOLFluiD_Physics_StructMech_MaterialLibSteel_hh

//////////////////////////////////////////////////////////////////////////////

#include "MaterialPropertyLib.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the properties of Steel
 *
 * @author Thomas Wuilbaut
 *
 */
class MaterialLibSteel : public MaterialPropertyLib
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
  MaterialLibSteel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MaterialLibSteel();

  /**
   * Compute the Young Modulus
   */
  virtual CFreal computeYoungModulus();

  /**
   * Compute the Poisson Coeficient
   */
  virtual CFreal computePoissonCoef();

  /**
   * Compute the Density
   */
  virtual CFreal computeDensity();

  /**
   * Is the material anisotropic
   */
  virtual bool isAnisotropic()
  {
    return false;
  }

  /**
   * Compute the thermal Expansion Coeficient
   */
  virtual CFreal computeThermalExpansionCoef();

  /**
   * Compute Anisotropic Matrix
   */
  virtual RealMatrix* computeStiffnessMatrix()
  {
    cf_assert(false);
    return CFNULL;
  }

protected:

  /// Anisotropy
  bool m_isAnisotropic;

  /// Consistitutive relation matrix
  /// stress = C * strain
  RealMatrix C;

  /// Young Modulus
  CFreal m_young;

  /// Poisson Coef
  CFreal m_poisson;

  /// Density
  CFreal m_density;

  /// Thermal Expansion Coef
  CFreal m_alpha;

}; // end of class MaterialLibSteel

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_MaterialLibSteel_hh
