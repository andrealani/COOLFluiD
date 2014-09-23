#ifndef COOLFluiD_Physics_StructMech_MaterialLibCSiC_hh
#define COOLFluiD_Physics_StructMech_MaterialLibCSiC_hh

//////////////////////////////////////////////////////////////////////////////

#include "MaterialPropertyLib.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the properties of a CSiC Material
 *
 * @author Thomas Wuilbaut
 *
 */
class MaterialLibCSiC : public MaterialPropertyLib
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
  MaterialLibCSiC(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MaterialLibCSiC();

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
  virtual CFreal computeDensity()
  {
    return m_density;
  }

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

  /// Young Modulus
  CFreal m_young;

  /// Poisson Coef
  CFreal m_poisson;

  /// Density
  CFreal m_density;

  /// Thermal Expansion Coef
  CFreal m_alpha;

}; // end of class MaterialLibCSiC

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_MaterialLibCSiC_hh
