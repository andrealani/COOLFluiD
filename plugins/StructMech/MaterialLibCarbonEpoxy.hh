#ifndef COOLFluiD_Physics_StructMech_MaterialLibCarbonEpoxy_hh
#define COOLFluiD_Physics_StructMech_MaterialLibCarbonEpoxy_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MatrixInverterT.hh"
#include "StructMech/MaterialPropertyLib.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the properties of CarbonEpoxy
 *
 * @author Thomas Wuilbaut
 *
 */
class MaterialLibCarbonEpoxy : public MaterialPropertyLib
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
  MaterialLibCarbonEpoxy(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MaterialLibCarbonEpoxy();

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

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
    return true;
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
    return &m_C;
  }


protected:

  /**
   * Transform the stiffness matrix with the correct angles
   */
  void stiffnessMatrixTransform();

protected:

  /// Anisotropy
  bool m_isAnisotropic;

  /// Consistitutive relation matrix
  /// stress = C * strain
  RealMatrix m_C;

  ///Transformation Matrices
  RealMatrix m_T;
  RealMatrix m_Tt;

  /// Young Modulus
  CFreal m_young;

  /// Poisson Coef
  CFreal m_poisson;

  /// Density
  CFreal m_density;

  /// Thermal Expansion Coef
  CFreal m_alpha;

  /// Angle between material and global (for anisotropic materials)
  RealVector m_materialAngles;

  /// matrix inverter size 3
  MathTools::MatrixInverterT<3> m_inverter3;

}; // end of class MaterialLibCarbonEpoxy

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_MaterialLibCarbonEpoxy_hh
