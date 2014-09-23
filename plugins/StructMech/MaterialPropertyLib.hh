#ifndef COOLFluiD_Physics_StructMech_MaterialPropertyLib_hh
#define COOLFluiD_Physics_StructMech_MaterialPropertyLib_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ConcreteProvider.hh"
#include "Common/OwnedObject.hh"
#include "Common/NonCopyable.hh"
#include "Config/ConfigObject.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class offers a basic interface to compute the characteristics of solid materials
 *
 * @author Thomas Wuilbaut
 *
 */
class MaterialPropertyLib : public Common::NonCopyable<MaterialPropertyLib>,
                            public Common::OwnedObject,
                            public Config::ConfigObject
{
public: // functions

  /// the provider of this type of classes
  typedef Environment::ConcreteProvider<MaterialPropertyLib,1> PROVIDER;

  /// the first argument in the creation should be the name
  typedef const std::string& ARG1;

  /**
   * Constructor
   */
  MaterialPropertyLib(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MaterialPropertyLib();

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args ){
    Config::ConfigObject::configure(args);
  }

  /**
   * Compute the Young Modulus
   */
  virtual CFreal computeYoungModulus() = 0;

  /**
   * Compute the Poisson Coeficient
   */
  virtual CFreal computePoissonCoef() = 0;

  /**
   * Compute the Density
   */
  virtual CFreal computeDensity() = 0;

  /**
   * Is the material anisotropic
   */
  virtual bool isAnisotropic() = 0;

  /**
   * Compute the thermal Expansion Coeficient
   */
  virtual CFreal computeThermalExpansionCoef() = 0;

  /**
   * Compute Anisotropic Matrix
   */
  virtual RealMatrix* computeStiffnessMatrix() = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "MaterialPropertyLib";
  }

protected: // data

}; // end of class MaterialPropertyLib

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_MaterialPropertyLib_hh
