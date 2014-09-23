#ifndef COOLFluiD_Physics_StructMech_MaterialData_hh
#define COOLFluiD_Physics_StructMech_MaterialData_hh


#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a library for MaterialData characteristics.
 *
 * @author Thomas Wuilbaut
 *
 */
class MaterialData{
public:

//   const MaterialData& operator=(const MaterialData& data)
//   {
//     young = data.young;
//     poisson = data.poisson;
//     density = data.density;
//     name = data.name;
//     C = data.C;
//     anisotropic = data.anisotropic;
//     return *this;
//   }

/// Name
std::string name;

/// Anisotropy
bool anisotropic;

/// Consistitutive relation matrix
/// stress = C * strain
RealMatrix C;

/// Young Modulus
CFreal young;

/// Poisson Coef
CFreal poisson;

/// Density
CFreal density;

/// Thermal Expansion Coef
CFreal alpha;

}; // end of class MaterialData

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_MaterialData_hh
