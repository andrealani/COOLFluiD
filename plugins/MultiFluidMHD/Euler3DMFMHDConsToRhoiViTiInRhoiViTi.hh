#ifndef COOLFluiD_Physics_MultiFluidMHD_Euler3DMFMHDConsToRhoiViTiInRhoiViTi_hh
#define COOLFluiD_Physics_MultiFluidMHD_Euler3DMFMHDConsToRhoiViTiInRhoiViTi_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

      class EulerMFMHDTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to primitive [p u v T] starting from primitive variables
 *
 * @author Alejandro Alvarez
 *
 */
class Euler3DMFMHDConsToRhoiViTiInRhoiViTi : public Framework::VarSetMatrixTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler3DMFMHDConsToRhoiViTiInRhoiViTi(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DMFMHDConsToRhoiViTiInRhoiViTi();

  /**
   * Set the transformation matrix from a given state
   */
  void setMatrix(const RealVector& state);

private:

  /**
   * Set the flag telling if the transformation is an identity one
   * @pre this method must be called during set up
   */
  bool getIsIdentityTransformation() const
  {
    return false;
  }

private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerMFMHDTerm> _model;

}; // end of class Euler3DMFMHDConsToRhoiViTiInRhoiViTi

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_Euler3DMFMHDConsToRhoiViTiInRhoiViTi_hh
