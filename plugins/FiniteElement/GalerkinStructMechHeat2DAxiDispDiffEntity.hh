#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DAxiDispDiffEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DAxiDispDiffEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiffusiveEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"
#include "StructMechHeat/StructMechHeat2DDiffusiveAxiDisp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Diffusive Term
   *
   * @author Thomas Wuilbaut
   */

class GalerkinStructMechHeat2DAxiDispDiffEntity : public DiffusiveEntity {

public:

  /// Default constructor without arguments
  GalerkinStructMechHeat2DAxiDispDiffEntity(const std::string& name);

  /// Default destructor
  ~GalerkinStructMechHeat2DAxiDispDiffEntity();

  /**
   * Setup the strategy
   */
  virtual void setup();

  /**
   * Overloading of operator()
   */
  virtual RealMatrix& operator()();

private: // helper functions

  /**
   * Modify the Stiffness for improving mesh movement
   */
  void modifyStiffness(Framework::GeometricEntity* const geo);

  /**
   * Resets the Stiffness coefs to their original value
   */
  void resetStiffness();

private:

  /// Variable Set
  Common::SafePtr<Physics::StructMechHeat::StructMechHeat2DDiffusiveAxiDisp> _structDiffVarSet;

  /// Anisotropy Flag
  bool _isAnisotropic;

  /// Non Linear Flag
  bool _isNonLinear;

  /// StructMechHeat is used for the movement of the mesh?
  bool _meshMovement;
  std::string _meshMovementMethod;

  /// Stiffness matrix
  RealMatrix _stiffness;

  /// What coordinate is R and what is Z
  CFuint m_radiusID;
  CFuint m_zID;

}; // end of class GalerkinStructMechHeat2DAxiDispDiffEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DAxiDispDiffEntity_hh
