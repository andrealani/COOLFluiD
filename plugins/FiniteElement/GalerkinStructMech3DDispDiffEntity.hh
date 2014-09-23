#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinStructMech3DDispDiffEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinStructMech3DDispDiffEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiffusiveEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"
#include "StructMech/StructMech3DDiffusiveDisp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Diffusive Term
   *
   * @author Tiago Quintino
   */

class GalerkinStructMech3DDispDiffEntity : public DiffusiveEntity {

public:

  /// Default constructor without arguments
  GalerkinStructMech3DDispDiffEntity(const std::string& name);

  /// Default destructor
  ~GalerkinStructMech3DDispDiffEntity();

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
  Common::SafePtr<Physics::StructMech::StructMech3DDiffusiveDisp> _structDiffVarSet;

  /// Anisotropy Flag
  bool _isAnisotropic;

  /// Non Linear Flag
  bool _isNonLinear;

  /// StructMech is used for the movement of the mesh?
  bool _meshMovement;
  std::string _meshMovementMethod;

  /// elements of the constitutive matrix
  CFreal c11;
  CFreal c12;
  CFreal c13;
  CFreal c21;
  CFreal c22;
  CFreal c23;
  CFreal c31;
  CFreal c32;
  CFreal c33;
  CFreal c41;
  CFreal c42;
  CFreal c43;
  CFreal c51;
  CFreal c52;
  CFreal c53;
  CFreal c61;
  CFreal c62;
  CFreal c63;
  CFreal c16;
  CFreal c14;
  CFreal c15;
  CFreal c26;
  CFreal c24;
  CFreal c25;
  CFreal c36;
  CFreal c34;
  CFreal c35;
  CFreal c46;
  CFreal c44;
  CFreal c45;
  CFreal c56;
  CFreal c54;
  CFreal c55;
  CFreal c66;
  CFreal c64;
  CFreal c65;

  /// constant lame coefs
  CFreal _lambda;
  CFreal _mu;

  /// Stiffness matrix
  RealMatrix _stiffness;

}; // end of class GalerkinStructMech3DDispDiffEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinStructMech3DDispDiffEntity_hh
