#ifndef COOLFluiD_Physics_StructMech_StructMech3DDiffusiveDisp_hh
#define COOLFluiD_Physics_StructMech_StructMech3DDiffusiveDisp_hh

//////////////////////////////////////////////////////////////////////////////

#include "StructMech3DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a StructMech physical model 3D for Dispitive
   * variables
   *
   * @author Tiago Quintino
   */
class StructMech3DDiffusiveDisp : public StructMech3DDiffusiveVarSet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  StructMech3DDiffusiveDisp(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~StructMech3DDiffusiveDisp();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Gets if the model is used for mesh movement
   */
  bool isMeshMovement()
  {
    return _meshMovement;
  }

  /**
   * Gets if the model is linear elastic or nonlinear (geometrically)
   */
  bool isNonLinear()
  {
    return _isNonLinear;
  }

  /**
   * Gets if the model is used for mesh movement
   */
  std::string getMeshMovementMethod()
  {
    return _meshMovementMethod;
  }

  /**
   * Gets the stiffness matrix
   */
  RealMatrix& getStiffnessMat()
  {
    return _stiffness;
  }

  /**
   * Gets the lambda coef
   * To be removed!!!
   */
  CFreal getLambdaCoef()
  {
    return _lambda;
  }

  /**
   * Gets the mu coef
   * To be removed!!!
   */
  CFreal getMuCoef()
  {
    return _mu;
  }


private:

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

  RealMatrix _stiffness;

}; // end of class StructMech3DDiffusiveDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech3DDiffusive_hh
