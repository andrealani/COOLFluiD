#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDiffusiveDisp_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDiffusiveDisp_hh

//////////////////////////////////////////////////////////////////////////////

#include "StructMechHeat2DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a StructMechHeat physical model 2D for Dispitive
   * variables
   *
   * @author Thomas Wuilbaut
   */
class StructMechHeat2DDiffusiveDisp : public StructMechHeat2DDiffusiveVarSet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  StructMechHeat2DDiffusiveDisp(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~StructMechHeat2DDiffusiveDisp();

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
   * Gets if the model is used for mesh movement
   */
  std::string getMeshMovementMethod()
  {
    return _meshMovementMethod;
  }

  /**
   * Gets if the model is used for nonlinear struct mech
   */
  bool isNonLinear()
  {
    return _isNonLinear;
  }

  /**
   * Gets the Stiffness Matrix
   */
  RealMatrix& getStiffnessMat();

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

private: //member data

  /// plane strain model
  bool _isPlaneStress;

  /// constant of the constitutive matrix
  CFreal _c11;
  CFreal _c12;
  CFreal _c21;
  CFreal _c22;
  CFreal _c16;
  CFreal _c26;
  CFreal _c61;
  CFreal _c62;
  CFreal _c66;

  /// Stiffness matrix
  RealMatrix _stiffness;

  /// constant lame coefs
  CFreal _lambda;
  CFreal _mu;

  /// Non Linear Flag
  bool _isNonLinear;

  /// StructMechHeat is used for the movement of the mesh?
  bool _meshMovement;
  std::string _meshMovementMethod;

}; // end of class StructMechHeat2DDiffusiveDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMechHeat_StructMechHeat2DDiffusive_hh
