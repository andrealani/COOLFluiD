#ifndef COOLFluiD_Physics_StructMech_StructMech2DDiffusiveAxiDisp_hh
#define COOLFluiD_Physics_StructMech_StructMech2DDiffusiveAxiDisp_hh

//////////////////////////////////////////////////////////////////////////////

#include "StructMech2DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a StructMech physical model 2D for Dispitive
   * variables
   *
   * @author Thomas Wuilbaut
   */
class StructMech2DDiffusiveAxiDisp : public StructMech2DDiffusiveVarSet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  StructMech2DDiffusiveAxiDisp(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~StructMech2DDiffusiveAxiDisp();

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


  /**
   * Gets the axis of symmetry
   */
  std::string getAxisymmetryAxis()
  {
    cf_assert(false); return "";
  }

private: //member data

  /// plane strain model
  bool _isPlaneStress;

  /// Stiffness matrix
  RealMatrix _stiffness;

  /// constant lame coefs
  CFreal _lambda;
  CFreal _mu;

  /// Non Linear Flag
  bool _isNonLinear;

  /// StructMech is used for the movement of the mesh?
  bool _meshMovement;
  std::string _meshMovementMethod;

  //Which axis is the symmetry axis
  std::string m_axisymmetryAxis;

}; // end of class StructMech2DDiffusiveAxiDisp

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMech2DDiffusive_hh
