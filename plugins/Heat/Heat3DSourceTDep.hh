#ifndef COOLFluiD_Physics_Heat_Heat3DSourceTDep_hh
#define COOLFluiD_Physics_Heat_Heat3DSourceTDep_hh

//////////////////////////////////////////////////////////////////////////////

#include "Heat3DSourceVarSet.hh"
#include "Framework/VectorialFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Source term definied relative to the Heat3DPrim varset
   *
   * @author Tiago Quintino
   */
class Heat3DSourceTDep : public Heat3DSourceVarSet {
public: // classes

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   * @see Heat3DSourceFlux
   */
  Heat3DSourceTDep(const std::string& name);

  /**
   * Default destructor
   */
  ~Heat3DSourceTDep();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Setup the member data
   */
  virtual void setup();

  /**
   * Evaluate the coefficients for the Linear VectorialFunction Source term
   */
  void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef);

  /**
   * Evaluate the coefficients for the Independent VectorialFunction Source term
   */
  void getIndepSourceCoefs(const Framework::State& state, RealVector& coef);

  /**
   * Checks if the source term has a part independent from the solution
   */
  bool hasIndepCoef() const
  {
    return !_functionsIndep.empty();
  }

  /**
   * Checks if the source term has a linear part
   */
  bool hasLinearCoef() const
  {
    return !_functionsLinear.empty();
  }

private: // member data

  /// a vector of string to hold the independent coefficient functions
  std::vector<std::string> _functionsIndep;

  /// a vector of string to hold the linear coefficient functions
  std::vector<std::string> _functionsLinear;

  /// a vector of string to hold the variables
  std::vector<std::string> _varsDescript;

  /// the VectorialFunction for independent coefficient
  Framework::VectorialFunction _vFunctionIndep;

  /// the VectorialFunction for linear coefficient
  Framework::VectorialFunction _vFunctionLinear;

  /// placeholder for variable values
  RealVector _vars;

  /// placeholder for evaluate coefficients which are matrices
  /// needed because VectorialFunction only uses RealVector
  RealVector _vectorialCoef;

}; // end of class Heat3DSourceTDep

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat3DSourceFlux_hh
