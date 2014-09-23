#ifndef COOLFluiD_Physics_Heat_Heat3DSourceTConst_hh
#define COOLFluiD_Physics_Heat_Heat3DSourceTConst_hh

//////////////////////////////////////////////////////////////////////////////

#include "Heat3DSourceVarSet.hh"
#include "MathTools/MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Heat3D source term depending on T
 *
 * @author Tiago Quintino
 */
class Heat3DSourceTConst : public Heat3DSourceVarSet {
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
  Heat3DSourceTConst(const std::string& name);

  /**
   * Default destructor
   */
  ~Heat3DSourceTConst();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    Heat3DSourceVarSet::configure(args);

    // add here configuration, specific of this class
  }

  /**
   * Checks if the source term has a part independent from the solution
   */
  bool hasIndepCoef() const
  {
    return MathTools::MathChecks::isNotZero(_indepCoef);
  }

  /**
   * Checks if the source term has a linear part
   */
  bool hasLinearCoef() const
  {
    return MathTools::MathChecks::isNotZero(_linearCoef);
  }

 /**
  * Gets the SourceFlux Coeficients
  */
  void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef);

 /**
  * Gets the SourceFlux Coeficients
  */
  void getIndepSourceCoefs(const Framework::State& state, RealVector& coef);

private: // member data

  /// independent coefficient
  CFreal _indepCoef;

  /// linear coefficient
  CFreal _linearCoef;

}; // end of class Heat3DSourceTConst

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat3DSourceFlux_hh
