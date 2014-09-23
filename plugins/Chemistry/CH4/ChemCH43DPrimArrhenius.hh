#ifndef COOLFluiD_Physics_Chemistry_CH4_ChemCH43DPrimArrhenius_hh
#define COOLFluiD_Physics_Chemistry_CH4_ChemCH43DPrimArrhenius_hh

//////////////////////////////////////////////////////////////////////////////

#include "ChemCH43DSourceVarSet.hh"
#include "ChemCH43DPrim.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Source term definied relative to the ChemCH43DPrim varset
   *
   * @author Tiago Quintino
   */
class ChemCH43DPrimArrhenius : public ChemCH43DSourceVarSet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  ChemCH43DPrimArrhenius(const std::string& name);

  /**
   * Default destructor
   */
  ~ChemCH43DPrimArrhenius();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

 /**
  * Gets the SourceFlux Coeficients
  */
  virtual void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef);

 /**
  * Gets the SourceFlux Coeficients
  */
  virtual void getIndepSourceCoefs(const Framework::State& state, RealVector& coef);

  /**
   * Checks if the source term has a part independent from the solution
   */
  bool hasIndepCoef() const
  {
    return true;
  }
  
  /**
   * Checks if the source term has a linear part
   */
  bool hasLinearCoef() const
  {
    return false;
  } 
  
private: // member data

  /// the Arrhenius constant
  CFreal _A;

  /// the exponent coefficient in the Arrhenius expression, divided by R
  CFreal _EaR;

  /// the exponent coefficient in the Arrhenius expression
  CFreal _Ea;

  /// the exponent of the CH4 concentration
  CFreal _m;

  /// the exponent of the O2 concentration
  CFreal _n;

  /// the parameters multiplying the Arrhenius terms
  RealVector _multP;

}; // end of class ChemCH43DPrimArrhenius

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Chemistry_CH4_ChemCH43DPrimArrhenius_hh
