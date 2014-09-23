#ifndef COOLFluiD_Physics_Chemistry_CH4_ChemCH43DPrimDiffusive_hh
#define COOLFluiD_Physics_Chemistry_CH4_ChemCH43DPrimDiffusive_hh

//////////////////////////////////////////////////////////////////////////////

#include "ChemCH43DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Chemistry physical model 2D for primitive
   * variables
   *
   * @author Tiago Quintino
   */
class ChemCH43DPrimDiffusive : public ChemCH43DDiffusiveVarSet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  ChemCH43DPrimDiffusive(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~ChemCH43DPrimDiffusive();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the diffusive parameter
   */
  std::vector<CFreal>& getDiffusiveCoefs()
  {
    return _vecD;
  }



private:

  /// vector to store the configuration of the diffusive parameters
  std::vector<CFreal> _vecD;

}; // end of class ChemCH43DPrimDiffusive

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Chemistry_CH4_ChemCH43DPrimDiffusive_hh
