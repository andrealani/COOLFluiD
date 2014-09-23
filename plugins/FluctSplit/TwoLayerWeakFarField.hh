#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerWeakFarField_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerWeakFarField_hh

//////////////////////////////////////////////////////////////////////////////

#include "TwoLayerWeakBC2D.hh"
#include "Framework/VectorialFunction.hh"
#include "Framework/VarSetTransformer.hh"
#include "Common/SelfRegistPtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak farfield bc for 2Layer SpaceTime Schemes
 *
 * @author Thomas Wuilbaut
 *
 */
class TwoLayerWeakFarField : public TwoLayerWeakBC2D {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::VarSetTransformer VectorTransformer;

  /**
   * Constructor.
   */
  TwoLayerWeakFarField(const std::string& name);

  /**
   * Default destructor
   */
  ~TwoLayerWeakFarField();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up the member data
   */
  void setup();

  /**
   * Unset up the member data
   */
  void unsetup();

 protected:

  /**
   * Set the state vector in the ghost State's
   */
  void setGhostState(const Framework::State& state,
                     Framework::State& gstate);

 private:

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<VectorTransformer> _inputToUpdateVar;

  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  /// RealVector holding the BC variables
  RealVector _variables;

  /// State holding the input (ghost) variables
  Framework::State* _input;

  /// State holding for the adimesnionalization
  Framework::State* _dimState;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// a string to hold the name of the input variables
  std::string _inputVarStr;

  /// a string to hold the name of the update variables
  std::string _updateVarStr;

 }; // end of class TwoLayerWeakFarField

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerWeakFarField_hh
