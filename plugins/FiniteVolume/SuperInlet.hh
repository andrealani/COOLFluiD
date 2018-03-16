#ifndef COOLFluiD_Numerics_FiniteVolume_SuperInlet_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperInlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic inlet command
   *
   * @author Andrea Lani
   *
   */
class SuperInlet : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SuperInlet(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SuperInlet();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

protected: // data

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> _varSet;

  /// storage for the temporary boundary point coordinates
  RealVector _bCoord;
  
  /// RealVector holding [x,y,z] and iteration number
  RealVector _xyzIter;
  
  /// Temporary state for the dimensional values
  Framework::State* _dimState;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// flag to know if the inputted values are adimensional
  bool _inputAdimensionalValues;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> _inputToUpdateVar;
  
  /// State holding the input variables
  Framework::State* _input;

  /// a string to hold the name of the input variables
  std::string _inputVarStr;

  /// a string to hold the name of the update variables
  std::string _updateVarStr;
  
  /// IDs of the variables that will be changed interactively
  std::vector<CFuint> _interVarIDs;
  
  /// Factor to multiply the selected InteractiveVarIDs (should be < 1)
  CFreal _interFactor;
  
}; // end of class SuperInlet

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperInlet_hh
