#ifndef COOLFluiD_Numerics_FiniteVolume_InitStateAddVar_hh
#define COOLFluiD_Numerics_FiniteVolume_InitStateAddVar_hh

//////////////////////////////////////////////////////////////////////////////

#include "InitState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a initalizing solution command
   *
   * @author Andrea Lani
   *
   */
class InitStateAddVar : public InitState {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor.
   */
  explicit InitStateAddVar(const std::string& name);
  
  /**
   * Destructor.
   */
  virtual ~InitStateAddVar();
  
  /**
   * Set up private data
   */
  virtual void setup();
  
protected:

  /**
   * Execute Processing actions
   */
  void executeOnTrs();
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // data
  
  /// array for storing temporary function results 
  RealVector _tmpFun;
  
  /// array for storing temporary variables to pass to the 
  /// function parser 
  RealVector _tmpVars;
  
  /// a vector of string to hold the functions
  std::vector<std::string> _initFunctions;
  
  /// a vector of string to hold the functions
  std::vector<std::string> _initVars;
  
  // the VectorialFunction to use
  Framework::VectorialFunction _vInitFunction;
  
}; // class InitStateAddVar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_InitStateAddVar_hh

