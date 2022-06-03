#ifndef COOLFluiD_FluxReconstructionMethod_InitStateAddVar_hh
#define COOLFluiD_FluxReconstructionMethod_InitStateAddVar_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/StdInitState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * A initalizing solution command for the flux reconstruction method with added 
 * variables
 *
 * @author Ray Vandenhoeck
 *
 */
class InitStateAddVar : public StdInitState {
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
  ~InitStateAddVar();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

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
  RealVector m_tmpFun;
  
  /// array for storing temporary variables to pass to the 
  /// function parser 
  RealVector m_tmpVars;
  
  /// a vector of string to hold the functions
  std::vector<std::string> m_initFunctions;
  
  /// a vector of string to hold the functions
  std::vector<std::string> m_initVars;
  
  // the VectorialFunction to use
  Framework::VectorialFunction m_vInitFunction;

}; // class InitStateAddVar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_InitStateAddVar_hh

