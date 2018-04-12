#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_AnalyticalSourceTerm_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_AnalyticalSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionMethod/StdSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * A command for adding an analytically defined source term
 *
 * @author Ray Vandenhoeck
 */
class AnalyticalSourceTerm : public StdSourceTerm {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit AnalyticalSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~AnalyticalSourceTerm();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * get data required for source term computation
   */
  void getSourceTermData();

  /**
   * add the source term
   */
  void addSourceTerm();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // data

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the variables
  std::vector<std::string> m_vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// the input variables needed to evaluate the function
  RealVector m_vFunctionInputVars;

  /// the source term for one state
  RealVector m_srcTerm;

  /// dimensionality
  CFuint m_dim;

}; // class AnalyticalSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_AnalyticalSourceTerm_hh

