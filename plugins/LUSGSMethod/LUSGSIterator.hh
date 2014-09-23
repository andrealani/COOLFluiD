#ifndef COOLFluiD_Numerics_LUSGSMethod_LUSGSIterator_hh
#define COOLFluiD_Numerics_LUSGSMethod_LUSGSIterator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethod.hh"

#include "LUSGSMethod/LUSGSIteratorData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a ConvergenceMethod that implements the (nonlinear)
 * LU-SGS method.
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class LUSGSIterator : public Framework::ConvergenceMethod {
public:

  /**
   * Defines the Config Option's of this class.
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   *
   * @param name missing documentation
   */
  explicit LUSGSIterator(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LUSGSIterator();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Defined the strategy list of this Method
  /// @todo Remove this function and configure strategies using the configureStrategy template function
  ///       provided by MethodData
  std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategyList() const;

protected: // helper functions

  /**
   * Gets the Data aggregator of this method.
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Gets the Data aggregator of this method.
   * @return SafePtr to the ConvergenceMethodData
   */
  virtual Common::SafePtr<Framework::ConvergenceMethodData> getConvergenceMethodData();

protected: // abstract interface implementations

  /**
   * Take one timestep.
   * @see ConvergenceMethod::takeStep()
   */
  virtual void takeStepImpl();

  /**
   * Sets up the data for the method commands to be applied.
   * @see  Method::setMethod()
   */
  virtual void setMethodImpl();

  /**
   *  Un set up the data of the method.
   * @see Method::unsetMethod()
   */
  virtual void unsetMethodImpl();

  /// Perform the prepare phase before any iteration
  virtual void prepare ();

protected: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<LUSGSIteratorCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<LUSGSIteratorCom> m_unSetup;

  ///The Prepare command to use
  Common::SelfRegistPtr<LUSGSIteratorCom> m_prepare;

  /// The Intermediate command to use between computing the
  /// space residual and the time residual
  Common::SelfRegistPtr<LUSGSIteratorCom> m_intermediate;

  ///The Initial command to use at beginning of each iteration
  Common::SelfRegistPtr<LUSGSIteratorCom> m_init;

  ///The UpdateSolution command to use
  Common::SelfRegistPtr<LUSGSIteratorCom> m_updateSol;

  ///The AleUpdate command to use
  Common::SelfRegistPtr<LUSGSIteratorCom> m_aleUpdate;

  ///The UpdateStatesSetIndex command to use
  Common::SelfRegistPtr<LUSGSIteratorCom> m_updateStatesSetIndex;

  ///The LUFactorization command to use
  Common::SelfRegistPtr<LUSGSIteratorCom> m_luFactorization;

  ///The TriangularSystemSolver command to use
  Common::SelfRegistPtr<LUSGSIteratorCom> m_computeStatesSetUpdate;

  ///The command that computes the diagonal block Jacobians by perturbation of the states
  Common::SelfRegistPtr<LUSGSIteratorCom> m_diagBlockJacobComputer;

  ///The string for configuration of m_setup command
  std::string m_setupStr;

  ///The string for configuration of m_unSetup command
  std::string m_unSetupStr;

  ///The string for configuration of m_prepare command
  std::string m_prepareStr;

  ///The string for configuration of m_intermidiate command
  std::string m_intermediateStr;

  ///The string for configuration of m_init command
  std::string m_initStr;

  ///The string for configuration of m_updateSol command
  std::string m_updateSolStr;

  ///The string for configuration of m_aleUpdate command
  std::string m_aleUpdateStr;

  ///The string for configuration of m_updateSol command
  std::string m_updateStatesSetIndexStr;

  ///The string for configuration of m_aleUpdate command
  std::string m_luFactorizationStr;

  ///The string for configuration of m_aleUpdate command
  std::string m_computeStatesSetUpdateStr;

  ///The string for configuration of m_diagBlockJacobComputer command
  std::string m_diagBlockJacobComputerStr;

  ///The data to share between LUSGSMethodMethod commands
  Common::SharedPtr<LUSGSIteratorData> m_data;

}; // class LUSGSIterator

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_hh
