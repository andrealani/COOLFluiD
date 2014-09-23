#ifndef COOLFluiD_Numerics_LUSGSMethod_LUSGSCrankNich_hh
#define COOLFluiD_Numerics_LUSGSMethod_LUSGSCrankNich_hh

//////////////////////////////////////////////////////////////////////////////

#include "LUSGSMethod/LUSGSIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a ConvergenceMethod that implements the (nonlinear)
 * LU-SGS method applied to a Crank-Nicholson time marching scheme.
 *
 * @author Kris Van den Abeele
 */
class LUSGSCrankNich : public LUSGSIterator {
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
  explicit LUSGSCrankNich(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LUSGSCrankNich();

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

protected: // abstract interface implementations

  /**
   * Take one timestep.
   * @see ConvergenceMethod::takeStep()
   */
  virtual void takeStepImpl();

  /// Perform the prepare phase before any iteration
  virtual void prepare ();

protected: // member data

  /// The command that backs up the past rhs
  Common::SelfRegistPtr<LUSGSIteratorCom> m_backupPastRhs;

  /// The command that adds the past rhs to the current rhs
  Common::SelfRegistPtr<LUSGSIteratorCom> m_addPastRhs;

  ///The string for configuration of m_backupPastRhs command
  std::string m_backupPastRhsStr;

  ///The string for configuration of m_addPastRhs command
  std::string m_addPastRhsStr;

}; // class LUSGSCrankNich

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_hh
