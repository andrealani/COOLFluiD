#ifndef COOLFluiD_Numerics_LUSGSMethod_LUSGSIteratorComputDiagJacob_hh
#define COOLFluiD_Numerics_LUSGSMethod_LUSGSIteratorComputDiagJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "LUSGSMethod/LUSGSIterator.hh"
#include "LUSGSMethod/LUSGSIteratorData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a ConvergenceMethod that implements the (nonlinear)
 * LU-SGS method. The block diagonal Jacobian matrices are computed by the convergence method
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class LUSGSIteratorComputDiagJacob : public LUSGSIterator {
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
  explicit LUSGSIteratorComputDiagJacob(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LUSGSIteratorComputDiagJacob();

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

protected: // member data

}; // class LUSGSIteratorComputDiagJacob

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_hh
