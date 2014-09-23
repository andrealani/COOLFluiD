#ifndef COOLFluiD_Numerics_MatrixStabilityMethod_hh
#define COOLFluiD_Numerics_MatrixStabilityMethod_hh

//////////////////////////////////////////////////////////////////////////////



#include "Framework/ConvergenceMethod.hh"
#include "MatrixStabilityMethodWriter/MatrixStabilityMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a ConvergenceMethod that creates the matrix corresponding to
 * a spatial method.
 * @warning this is not actually a convergence method!
 *
 * @author Kris Van den Abeele
 *
 */
class MatrixStabilityMethod : public Framework::ConvergenceMethod {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   *
   * @param name missing documentation
   */
  explicit MatrixStabilityMethod(const std::string& name);

  /**
   * Default destructor
   */
  ~MatrixStabilityMethod();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // helper functions

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the ConvergenceMethodData
   */
  virtual Common::SafePtr<Framework::ConvergenceMethodData> getConvergenceMethodData();

protected: // abstract interface implementations

  /**
   * Take one timestep
   * @see ConvergenceMethod::takeStep()
   */
  virtual void takeStepImpl();

  /**
   * Sets up the data for the method commands to be applied.
   * @see Method::unsetMethod()
   */
  virtual void unsetMethodImpl();

  /**
   * UnSets the data of the method.
   * @see Method::setMethod()
   */
  virtual void setMethodImpl();

protected: // member data

  /// The Setup command to use
  Common::SelfRegistPtr<MatrixStabilityMethodCom> m_setup;

  /// The UnSetup command to use
  Common::SelfRegistPtr<MatrixStabilityMethodCom> m_unSetup;

  /// The SetStates command to use
  Common::SelfRegistPtr<MatrixStabilityMethodCom> m_setStates;

  /// The AddMatrixColumnToFile command to use
  Common::SelfRegistPtr<MatrixStabilityMethodCom> m_addMatrixColumnToFile;

  /// The Setup string for configuration
  std::string m_setupStr;

  /// The UnSetup string for configuration
  std::string m_unSetupStr;

  /// The SetStates string for configuration
  std::string m_setStatesStr;

  /// The AddMatrixColumnToFile string for configuration
  std::string m_addMatrixColumnToFileStr;

  /// The data to share between MatrixStabilityMethod commands
  Common::SharedPtr<MatrixStabilityMethodData> m_data;

}; // class MatrixStabilityMethod

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MatrixStabilityMethod_hh
