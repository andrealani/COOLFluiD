#ifndef COOLFluiD_Numerics_RKRD_RKRD_hh
#define COOLFluiD_Numerics_RKRD_RKRD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethod.hh"
#include "RKRD/RKRDData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a RungeKutta iteration to use with RKRD FSM method
/// @author Mario Ricchiuto
/// @author Tiago Quintino
class RKRD : public Framework::ConvergenceMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  /// @param name missing documentation
  explicit RKRD(const std::string& name);

  /// Default destructor
  virtual ~RKRD();

  /// Configures the method, by allocating the it's dynamic members.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

protected: // helper functions

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the ConvergenceMethodData
  virtual Common::SafePtr<Framework::ConvergenceMethodData> getConvergenceMethodData();

protected: // abstract interface implementations

  /// Take one timestep
  /// @see ConvergenceMethod::takeStep()
  virtual void takeStepImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

protected: // member data

  /// The Setup command to use
  Common::SelfRegistPtr<RKRDCom> m_setup;

  /// The UnSetup command to use
  Common::SelfRegistPtr<RKRDCom> m_unSetup;

  /// The Backup command to use
  Common::SelfRegistPtr<RKRDCom> m_backup;

  /// The Update command to use
  Common::SelfRegistPtr<RKRDCom> m_update;

  /// The Shift command to use
  Common::SelfRegistPtr<RKRDCom> m_shift;

  /// string for configuration of setup  command
  std::string m_setupStr;

  /// string for configuration of unsetup command
  std::string m_unSetupStr;

  /// string for configuration of backup command
  std::string m_backupStr;

  /// string for configuration of update command
  std::string m_updateStr;

  /// string for configuration of shift command
  std::string m_shiftStr;

  /// The data to share between RKRDMethod commands
  Common::SharedPtr<RKRDData> m_data;

}; // class RKRD

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RKRD_RKRD_hh
