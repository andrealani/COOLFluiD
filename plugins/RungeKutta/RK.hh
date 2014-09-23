// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKutta_RK_hh
#define COOLFluiD_Numerics_RungeKutta_RK_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethod.hh"
#include "RungeKutta/RKData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace RungeKutta {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the Explicit Time
/// stepping Runge-Kutta algorithm.
/// @author Thomas Wuilbaut
class RK : public Framework::ConvergenceMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  /// @param name missing documentation
  explicit RK(const std::string& name);

  /// Default destructor
  virtual ~RK();

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

  ///The Setup command to use
  Common::SelfRegistPtr<RKCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<RKCom> m_unSetup;

  ///The backupSol command to use
  Common::SelfRegistPtr<RKCom> m_backupSol;

  ///The predictor step command to use
  Common::SelfRegistPtr<RKCom> m_rungeKuttaStep;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///The backup solution command string for configuration
  std::string m_backupSolStr;

  ///The corrector step command string for configuration
  std::string m_rungeKuttaStepStr;

  ///The data to share between RungeKuttaMethod commands
  Common::SharedPtr<RKData> m_data;

}; // class RK

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKutta_RK_hh
