// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_FwdEuler_hh
#define COOLFluiD_Numerics_ForwardEuler_FwdEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethod.hh"
#include "FwdEulerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {  class NumericalCommand;  }

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the Explicit
/// time stepping algorithm of first order.
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class  ForwardEuler_API FwdEuler : public Framework::ConvergenceMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit FwdEuler(const std::string& name);

  /// Default destructor
  ~FwdEuler();

  /// Configures the method, by allocating the it's dynamic members.
  /// @param args arguments from where to read the configuration
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
  Common::SelfRegistPtr<FwdEulerCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<FwdEulerCom> m_unSetup;

  ///The Prepare command to use
  Common::SelfRegistPtr<FwdEulerCom> m_prepare;

  ///The Intermediate command to use
  Common::SelfRegistPtr<FwdEulerCom> m_intermediate;

  ///The StdUpdateSolution command to use
  Common::SelfRegistPtr<FwdEulerCom> m_updateSol;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///The Prepare string for configuration
  std::string m_prepareStr;

  ///The Intermediate string for configuration
  std::string m_intermediateStr;

  ///The StdUpdateSolution for configuration
  std::string m_updateSolStr;

  ///The data to share between ForwardEulerMethod commands
  Common::SharedPtr<FwdEulerData> m_data;

}; // class FwdEuler

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ForwardEuler_hh
