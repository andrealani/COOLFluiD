// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_BackwardEuler_BwdEuler_hh
#define COOLFluiD_Numerics_BackwardEuler_BwdEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethod.hh"
#include "BwdEulerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the Explicit
/// time stepping algorithm of first order.
/// @author Andrea Lani
class BwdEuler : public Framework::ConvergenceMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit BwdEuler(const std::string& name);

  /// Default destructor
  ~BwdEuler();

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

  ///The data to share between BackwardEulerMethod commands
  Common::SharedPtr<BwdEulerData> m_data;

  ///The Setup command to use
  Common::SelfRegistPtr<BwdEulerCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<BwdEulerCom> m_unSetup;

  ///The updateSolution command to use
  Common::SelfRegistPtr<BwdEulerCom> m_updateSol;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///The updateSolution for configuration
  std::string m_updateSolStr;

}; // class BwdEuler

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_BackwardEuler_hh
