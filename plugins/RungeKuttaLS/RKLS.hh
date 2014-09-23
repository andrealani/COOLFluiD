// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKuttaLS_RKLS_hh
#define COOLFluiD_Numerics_RungeKuttaLS_RKLS_hh

//////////////////////////////////////////////////////////////////////////////



#include "Framework/ConvergenceMethod.hh"
#include "RungeKuttaLS/RKLSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace RungeKuttaLS {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a ConvergenceMethod that implements the explicit low storage
 * Runge-Kutta algorithms
 *
 * @author Kris Van den Abeele
 *
 */
class RKLS : public Framework::ConvergenceMethod {
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
  explicit RKLS(const std::string& name);

  /**
   * Default destructor
   */
  ~RKLS();

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


  ///The Setup command to use
  Common::SelfRegistPtr<RKLSCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<RKLSCom> m_unSetup;

  ///The backupSol command to use
  Common::SelfRegistPtr<RKLSCom> m_backupSol;

  ///The predictor step command to use
  Common::SelfRegistPtr<RKLSCom> m_rungeKuttaStep;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///The backup solution command string for configuration
  std::string m_backupSolStr;

  ///The corrector step command string for configuration
  std::string m_rungeKuttaStepStr;

  ///The data to share between RungeKuttaMethod commands
  Common::SharedPtr<RKLSData> m_data;

}; // class RKLS

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKuttaLS

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RungeKuttaLS_RKLS_hh
