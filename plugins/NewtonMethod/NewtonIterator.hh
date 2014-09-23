// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_NewtonIterator_hh
#define COOLFluiD_Numerics_NewtonMethod_NewtonIterator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethod.hh"

#include "NewtonMethod/NewtonIteratorData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the Explicit
/// time stepping algorithm of first order.
/// @author Tiago Quintino
/// @author Andrea Lani
class NewtonIterator : public Framework::ConvergenceMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit NewtonIterator(const std::string& name);

  /// Default destructor
  virtual ~NewtonIterator();

  /// Configures the method, by allocating the it's dynamic members.
  /// @param args arguments from where to read the configuration
  virtual void configure ( Config::ConfigArgs& args );
 
  /// Sets the SpaceMethod
  virtual void setCollaborator(Framework::MultiMethodHandle<Framework::SpaceMethod> sm) {m_data->setSpaceMethod(sm);}
  
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

  /// Perform the prepare phase before any iteration
  virtual void prepare ();

protected: // member data

///The data to share between NewtonMethodMethod commands
  Common::SharedPtr<NewtonIteratorData> m_data;

  ///The Setup command to use
  Common::SelfRegistPtr<NewtonIteratorCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<NewtonIteratorCom> m_unSetup;

  ///The Prepare command to use
  Common::SelfRegistPtr<NewtonIteratorCom> m_prepare;

  /// The Intermediate command to use between computing the
  /// space residual and the time residual
  Common::SelfRegistPtr<NewtonIteratorCom> m_intermediate;

  ///The Initial command to use at beginning of each iteration
  Common::SelfRegistPtr<NewtonIteratorCom> m_init;

  ///The UpdateSolution command to use
  Common::SelfRegistPtr<NewtonIteratorCom> m_updateSol;

  ///The AleUpdate command to use
  Common::SelfRegistPtr<NewtonIteratorCom> m_aleUpdate;

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

}; // class NewtonIterator

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_hh
