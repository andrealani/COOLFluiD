// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_MarcoTest_MarcoTestMethod_hh
#define COOLFluiD_Numerics_MarcoTest_MarcoTestMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethod.hh"
#include "MarcoTest/MarcoTestMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {  class NumericalCommand;  }

    namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the Explicit
/// time stepping algorithm of first order.
/// @author Andrea Lani
/// @author Marco Panesi
class  MarcoTest_API MarcoTestMethod : public Framework::ConvergenceMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit MarcoTestMethod(const std::string& name);

  /// Default destructor
  ~MarcoTestMethod();

  /// Configures the method, by allocating the it's dynamic members.
  /// @param args arguments from where to read the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets a vector with all the NumericalStrategy's this method will use.
  /// @return vector with the strategy pointers.
  virtual std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategyList () const;
  
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
  
  ///The data to share between MarcoTestMethod commands
  Common::SharedPtr<MarcoTestMethodData> m_data;
  
  /// Command to setup the data structure
  Common::SelfRegistPtr<MarcoTestMethodCom> m_setup;
  
  /// Command to implement the whole algorithm
  Common::SelfRegistPtr<MarcoTestMethodCom> m_algo;
  
  /// Command to unsetup the data structure
  Common::SelfRegistPtr<MarcoTestMethodCom> m_unsetup;
  
  /// The setup string for configuration
  std::string m_setupStr;
  
  /// The algo string for configuration
  std::string m_algoStr;
  
  /// The unsetup string for configuration
  std::string m_unsetupStr;
  
  
}; // class MarcoTestMethod

//////////////////////////////////////////////////////////////////////////////

    } // namespace MarcoTest



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MarcoTest_hh
