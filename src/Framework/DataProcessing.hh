// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataProcessing_hh
#define COOLFluiD_Framework_DataProcessing_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingMethod.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/DataProcessingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a DataProcessing.
/// @author Thomas Wuilbaut
/// @author Tiago Quintino
class Framework_API DataProcessing : public DataProcessingMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  DataProcessing(const std::string& name);

  /// Default destructor
  virtual ~DataProcessing();

  /// Configure
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName() {  return "DataProcessing";  }

protected: // abstract interface implementations

  /// Execute the data processing on the trs
  /// This is the abstract function that the concrete methods must implement.
  virtual void processDataImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

protected: // helper functions

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Clear the commands
  void clearComs();

private: // member data
  
  ///The data to share between DataProcessing commands
  Common::SharedPtr<DataProcessingData> m_data;
  
  ///The command to use for initializing the solution.
  std::vector<Common::SelfRegistPtr<DataProcessingCom> > m_dataprocessing;
  
  ///The data processing command types
  std::vector<std::string> m_dataprocessTypeStr;

  ///The data processing command names for configuration
  std::vector<std::string> m_dataprocessNameStr;
  
  /// flag telling to skip the first iteration
  bool m_skipFirstIteration;
  
  /// flag telling to run at setup time
  bool m_runAtSetup;

  /// flag telling to run at setup time and after use process rate
  bool m_runAtSetupAndAfter; 
  
}; // class DataProcessing

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataProcessing_hh
