// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_CFEnv_hh
#define COOLFluiD_Environment_CFEnv_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Common/COOLFluiD.hh"
#include "Common/SafePtr.hh"
#include "Common/SetupObject.hh"
#include "Common/SharedPtr.hh"
#include "Config/ConfigObject.hh"

#include "Environment/Environment.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common
  {
    class EventHandler;
    class FactoryRegistry;
  }

  namespace Environment {
    
    class ModuleRegistry;
    class CFEnvVars;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a singleton object where
/// which is used to initialize the COOLFluiD runtime environment.
/// @author Tiago Quintino
class Environment_API  CFEnv :
      public Config::ConfigObject,
      public Common::SetupObject,
      public Common::NonCopyable<CFEnv> {
  
 public: // methods

  /// @return the instance of this singleton
  static CFEnv& getInstance();
  
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configures this Simulator.
  /// Sets up the data for this Object.
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup the environment
  /// @pre called after configure
  virtual void setup();
  /// Unsetup the object
  /// @pre called after setup
  virtual void unsetup();

  /// Initializes the COOLFluiD runtime enviroment.
  /// @pre Must be called prior to any COOLFluiD runtime function,
  ///      only module registration procedures are allowed beforehand.
  /// @TODO This is still broken for mpi as it doesn't allow modification
  ///  of argc & argv
  void initiate(int argc, char** argv);

  /// Closes the COOLFluiD runtime environment.
  /// @post Must not call any COOLFluiD runtime functions after,
  ///       only destruction procedures ar allowed afterwards.
  void terminate();

  /// Gets the ModuleRegistry
  /// @note Does not need to be initialized before
  Common::SafePtr<Environment::ModuleRegistry> getModuleRegistry();

  /// Gets the FactoryRegistry
  /// @note Does not need to be initialized before
  Common::SafePtr<Common::FactoryRegistry> getFactoryRegistry();
  
  /// Gets the EventHandler of the COOLFluiD runtime environment
  /// @note Does not need to be initialized before
  Common::SafePtr<Common::EventHandler> getEventHandler();
  
  /// Gets the CFEnvVars
  Common::SafePtr<Environment::CFEnvVars> getVars();

  /// Gets the local CPU rank
  /// @post return 0 if parallel environment hasn't been initiated yet
  /// @return CPU rank
  CFuint getCPURank();

  /// Initializes the loggers
  void initLoggers();

  /// Calls initiate() on all registered modules.
  /// Mind that some modules might already have been initiated.
  /// It is up to the modules to track if they have or not been initiated.
  /// @see ModuleRegisterBase
  void initiateModules();

  /// Calls terminate() on all registered modules
  /// Mind that some modules might already have been terminated.
  /// It is up to the modules to track if they have or not been terminated.
  /// @see ModuleRegisterBase
  void terminateModules();

  /// Return the version string of this build
  std::string getVersionString () const;
  /// Return the subversion version string of this build
  std::string getSvnVersion () const;
  /// Return the COOLFluiD version string
  std::string getCFVersion () const;
  /// Return the COOLFluiD Kernel version string
  std::string getKernelVersion () const;
  /// Return the COOLFluiD build type
  std::string getBuildType () const;
  /// Return the CMake version
  std::string getBuildSystem () const;
  /// Return the build processor
  std::string getBuildProcessor () const;
  /// OS short name. Examples: "Linux" or "Windows"
  /// @return string with short OS name
  std::string getSystemName() const;
  /// OS short name. Examples: "Linux-2.6.23" or "Windows 5.1"
  /// @return string with long OS name
  std::string getLongSystemName() const;
  /// OS version. Examples: "2.6.23" or "5.1"
  /// @return string with OS version
  std::string getSystemVersion() const;
  /// OS bits. Examples: "32" or "64"
  /// @post should be equal to 8 * size_of(void*) but it given by the build system
  /// @return string with OS addressing size
  std::string getSystemBits() const;

private: // methods

  /// Default contructor
  CFEnv();
  /// Default destructor
  ~CFEnv();

private: // data

  /// the EventHandler object is only held by the CFEnv singleton object
  Common::SharedPtr<Common::EventHandler> m_eventHandler;

  /// the ModuleRegistry singleton object is only held by the CFEnv singleton object
  Environment::ModuleRegistry* m_moduleRegistry;

  /// the FactoryRegistry singleton object is only held by the CFEnv singleton object
  Common::FactoryRegistry* m_factoryRegistry;
  
  /// @brief Static environment variables
  /// pointer to a struct of variables that always exist
  CFEnvVars * m_env_vars;

}; // end of class CFEnv

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_CFEnv_hh
