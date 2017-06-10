// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Simulator_hh
#define COOLFluiD_Simulator_hh

//////////////////////////////////////////////////////////////////////////////


#include <boost/filesystem/path.hpp>

#include "Common/NonCopyable.hh"
#include "Common/SelfRegistPtr.hh"
#include "Common/OwnedObject.hh"
#include "Common/DynamicObject.hh"

#include "Config/ConfigObject.hh"

#ifndef CF_HAVE_ALLSTATIC
#include "Environment/ModuleLoader.hh"
#endif

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

    class SubSystem;

//////////////////////////////////////////////////////////////////////////////

class Framework_API SimulatorPaths : public Common::OwnedObject,
                                     public Config::ConfigObject,
                                     public Common::NonCopyable<SimulatorPaths>
{
public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor with given argument.
  /// @param startFile name of the file to read the parameters from.
  explicit SimulatorPaths (const std::string& name);

  /// Configures this Simulator.
  /// Sets up the data for this Object.
  virtual void configure ( Config::ConfigArgs& args );

private: // data

  /// dir path string for configuration
  std::string m_base_dir;
  /// strings for configuration of the modules directories
  std::vector< std::string > m_modules_dirs;
  /// dir path string for configuration
  std::string m_working_dir;
  /// dir path string for configuration
  std::string m_results_dir;
  /// the url for file repository
  std::string m_repository_url;

};

//////////////////////////////////////////////////////////////////////////////

/// This class represents a simulator.
/// A Simulator controls a Simulation.
/// @see SubSystem
class Framework_API Simulator :
    public Common::OwnedObject,
    public Config::ConfigObject,
    public Common::NonCopyable<Simulator> {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor with given argument.
  /// @param name name of the Simualtor object
  explicit Simulator(const std::string& name);

  /// Destructor.
  ~Simulator();

  /// Sets a single CFcase file for simulation
  void openCaseFile (const std::string& sCFcaseFile);

  /// Adds the ActionListener's of this EventListener to the EventHandler
  void registActionListeners();

  /// Access the existing sub systems types in the Simulator
  std::vector< std::string > getSubSystemTypes() const;

  /// Access the existing sub systems names in the Simulator
  std::vector< std::string > getSubSystemNames() const;
  
  /// Tell if the given rank belongs to the given subsystem 
  bool isSubSystemRank(const CFuint rank, const std::string& subSystemName) const;
  
protected: // functions

  /// Processes the argument list from the configuration arguments
  /// Overriden to propagate the file configuration down the object hierarchy
  virtual void processOptions ( Config::ConfigArgs& args );

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_CONFIGSUBSYSTEM" Event
  /// @param eConfig the event which provoked this action
  /// @return an Event with a message in its body replying to the action request
  Common::Signal::return_t configSubSystem(Common::Signal::arg_t eConfig);

  /// Build the SubSystem according to the name specified in the Event
  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_BUILDSUBSYSTEM" Event
  /// @param eBuild the event which provoked this action
  /// @return an Event with a message in its body replying to the action request
  Common::Signal::return_t buildSubSystem(Common::Signal::arg_t eBuild);

  /// Destroys the SubSystem according to the name specified in the Event
  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_DESTROYSUBSYSTEM" Event
  /// @param eDestroy the event which provoked this action
  /// @return an Event with a message in its body replying to the action request
  Common::Signal::return_t destroySubSystem(Common::Signal::arg_t eDestroy);

  /// Configures this Simulator.
  /// Sets up the data for this Object.
  virtual void configure ( Config::ConfigArgs& args );
    
 private: // data

#ifndef CF_HAVE_ALLSTATIC
  /// Module loader
  Environment::ModuleLoader m_moduleLoader;
#endif

  /// the names of the SubSytem's to create in this simulation
  std::vector<std::string> m_subSystemNames;
  /// the types of the SubSytem's to create in this simulation
  std::vector<std::string> m_subSystemTypes;
  /// the current SubSystem being handled in this processor
  Common::SelfRegistPtr<SubSystem> m_subSys;
  /// Simulation configuration arguments
  Config::ConfigArgs m_sim_args;
  
  /// map the subsystem names to the corresponding ranks
  std::map<std::string, std::string> m_mapSS2Ranks;
  
  /// Paths object
  Common::SharedPtr<SimulatorPaths> m_paths;
  
  /// ranks associated to each namespace in the form
  /// START0:END0 START1:END1 etc.
  std::vector<std::string> m_ranksString;
  
}; // class Simulator

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Simulator_hh
