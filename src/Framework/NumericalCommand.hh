// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NumericalCommand_hh
#define COOLFluiD_Framework_NumericalCommand_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/OwnedObject.hh"
#include "Common/SetupObject.hh"
#include "Config/ConfigObject.hh"

#include "Environment/ConcreteProvider.hh"
#include "Framework/TopologicalRegionSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class CommandGroup;
    class BaseDataSocketSource;
    class BaseDataSocketSink;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand.
/// It is the abstract class definning the interface to the derived NumericalCommands
/// NumericalCommand allow a given Method to be costumized through the setup
/// of the commands without subclassing. A given Method might use a several ways to compute
/// something. Instead of multiclassing, the user chooses a suitable Command and configures
/// the method with it.
/// @author Tiago Quintino
/// @author Andrea Lani
class Framework_API NumericalCommand :
    public Common::OwnedObject,
    public Common::SetupObject,
    public Config::ConfigObject,
    public Common::NonCopyable<NumericalCommand> {

public: // typedefs

  typedef Environment::ConcreteProvider<NumericalCommand> PROVIDER;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor
  explicit NumericalCommand(const std::string& name);

  /// Default destructor pure virtual
  virtual ~NumericalCommand() = 0;

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unset up private data and data of the aggregated classes
  /// in this command after the processing phase
  virtual void unsetup();

  /// Execute the command
  /// @pre this method has to be overridden by concrete NumericalCommand's
  ///      that don't need TRS
  virtual void execute();

  /// Set the TRS list
  void setTrsList(const std::vector< Common::SafePtr<TopologicalRegionSet> >& trsList) { m_trsList = trsList; }

  /// Gets the names of the TopologicalRegionSet's that this NumericalCommand
  /// applies to. These names are configured at startup time.
  /// @post return != CFNULL
  const std::vector<std::string>& getTrsNames() const { return m_trsNames; }

  /// Set CommandGroup
  void setCommandGroup(Common::SafePtr<CommandGroup> commandGroup);

  /// Get CommandGroup
  Common::SafePtr<CommandGroup> getCommandGroup() const {  return m_group;  }

  /// Get CommandGroupName
  std::string getCommandGroupName() const {  return m_groupName;  }

  /// Get CommandGroupName
  void setCommandGroupName(std::string groupName) { m_groupName = groupName; }

  /// Gets the Class name
  static std::string getClassName() { return "NumericalCommand"; }

  /// Configure the nested sockets in this NumericalCommand
  void configureNestedSockets ( Config::ConfigArgs& args );

  /// Allocates all sockets in this NumericalCommand
  void allocateCommandSockets();

  /// Deallocates all sockets in this NumericalCommand
  void deallocateCommandSockets();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<BaseDataSocketSink> > needsSockets();
  
protected: // functions

  /// Execute on Trs
  /// @pre this method has to be overridden by all the concrete
  ///      NumericalCommand's that compute on TRSs
  virtual void executeOnTrs();

  /// Sets the current TopologicalRegionSet to work on
  /// @param iTrs supplied TRS to apply this command.
  void setCurrentTrsID(const CFuint iTrs);

  /// Get the current trs to work on
  Common::SafePtr<TopologicalRegionSet> getCurrentTRS() const;

  /// Get the ID of the current trs to work on
  CFuint getCurrentTrsID() const { return m_iTrs; }

  /// Get the process rate
  CFuint getProcessRate() const { return m_processRate; }
  
  /// Get the TRS list
  std::vector< Common::SafePtr<TopologicalRegionSet> >& getTrsList();
  
  /// Get the TRS name
  const std::string getTrsName(const CFuint iTrs);

private: // data

  /// ID of the current TRS to work on
  CFuint                                    m_iTrs;

  /// list of the TRSs on which this command applies
  std::vector< Common::SafePtr<TopologicalRegionSet> > m_trsList;

  /// the names of the TopologicalRegionSet's to which the command applies to
  std::vector<std::string>                     m_trsNames;

  /// the name of the CommandGroup to which the command belongs
  std::string                                  m_groupName;

  /// pointer to the CommandGroup to which this command belongs
  Common::SafePtr<CommandGroup>              m_group;

  /// rate at which the processing has to be done
  CFuint m_processRate;
  
}; // class NumericalCommand

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(NumericalCommand) // declare this factory instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NumericalCommand_hh
