// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MethodCommand_hh
#define COOLFluiD_Framework_MethodCommand_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SharedPtr.hh"
#include "Common/NullableObject.hh"

#include "Common/CFLog.hh"

#include "Framework/BaseMethodCommandProvider.hh"
#include "Framework/NumericalCommand.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand to be sent by a Method.
/// It is a template intended to function as the command that it is
/// parametrized with the DATA object shared between all MethodCommands
/// belonging to the same Method.
/// @author Tiago Quintino
/// @see NumericalCommand
template <class DATA>
class MethodCommand : public NumericalCommand,
                      public Common::NullableObject
{
public:

  typedef MethodCommand<DATA> MET_COMMAND;
  typedef BaseMethodCommandProvider<DATA, MET_COMMAND> PROVIDER;

  /// Constructor.
  /// @see NumericalCommand()
  explicit MethodCommand(const std::string& name) : NumericalCommand(name) {}

  /// Default destructor
  virtual ~MethodCommand() {}

  /// Gets the shared DATA object for the Method.
  /// @return reference to the Method's Data
  DATA& getMethodData() { return *m_data_ptr; }

  /// Sets the pointer to the shared DATA object
  /// @param dataPtr is the pointer to the Data to be set.
  void setMethodData ( const Common::SharedPtr<DATA>& ptr ) {  m_data_ptr = ptr; }

  /// Checks validity of _dataPtr.
  bool hasMethodData() const { return m_data_ptr.isNull(); }

  /// Gets the Class name
  static std::string getClassName() { return DATA::getClassName() + "Command"; }

private:

  /// Pointer to the data object
  Common::SharedPtr<DATA> m_data_ptr;

}; // end of class MethodCommand

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Null NumericalCommand to be sent by a Method.
/// @author Tiago Quintino
/// @see MethodCommand
template <class DATA>
class NullMethodCommand : public MethodCommand<DATA>    {
public:

  /// Constructor.
  /// @see NumericalCommand()
  explicit NullMethodCommand(const std::string& name) : MethodCommand<DATA>(name) {}

  /// Default destructor
  virtual ~NullMethodCommand() {}

  /// Checks if this object is a Null object.
  /// Sincethis is a Null command, it returns true
  /// @return true if Null and false otherwise
  virtual bool isNull() const { return true; }

  /// Check if this command needs a TRS to be supplied
  /// before being executed.
  virtual bool usesTRS() const { return false; }

  /// Execute Processing actions
  virtual void execute()
  {
    CFLogDebugMin ( MethodCommand<DATA>::getClassName() << "::execute() called on a Null command!" << "\n");
  }

}; // end of class NullMethodCommand

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MethodCommand_hh
