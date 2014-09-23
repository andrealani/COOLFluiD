// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullCommand_hh
#define COOLFluiD_Framework_NullCommand_hh

//////////////////////////////////////////////////////////////////////////////

#include "NumericalCommand.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullCommand.
/// Its simple intent on testing if the Command its applied to the correct
/// TRS
/// @author Tiago Quintino
class Framework_API NullCommand : public NumericalCommand {
public:

  /// Default constructor
  explicit NullCommand(const std::string& name);

  /// Default destructor
  ~NullCommand();

  /// Execute the command on a TopologicalRegionSet
  void execute();

  /// Check if this command needs a TRS to be supplied
  /// before being executed.
  bool usesTRS() const;

}; // class NullCommand

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullCommand_hh
