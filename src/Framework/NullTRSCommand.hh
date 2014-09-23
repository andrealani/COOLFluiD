// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullTRSCommand_hh
#define COOLFluiD_Framework_NullTRSCommand_hh

//////////////////////////////////////////////////////////////////////////////

#include "NumericalCommand.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullTRSCommand.
/// Its simple intent on testing if the Command its applied to the correct
/// TRS
/// @author Tiago Quintino
class Framework_API NullTRSCommand : public NumericalCommand {
public:

  /// Default constructor
  explicit NullTRSCommand(const std::string& name);

  /// Default destructor
  virtual ~NullTRSCommand();

  /// Execute the command on a TopologicalRegionSet
  virtual void execute();

  /// Check if this command needs a TRS to be supplied
  /// before being executed.
  bool usesTRS() const;

}; // class NullTRSCommand

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullTRSCommand_hh
