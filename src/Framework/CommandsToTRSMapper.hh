// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CommandsToTRSMapper_hh
#define COOLFluiD_Framework_CommandsToTRSMapper_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

  class NumericalCommand;
  class TopologicalRegionSet;

//////////////////////////////////////////////////////////////////////////////

/// This class is capable of mapping commands to their TopologicalRegionSet's
/// @author Tiago Quintino
class Framework_API CommandsToTRSMapper {
public: // functions

  /// Default constructor without arguments.
  CommandsToTRSMapper(std::vector< Common::SafePtr<NumericalCommand> >& coms,
                      std::vector< Common::SafePtr<TopologicalRegionSet> >& trs);

  /// Default destructor.
  ~CommandsToTRSMapper();

  /// Map the Commands to the TRS using the configurated list of TRS names
  void mapComsToTrs() const;

private: // helper functions

  /// Helper function that actually does the job for each NumericalCommand
  void processCommand( Common::SafePtr<NumericalCommand> comPtr) const;

private: // data

  /// Reference to the list of Commands to map
  std::vector< Common::SafePtr<NumericalCommand> >& m_coms;

  /// Reference to the list of Trs's to map
  std::vector< Common::SafePtr<TopologicalRegionSet> >& m_trs;

}; // end of class CommandsToTRSMapper

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CommandsToTRSMapper_hh
