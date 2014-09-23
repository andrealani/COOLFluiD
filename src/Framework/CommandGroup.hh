// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CommandGroup_hh
#define COOLFluiD_Framework_CommandGroup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"
#include "TopologicalRegionSet.hh"
#include "Environment/ConcreteProvider.hh"
#include "Common/SafePtr.hh"
#include "Common/OwnedObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a CommandGroup.
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class Framework_API CommandGroup : public Common::OwnedObject,
                     public Config::ConfigObject {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor
  explicit CommandGroup(const std::string& name);

  /// Default destructor
  ~CommandGroup();

  /// Gets the names of the TopologicalRegionSet's belonging to this CommandGroup
  /// These names are configured at startup time.
  const std::vector<std::string>& getTrsNames() const
  {
    return _trsNames;
  }

  /// Gets the names of the NumericalCommand's belonging to this CommandGroup
  /// These names are configured at startup time.
  const std::vector<std::string>& getComNames() const
  {
    return _comsNames;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "CommandGroup";
  }

protected:

  /// Get the TRS name
  const std::string getTrsName(const CFuint i)
  {
    return _trsNames[i];
  }

  /// Get the NumericalCommand name
  const std::string getComName(const CFuint i)
  {
    return _comsNames[i];
  }

private:

  /// the names of the TopologicalRegionSet's belonging to the CommandGroup
  std::vector<std::string>                     _trsNames;

  /// the names of the NumericalCommand's belonging to the CommandGroup
  std::vector<std::string>                     _comsNames;

}; // class CommandGroup

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CommandGroup_hh
