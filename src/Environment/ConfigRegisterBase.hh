// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOFluiD_Environment_ConfigRegisterBase_hh
#define COOFluiD_Environment_ConfigRegisterBase_hh

#include "Common/NonCopyable.hh"
#include "Common/COOLFluiD.hh"
#include "Environment/EnvironmentAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// This class servers as a base class for ConfigRegister
/// @author Tiago Quintino
class Environment_API ConfigRegisterBase : public Common::NonCopyable <ConfigRegisterBase> {
public:

  /// @return the name of the type of this factory
  virtual std::string getTypeName() const = 0;

protected:

  /// Constructor
  ConfigRegisterBase();

  /// virtual destructor
  virtual ~ConfigRegisterBase();

}; // end class ConfigRegisterBase

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Environment_ConfigRegisterBase_hh
