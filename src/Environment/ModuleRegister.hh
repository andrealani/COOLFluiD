// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOFluiD_Environment_ModuleRegister_hh
#define COOFluiD_Environment_ModuleRegister_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegisterBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a module register template.
/// @author Tiago Quintino
template < typename MODULE >
class ModuleRegister : public Environment::ModuleRegisterBase {
public: // methods

    /// Acessor to the singleton
  static MODULE& getInstance();

  /// Returns the description of the module.
  /// Must be implemented by the ModuleRegister
  /// @return descripton of the module
  virtual std::string getDescription() const;

protected: // methods

    /// Constructor
  ModuleRegister();

    /// Destructor
  ~ModuleRegister();

}; // end class ModuleRegister

//////////////////////////////////////////////////////////////////////////////

template < typename MODULE >
MODULE& ModuleRegister<MODULE>::getInstance()
{
  static MODULE instance;
  return instance;
}

//////////////////////////////////////////////////////////////////////////////

template < typename MODULE >
ModuleRegister<MODULE>::ModuleRegister() :
 ModuleRegisterBase(MODULE::getModuleName())
{
}

//////////////////////////////////////////////////////////////////////////////

template < typename MODULE >
ModuleRegister<MODULE>::~ModuleRegister()
{
}

//////////////////////////////////////////////////////////////////////////////

template < typename MODULE >
std::string ModuleRegister<MODULE>::getDescription() const
{
  return MODULE::getModuleDescription();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Environment_ModuleRegister_hh
