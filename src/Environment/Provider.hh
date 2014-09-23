// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_Provider_hh
#define COOLFluiD_Environment_Provider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NamedObject.hh"
#include "Common/ProviderBase.hh"
#include "Environment/Factory.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

  template < class BASE > class Factory;

//////////////////////////////////////////////////////////////////////////////

/// @brief Templated class provider
/// Used in the template abstract factory with Self Registering Objects.
/// @author Andrea Lani
/// @author Tiago Quintino
template <class BASE>
class Provider : public Common::NamedObject,
                 public Common::ProviderBase
{
public: // methods

  /// Constructor
  /// @param name provider registration name
  Provider(const std::string& name) :
    Common::NamedObject(name),
    Common::ProviderBase()
  {
    Environment::Factory<BASE>::getInstance().regist(this);
  }

  /// Virtual destructor
  virtual ~Provider() {}

  /// @return the name of this provider
  virtual std::string getProviderName () const { return getName(); }

  /// @return the BASE of this provider
  virtual std::string getProviderType () const { return BASE::getClassName(); }

}; // end of class Provider

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
#endif // COOLFluiD_Environment_Provider_hh
