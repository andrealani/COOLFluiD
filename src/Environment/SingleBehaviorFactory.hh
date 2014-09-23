// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_SingleBehaviorFactory_hh
#define COOLFluiD_Common_SingleBehaviorFactory_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SelfRegistPtr.hh"
#include "Environment/Factory.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

    template <class T> class Provider;

//////////////////////////////////////////////////////////////////////////////

/// This class provides the creation of objects of a certain abstract
/// TYPE polymorphically.
/// It behaves as a Factory but it should be used when the concrete
/// behavior should be exposed to the end-user. It should be used when
/// alternative algorithms are deterministic and produce the same exact result.
/// This is mostly convenient for the instantiation of template objects
/// configurered with different trait classes.
/// @author Tiago Quintino
template <class TYPE>
class SingleBehaviorFactory
{
public: // methods

  /// @return the instance of this singleton
  static SingleBehaviorFactory<TYPE>& getInstance();

  /// Function to create an instance polymorphically
  /// @return a pointer to the instantiated value
  Common::SelfRegistPtr<TYPE> create();

  /// Sets th default behavior to be constructed by this factory
  /// @param name of the provider which is default
  void setDefaultBehavior(const std::string& name);

  /// Sets th default behavior to be constructed by this factory
  /// @param name of the provider which is default
  std::string getDefaultBehavior();

    /// @return the name of the type of this factory
  virtual std::string getTypeName() const
  {
      return TYPE::getClassName();
  }

protected: // methods

  /// Constructor
  SingleBehaviorFactory();

  /// Default destructor
  virtual ~SingleBehaviorFactory();

private: // data

  std::string m_default;

}; // end of class SingleBehaviorFactory

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SingleBehaviorFactory.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Common_SingleBehaviorFactory_hh
