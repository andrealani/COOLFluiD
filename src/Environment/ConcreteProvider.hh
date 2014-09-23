// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_ConcreteProvider_hh
#define COOLFluiD_Environment_ConcreteProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/SelfRegistPtr.hh"

#include "Environment/Provider.hh"
#include "Environment/Factory.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// @brief Concrete provider class for all types which the constructor takes zero arguments
/// @author Andrea Lani
/// @author Tiago Quintino
template <class BASE, CFint NBARG = 0>
class ConcreteProvider : public Environment::Provider<BASE> {
public:

  /// Constructor
  explicit ConcreteProvider(const std::string& name) :  Environment::Provider<BASE>(name) {}

  /// Virtual destructor
  virtual ~ConcreteProvider() {}

  /// Polymorphic function to create objects of dynamical type BASE
  /// @return SelfRegistPtr olding the created object
  virtual Common::SelfRegistPtr<BASE> create() = 0;

}; // end of class ConcreteProvider

//////////////////////////////////////////////////////////////////////////////

/// @brief Concrete provider class for all types which the constructor takes one argument
/// @author Andrea Lani
/// @author Tiago Quintino
template <class BASE>
class ConcreteProvider<BASE,1> : public Environment::Provider<BASE> {
public:

  typedef BASE BASE_TYPE;
  typedef typename BASE::ARG1 BASE_ARG1;

  /// Constructor
  explicit ConcreteProvider(const std::string& name) :  Environment::Provider<BASE>(name)  {}

  /// Virtual destructor
  virtual ~ConcreteProvider() {}

  /// Polymorphic function to create objects of dynamical type BASE
  /// @param arg1 first parameter
  /// @return SelfRegistPtr olding the created object
  virtual Common::SelfRegistPtr<BASE> create(BASE_ARG1 arg1) = 0;

}; // end of class ConcreteProvider

//////////////////////////////////////////////////////////////////////////////

/// @brief Concrete provider class for all types which the constructor takes two argument
/// @author Andrea Lani
/// @author Tiago Quintino
template <class BASE>
class ConcreteProvider<BASE,2> : public Environment::Provider<BASE> {
public:

  typedef BASE BASE_TYPE;
  typedef typename BASE::ARG1 BASE_ARG1;
  typedef typename BASE::ARG2 BASE_ARG2;

  /// Constructor
  explicit ConcreteProvider(const std::string& name) : Environment::Provider<BASE>(name) {}

  /// Virtual destructor
  virtual ~ConcreteProvider() {}

  /// Polymorphic function to create objects of dynamical type BASE
  /// @param arg1 first parameter
  /// @param arg2 first parameter
  /// @return SelfRegistPtr olding the created object
  virtual Common::SelfRegistPtr<BASE> create(BASE_ARG1 arg1, BASE_ARG2 arg2) = 0;

}; // end of class ConcreteProvider

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_ConcreteProvider_hh
