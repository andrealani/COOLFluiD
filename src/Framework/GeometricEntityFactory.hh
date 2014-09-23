// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GeometricEntityFactory_hh
#define COOLFluiD_Framework_GeometricEntityFactory_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"

#include "Framework/Framework.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

  class GeometricEntity;
  class State;
  class Node;
  class BaseGeometricEntityProvider;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a factory of geometric entities.
/// The Factory get information about the concrete GeometricEntityProvider
/// and register's it in the database.
/// When requested it will return the provider by its name.
/// The request of creating a GeometricEntity is should be sent
/// to the concrete provider.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API GeometricEntityFactory {
public:

  /// Register a BaseGeometricEntityProvider
  /// @param provider pointer to factory
  static void add(BaseGeometricEntityProvider* provider);

  /// Remove a registered factory
  static void remove(const std::string& providerName);

  /// Get a given Provider
  static BaseGeometricEntityProvider* getProvider(const std::string& providerName);

  /// Create a GeometricEntity
  static GeometricEntity* create(const std::string& providerName);

private:

  /// Constructor
  GeometricEntityFactory();

  /// Copy Constructor
  GeometricEntityFactory(const GeometricEntityFactory&);

  /// Destructor
  ~GeometricEntityFactory();

  /// Storage of pairs of type <name, provider>
  static std::map<std::string, BaseGeometricEntityProvider*>& getProviderMap();

}; // end of class GeometricEntityFactory

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GeometricEntityFactory_hh
