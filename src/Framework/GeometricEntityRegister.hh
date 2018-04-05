// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GeometricEntityRegister_hh
#define COOLFluiD_Framework_GeometricEntityRegister_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/NoSuchValueException.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {
    class FactoryRegistry;
  }

  namespace Framework {

  class BaseGeometricEntityProvider;

//////////////////////////////////////////////////////////////////////////////

/// This class provides a class that has knowledge about
/// which ShapeSunction's are in use in the current SubSystem.
/// This differs from a Factory of Provider's as the information
/// stored is changes with each different simulation options.
/// This class is a Singleton.
/// @author Tiago Quintino
class Framework_API GeometricEntityRegister {
private: // typedefs

  /// entry type in database
  typedef BaseGeometricEntityProvider* EntryType;

  /// Database type
  typedef std::vector<EntryType> DatabaseType;

public: // methods

  /// Gets the instance of the Singleton
  static GeometricEntityRegister& getInstance();

  /// Set the factory registry
  void setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr);
 
  /// Get the factory registry
  Common::SafePtr<Common::FactoryRegistry> getFactoryRegistry();
  
  /// Regists this type of shape function as being used
  /// The returned index can be later used in the getProvider() function
  /// @return the index which corresponds to the registered type
  CFuint regist(const std::string& providerName);

  /// Gets the Provider of type of shape function by the its name
  /// @return the provider pointer
  /// @post returns CFNULL if it is not on the database
  BaseGeometricEntityProvider* getProvider(const std::string& providerName);

  /// Gets the Provider of type of shape function by the its index
  /// @return the provider pointer
  BaseGeometricEntityProvider* getProvider(const CFuint idx);

  /// Checks the name of Provider of shape function if exists in database
  bool exists(const std::string& providerName);

  /// Checks the Provider of shape function if exists in database
  bool exists(BaseGeometricEntityProvider* provider);

  /// Get he size of the register
  CFuint getSize() const
  {
    return _database.size();
  }

  /// Clears the registry
  void clear();

private: // methods

  /// Constructor.
  /// Private to make it a Singleton.
  GeometricEntityRegister();

  /// Destructor
  /// Private to make it a Singleton.
  ~GeometricEntityRegister();

  /// Copy Constructor.
  /// Private to make it a Singleton.
  GeometricEntityRegister(const GeometricEntityRegister&);

  /// Copy Operator.
  /// Private to make it a Singleton.
  GeometricEntityRegister& operator=(const GeometricEntityRegister&);

private: // data

  // the shape function register database
  DatabaseType _database;

  /// factory registry to allow polymorphic creation of objects
  Common::SafePtr<Common::FactoryRegistry> m_fr;

}; // end class GeometricEntityRegister

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GeometricEntityRegister_hh
