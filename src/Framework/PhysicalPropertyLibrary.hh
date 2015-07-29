// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PhysicalPropertyLibrary_hh
#define COOLFluiD_Framework_PhysicalPropertyLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "Common/NonCopyable.hh"
#include "Common/OwnedObject.hh"

#include "Config/ConfigObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Provides an abstract interface for libraries that compute the physical properties.
/// @author Tiago Quintino
class Framework_API PhysicalPropertyLibrary :
      public Common::NonCopyable<PhysicalPropertyLibrary>,
      public Common::OwnedObject,
      public Common::SetupObject,
      public Common::NullableObject,
      public Config::ConfigObject  {

public: // typedefs

  typedef Environment::ConcreteProvider<PhysicalPropertyLibrary,1> PROVIDER;
  typedef const std::string& ARG1;
 public: // functions
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Constructor
  PhysicalPropertyLibrary(const std::string& name);

  /// Default destructor
  virtual ~PhysicalPropertyLibrary();
  
  /// Configures this configurable object.
  virtual void configure ( Config::ConfigArgs& args );
  
  /// Gets the Class name
  static std::string getClassName() { return "PhysicalPropertyLibrary"; }
  
  /// Setups the path name
  void setLibPathName(const std::string& libPathName)
  {
    m_libPath = libPathName;
  }
  
protected:
  
  /// path to some input library data
  std::string m_libPath;
  
}; // end of class PhysicalPropertyLibrary

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(PhysicalPropertyLibrary)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PhysicalPropertyLibrary_hh
