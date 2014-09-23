// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FileFormatChecker_hh
#define COOLFluiD_Framework_FileFormatChecker_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <boost/filesystem/path.hpp>

#include "Common/NonCopyable.hh"
#include "Config/ConfigObject.hh"
#include "Common/NullableObject.hh"
#include "Common/OwnedObject.hh"
#include "Common/NullableObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides a basic interface for file readers
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API FileFormatChecker : public Common::OwnedObject,
                          public Config::ConfigObject,
                          public Common::NonCopyable<FileFormatChecker>,
                          public Common::NullableObject {
public:

  typedef Environment::ConcreteProvider<FileFormatChecker,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Default constructor without arguments
  FileFormatChecker(const std::string& name);

  /// Default destructor
  virtual ~FileFormatChecker();

  /// Checks the given file. This is a template method
  /// @throw Common::FilesystemException
  virtual void check(const boost::filesystem::path& filepath) = 0;

  /// Gets the Class name
  static std::string getClassName()
  {
    return "FileFormatChecker";
  }

}; // end of class FileFormatChecker

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(FileFormatChecker) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FileFormatChecker_hh
