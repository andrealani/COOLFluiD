// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_FileHandlerOutput_hh
#define COOLFluiD_Environment_FileHandlerOutput_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include <boost/filesystem/path.hpp>

#include "Common/NonCopyable.hh"
#include "Common/OwnedObject.hh"

#include "Environment/ConcreteProvider.hh"
#include "Environment/Factory.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// A class to handle files in a portable way.
/// @author Tiago Quintino
class Environment_API FileHandlerOutput :
    public Common::NonCopyable<Environment::FileHandlerOutput>,
    public Common::OwnedObject
{
public: // typedef

  /// typedef needed for the registration into the factory
  typedef Environment::ConcreteProvider<Environment::FileHandlerOutput> PROVIDER;

public: // methods

  /// Gets the Class name
  static std::string getClassName() { return "FileHandlerOutput"; }

  /// Constructor
  FileHandlerOutput();

  /// Destructor
  virtual ~FileHandlerOutput();

  /// Opens the file stream and returns the handle
  /// @pre isopen == false, no file should be open.
  /// @param filepath file name with path to be open
  /// @return a standard std::ofstream file handle opened
  /// @throw  FilesystemException if an error occurs while accessing the filesystem
  std::ofstream& open(const std::string& filepath, std::ios_base::openmode mode = std::ios_base::out);

  /// Opens the file stream and returns the handle
  /// @pre isopen == false, no file should be open.
  /// @param filepath file name with path to be open
  /// @return a standard std::ofstream file handle opened
  /// @throw  FilesystemException if an error occurs while accessing the filesystem
  virtual std::ofstream& open(const boost::filesystem::path& filepath, std::ios_base::openmode mode = std::ios_base::out) = 0;

  /// Closes the file stream.
  /// @post isopen == false, no file is open.
  virtual void close() = 0;

  /// Accesses the file stream.
  /// @pre isopen == true, file should be open.
  virtual std::ofstream& get() = 0;

  /// Checks if the file stream is open
  bool isopen() { return m_isopen; }

protected: // data

  /// flag if file is open
  bool m_isopen;

}; // class FileHandlerOutput

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Environment_Factory(FileHandlerOutput) // declare the factory of this type

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_FileHandlerOutput_hh
