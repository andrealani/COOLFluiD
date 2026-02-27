// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_DirPaths_hh
#define COOLFluiD_Environment_DirPaths_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
#ifdef CF_HAVE_BOOST_1_85
#include <boost/filesystem.hpp>
#else
#include <boost/filesystem/convenience.hpp>
#endif 

#include "Common/NonCopyable.hh"
#include "Common/FilesystemException.hh"

#include "Config/ConfigObject.hh"

#include "Environment/EnvironmentAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a singleton object where
/// the paths for certain directories are stored.
/// This class is a Singleton pattern implementation.
/// @author Tiago Quintino
class Environment_API DirPaths : public Common::NonCopyable<DirPaths> {

public: // methods

  /// Set the Base dir path
  /// @throw FilesystemException if dir does not exist
  void setBaseDir(const std::string& baseDir);

  /// Set the repository
  void setRepositoryURL(const std::string& url);

  /// Set the Modules dir paths
  /// @throw FilesystemException if one dir does not exist
  void addModuleDirs(const std::vector<std::string>& modulesDir);

  /// Set the Working dir path
  /// @throw FilesystemException if dir does not exist
  void setWorkingDir(const std::string& workingDir);

  /// Set the Results dir path
  /// @throw FilesystemException if dir does not exist
  void setResultsDir(const std::string& resultsDir);

  /// Get the repository URL
  std::string getRepositoryURL() const { return m_repURL; }

  /// Get the Working dir path
  boost::filesystem::path getWorkingDir() const { return m_workinDir; }

  /// Get the Results dir path
  boost::filesystem::path getResultsDir() const {  return m_resultsDir; }

  /// Get the Base dir path
  boost::filesystem::path getBaseDir() const {  return m_baseDir; }

  /// Get the Modules dir paths.
  std::vector< boost::filesystem::path > getModulesDir() const { return m_modulesDir; }

public: // methods for singleton

  /// @return the instance of this singleton
  static DirPaths& getInstance();

protected: // methods

  /// Sets a dir path
  /// @param dpath the path to be set (usually a private variable of the class)
  /// @param dstr the string to convert to the path
  /// @param createdir if directory does not exist, create it
  /// @throw FilesystemException if dir does not exist, and createdir is false
  void setDir(boost::filesystem::path& dpath, const std::string& dstr, const bool createdir=false);

  /// Throw exception for all the supplied paths
  void checkThrowMultipleBadDirException(const std::vector< boost::filesystem::path >& paths);

private: // methods

  /// Default constructor
  DirPaths();

  /// Default destructor
  ~DirPaths();

private: // data

  /// base dir path
  boost::filesystem::path m_baseDir;
  /// the modules directories paths where to search for module libs
  std::vector< boost::filesystem::path > m_modulesDir;
  /// the working dir path
  boost::filesystem::path m_workinDir;
  /// the results dir path
  boost::filesystem::path m_resultsDir;
  /// the url for file repository
  std::string m_repURL;

}; // end of class DirPaths

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_DirPaths_hh
