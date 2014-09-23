// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"

#include "Environment/DirPaths.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

DirPaths::DirPaths() :
  m_baseDir(),
  m_modulesDir(),
  m_workinDir(),
  m_resultsDir(),
  m_repURL()
{
}

//////////////////////////////////////////////////////////////////////////////

DirPaths::~DirPaths()
{
}

//////////////////////////////////////////////////////////////////////////////

DirPaths& DirPaths::getInstance()
{
  static DirPaths dirPaths;
  return dirPaths;
}

//////////////////////////////////////////////////////////////////////////////

void DirPaths::checkThrowMultipleBadDirException(const std::vector< boost::filesystem::path >& paths)
{
  if (!paths.empty()) {
    std::string totalPaths;
    vector< boost::filesystem::path >::const_iterator itr = paths.begin();
    for(; itr != paths.end(); ++itr) {
      totalPaths += "\'"+itr->string()+"\'";
      totalPaths += std::string(" ");
    }
    throw FilesystemException (FromHere(),"Following paths are not directories : " + totalPaths);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DirPaths::setDir(boost::filesystem::path& dpath, const std::string& dstr, const bool createdir)
{
  boost::filesystem::path d(dstr);

  if (!dstr.empty()) {
    if (createdir) {

      // The logic for directory creation is:
      // - path exists, not a directory: throw exception
      // - path exists, is a directory:  do nothing
      // - path does not exist:          create directory
      if (boost::filesystem::exists(d)) {
        if (!boost::filesystem::is_directory(d))  {
          throw FilesystemException(FromHere(),"Path \'"+dstr+"\'" + " exists and is not a directory");
        }
      } else {
        boost::filesystem::create_directory(d);
      }

    } else {

      if (!boost::filesystem::exists(d))  {
        throw FilesystemException(FromHere(),"Path \'"+dstr+"\'" + " does not exist");
      }

    }
  }
  dpath = d;
}

//////////////////////////////////////////////////////////////////////////////

void DirPaths::setBaseDir(const std::string& baseDir)
{
  CFLog(NOTICE,"Base Dir set to: \'" << baseDir << "\'\n");
  setDir(m_baseDir,baseDir);
}

//////////////////////////////////////////////////////////////////////////////

void DirPaths::addModuleDirs(const vector<std::string>& modulesDir)
{
  vector< boost::filesystem::path > badPaths;

  vector<std::string>::const_iterator itr = modulesDir.begin();
  for(; itr != modulesDir.end(); ++itr) {
    boost::filesystem::path p (*itr);
    if(boost::filesystem::exists(p))
    {
      if (boost::filesystem::is_directory(p)) {
        m_modulesDir.push_back(p);
        CFLog(VERBOSE,"Adding Module Dir: " << "\'" << p.string() << "\'\n");
      }
      else {
        badPaths.push_back(p);
      }
    }
  }

  checkThrowMultipleBadDirException(badPaths);
}

//////////////////////////////////////////////////////////////////////////////

void DirPaths::setWorkingDir(const std::string& workingDir)
{
  // if workingDir starts with "." or ".." then do not append m_baseDir
  boost::filesystem::path wDir;
  boost::filesystem::path workPath(workingDir);
  if (Common::StringOps::startsWith(workingDir,".")) {
    wDir = workPath;
  }
  else {
    wDir = m_baseDir / workPath;
  }

  CFLog(NOTICE,"Working Dir set to: \'" << wDir.string() << "\'\n");
  setDir(m_workinDir,wDir.string());
}

//////////////////////////////////////////////////////////////////////////////

void DirPaths::setRepositoryURL(const std::string& url)
{
  m_repURL = url;
}

//////////////////////////////////////////////////////////////////////////////

void DirPaths::setResultsDir(const std::string& resultsDir)
{
  // if resultsDir starts with "." or ".." then do not append m_baseDir
  boost::filesystem::path rDir;
  boost::filesystem::path resultsPath(resultsDir);
  if (Common::StringOps::startsWith(resultsDir,".")) {
    rDir = resultsPath;
  }
  else {
    rDir = m_baseDir / resultsPath;
  }

  CFLog(NOTICE,"Results Dir set to: \'" << rDir.string() << "\'\n");
  setDir(m_resultsDir,rDir.string(),true);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD
