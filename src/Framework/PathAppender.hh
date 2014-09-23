// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PathAppender_hh
#define COOLFluiD_Framework_PathAppender_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp>

#include "Common/NonCopyable.hh"
#include "Common/FilesystemException.hh"
#include "Common/CFLog.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a singleton object where
/// the paths get appended with specific simulation information.
/// This class is stateless and only possesses static functions.
/// This class is a Singleton pattern implementation.
/// @author Tiago Quintino
class Framework_API PathAppender : public Common::NonCopyable<PathAppender> {

public: // methods

  /// Appends all information to file name
  /// @param fpath path to be appended
  /// @return appended path
  boost::filesystem::path appendAllInfo(const boost::filesystem::path& fpath);

  /// Appends all information to file name
  /// @param fpath path to be appended
  /// @param appendIter force the iteration number to be appended to the file name (or not)
  /// @param appendTime force the time to be appended to the file name (or not)
  /// @return appended path
  boost::filesystem::path appendAllInfo(const boost::filesystem::path& fpath, 
                                        bool appendIterFlag, 
                                        bool appendTimeFlag,
                                        bool appendParallel = true);

  /// Appends the processor rank to the file path
  /// @param fpath path to be appended
  /// @return appended path
  boost::filesystem::path appendParallel(const boost::filesystem::path& fpath) const;

  /// Appends the current iteration to the file path
  /// @param fpath path to be appended
  /// @return appended path
  boost::filesystem::path appendIter(const boost::filesystem::path& fpath) const;

  /// Appends the current global iteration to the file path
  /// @param fpath path to be appended
  /// @return appended path
  boost::filesystem::path appendGlobalIter(const boost::filesystem::path& fpath) const;

  /// Appends the current time to the file path
  /// @param fpath path to be appended
  /// @return appended path
  boost::filesystem::path appendTime(const boost::filesystem::path& fpath) const;
  
  /// Appends a custom string to the file path
  /// @param fpath path to be appended
  /// @param custom A custom string to append
  /// @return appended path
  boost::filesystem::path appendCustom(const boost::filesystem::path& fpath, std::string& custom) const;

public: // methods for singleton

  /// @return the instance of this singleton
  static PathAppender& getInstance();

private: // methods

  /// Default constructor
  PathAppender();

  /// Default destructor
  ~PathAppender();

private: // no data define

  bool m_appendIterFlag;

  bool m_appendTimeFlag;
  
  bool m_appendParallelFlag; 

}; // end of class PathAppender

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PathAppender_hh
