// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_FileDownloader_hh
#define COOLFluiD_Common_FileDownloader_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NonCopyable.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// A class to handle file downloads
/// @author Tiago Quintino
class Common_API FileDownloader : public Common::NonCopyable<FileDownloader>
{
public: // methods

    /// Constructor
  FileDownloader();

    /// Destructor
  virtual ~FileDownloader();

    /// Download the file from the supplied URL and give it the supplied name
    /// @param url location of the file
    /// @param filepath file name to write
    /// @throw URLException if some problem occurs when accessing the URL
  virtual void download ( const std::string& url, const std::string& filepath ) = 0;

}; // class FileDownloader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_FileDownloader_hh
