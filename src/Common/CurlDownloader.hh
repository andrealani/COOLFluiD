// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CurlDownloader_hh
#define COOLFluiD_Common_CurlDownloader_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/progress.hpp>
#include <boost/filesystem/path.hpp>

#include "Common/FileDownloader.hh"
#include "Common/URLException.hh"

#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// A class to handle file downloads via the use of the Curl library
/// @author Tiago Quintino
class Common_API CurlDownloader : public FileDownloader
{
public: // methods

    /// Constructor
  CurlDownloader();

    /// Destructor
  virtual ~CurlDownloader();

    /// Download the file from the supplied URL and give it the supplied name.
    /// Uses the do_download() to do the actual job.
    /// Makes sure that in a parallel simulation only CPU with rank 0 does the downloading.
    /// @param url location of the file
    /// @param filepath file name to write
    /// @throw URLException if some problem occurs when accessing the URL
  virtual void download ( const std::string& url, const std::string& filepath );

private: // helper class

  /// File handle class to help interaction with Curl library.
  class FileHandle {
  public: // methods

      FileHandle() : fp(NULL), progress(NULL), curr(0) , tot(0), isopen(false) {}

      virtual ~FileHandle() { if (progress == NULL) delete progress; progress = NULL;}

  public: // data

      FILE * fp;
      boost::filesystem::path fpath;
      boost::progress_display * progress;
      unsigned int curr;
      unsigned int tot;
      bool isopen;

  }; // class FileHandle

private: // helper methods

    /// Actually do the work of downloading the file from the supplied URL and give it the supplied name
    /// @param url location of the file
    /// @param filepath file name to write
    /// @throw URLException if some problem occurs when accessing the URL
  void do_download ( const std::string& url, const std::string& filepath );

    /// Writes the data to the file handle as hte buffer is filled
  static int write_data(void * buffer, int size, int nmemb, void * data);

    /// Writes a progress bar as the file is doenloaded
  static int progress_func( class FileHandle * fh,
                     double dltotal,
                     double dlnow,
                     double ultotal,
                     double ulnow);

}; // class CurlDownloader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CurlDownloader_hh
