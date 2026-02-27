// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_CurlAccessRepository_hh
#define COOLFluiD_Environment_CurlAccessRepository_hh

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_BOOST_1_85
#include "boost/filesystem.hpp"
#else
#include "boost/filesystem/convenience.hpp"
#endif
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/filesystem/exception.hpp"

#include "Common/CFLog.hh"
#include "Common/NonInstantiable.hh"
#include "Common/CurlDownloader.hh"
#include "Common/URLException.hh"
#include "Common/COOLFluiD.hh"
#include "Environment/DirPaths.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// A trait class for FileHandlerConcrete to change the behavior of the
/// FileHandlerConcrete::open() method in order to try first to access
/// the file on the disk, and if not found to try to retrieve it from a
/// remote repository by using the curl library.
/// @author Tiago Quintino
class Environment_API CurlAccessRepository : public Common::NonInstantiable<CurlAccessRepository>
{
public: // methods
  
  template <typename TYPE>
  static TYPE& open(TYPE& fin,
		    const boost::filesystem::path& filepath,
		    std::ios_base::openmode mode);
  
}; // class CurlAccessRepository

//////////////////////////////////////////////////////////////////////////////

template <typename TYPE>
TYPE& CurlAccessRepository::open(TYPE& fin,
				 const boost::filesystem::path& fp,
				 std::ios_base::openmode mode)
{
  using namespace Common;
  
  // download if the file is not present
  if( !boost::filesystem::exists(fp) )
  {
#ifdef CF_HAVE_BOOST_1_85
    boost::filesystem::path bp = fp.parent_path();
#else 
    boost::filesystem::path bp = fp.branch_path();
#endif   
    if (!bp.string().empty())
      {
	// make sure directory exists
	if ( ! boost::filesystem::exists( bp ) )
	  {
	    boost::filesystem::create_directories(bp);
	  }
	
	// make sure it is a directory
	if ( ! boost::filesystem::is_directory( bp ) )
	  {
	    throw boost::filesystem::filesystem_error( (bp.string() + "is not a directory") , boost::system::error_code() );
	  }
      }
    
    // download the file
    FileDownloader * dl = new CurlDownloader();
    
    // remove the working dir from the file name
    std::string furl = fp.string();
    Common::StringOps::subst(Environment::DirPaths::getInstance().getBaseDir().string(), "", furl);
    
    std::string urlpath = Environment::DirPaths::getInstance().getRepositoryURL() + furl;
    CFLog(NOTICE, "File " << fp.string() << " was not found.\nRetrieving from " << urlpath << "\n");
    
    // add error control for download failure
    try {
      dl->download(urlpath,fp.string());
    }
    catch (URLException& e) // convert to an exception the client code is expecting
      {
        throw  boost::filesystem::filesystem_error( furl + " is not accessible in the repository\n",
                                                    boost::system::error_code() );
      }
    
    delete dl;
  }
  else
    {
      CFLog(VERBOSE, "File " <<  fp.string() << " exists. Not retrieving." << "\n");
    }
  
   // open it
  CFLog(VERBOSE, "Opening file " <<  fp.string() << "\n");
  fin.open(fp,mode);
  return fin;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_CurlAccessRepository_hh
