// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_DirectFileAccess_hh
#define COOLFluiD_Environment_DirectFileAccess_hh

//////////////////////////////////////////////////////////////////////////////

#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/filesystem/exception.hpp"

#include "Common/CFLog.hh"
#include "Common/NonInstantiable.hh"
#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// A trait class for FileHandlerConcrete to change the behavior of the
/// FileHandlerConcrete::open() method in order to access directly the
/// file on the filesystem.
/// @author Tiago Quintino
class Environment_API DirectFileAccess : public Common::NonInstantiable<DirectFileAccess>
{
public: // methods

  template <typename TYPE>
  static TYPE& open(TYPE& fin,
		    const boost::filesystem::path& filepath,
		    std::ios_base::openmode mode);
  
}; // class DirectFileAccess

//////////////////////////////////////////////////////////////////////////////

template <typename TYPE>
TYPE& DirectFileAccess::open(TYPE& fin,
			     const boost::filesystem::path& filepath,
			     std::ios_base::openmode mode)
{
  boost::filesystem::path fp (filepath);
  
  // if the file is present open it
  // does it fail if it is  a directory ?
  if( boost::filesystem::exists(fp) )
    {
      CFLog(VERBOSE, "DirectFileAccess::Opening file " <<  fp.string() << "\n");
      fin.open(fp,mode); // exists so open it
    }
  else // doesnt exist so throw exception
    {
      throw boost::filesystem::filesystem_error( fp.string() + " does not exist for DirectFileAccess",
						 boost::system::error_code() );
    }
  return fin;
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_DirectFileAccess_hh
