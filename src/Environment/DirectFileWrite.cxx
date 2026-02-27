// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifdef CF_HAVE_BOOST_1_85
#include "boost/filesystem.hpp"
#else
#include "boost/filesystem/convenience.hpp"
#endif
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/filesystem/exception.hpp"

#include "Environment/FileHandlerOutputConcrete.hh"

#include "Environment/Environment.hh"
#include "Environment/DirectFileWrite.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Environment {

//////////////////////////////////////////////////////////////////////////////

// provider for this behavior
Environment::ObjectProvider<
                Environment::FileHandlerOutputConcrete< DirectFileWrite >,
                Environment::FileHandlerOutput,
                EnvironmentModule >
Provider_DirectFileWrite("DirectFileWrite");

//////////////////////////////////////////////////////////////////////////////

std::ofstream& DirectFileWrite::open(boost::filesystem::ofstream& fout,
                                     const boost::filesystem::path& filepath,
                                     std::ios_base::openmode mode)
{
   boost::filesystem::path fp (filepath);

   CFLog(VERBOSE, "Opening file " <<  fp.string() << "\n");
   fout.open(fp,mode);
   if (!fout) // didn't open so throw exception
   {
      throw boost::filesystem::filesystem_error( fp.string() + " failed to open",
                                                 boost::system::error_code() );
   }
   return fout;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Environment

} // namespace COOLFluiD
