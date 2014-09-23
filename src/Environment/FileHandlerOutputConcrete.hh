// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_FileHandlerOutputConcrete_hh
#define COOLFluiD_Environment_FileHandlerOutputConcrete_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>

#include "Environment/FileHandlerOutput.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// A class to handle files in a portable way.
/// Makes use of boost::filesystem library.
/// @author Tiago Quintino
template < typename TRAIT_OPEN >
class FileHandlerOutputConcrete : public FileHandlerOutput
{
public: // methods

    /// Constructor
  FileHandlerOutputConcrete();

    /// Destructor
  virtual ~FileHandlerOutputConcrete();

    /// Opens the file stream and returns the handle
    /// @pre isopen == false, no file should be open.
    /// @param filepath file name with path to be open
    /// @return a standard std::ifstream file handle opened
    /// @throw  FilesystemException if an error occurs while accessing the filesystem
  virtual std::ofstream& open(const boost::filesystem::path& filepath,
                               std::ios_base::openmode mode);

    /// Closes the file stream.
    /// @post isopen == false, no file is open.
  virtual void close();

    /// Accesses the file stream.
    /// @pre isopen == true, file should be open.
  virtual std::ofstream& get();

private: // data

  /// internal file handle from boost library
  boost::filesystem::ofstream fout;

}; // class FileHandlerOutputConcrete

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FileHandlerOutputConcrete.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_FileHandlerOutputConcrete_hh
