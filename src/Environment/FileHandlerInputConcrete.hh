// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_FileHandlerInputConcrete_hh
#define COOLFluiD_Environment_FileHandlerInputConcrete_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/fstream.hpp>

#ifdef CF_HAVE_BOOST_1_42
#include <boost/filesystem/exception.hpp>
#else
#include <boost/filesystem/operations.hpp>
#endif

#include "Environment/FileHandlerInput.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// A class to handle files in a portable way.
/// Makes use of boost::filesystem library.
/// @author Tiago Quintino
template < typename TRAIT_OPEN >
class FileHandlerInputConcrete : public FileHandlerInput
{
public: // methods

    /// Constructor
  FileHandlerInputConcrete();

    /// Destructor
  virtual ~FileHandlerInputConcrete();
  
  /// Opens the file stream and returns the handle
  /// @pre isopen == false, no file should be open.
  /// @param filepath file name with path to be open
  /// @return a standard std::ifstream file handle opened
  /// @throw  FilesystemException if an error occurs while accessing the filesystem
  virtual std::ifstream& open(const boost::filesystem::path& filepath,
			      std::ios_base::openmode mode);

  /// Opens the file stream and returns the handle
  /// @pre isopen == false, no file should be open.
  /// @param filepath file name with path to be open
  /// @return a standard std::ifstream file handle opened
  /// @throw  FilesystemException if an error occurs while accessing the filesystem
  virtual std::fstream& openBinary(const boost::filesystem::path& filepath,
				   std::ios_base::openmode mode);
  
  /// Closes the file stream.
  /// @post isopen == false, no file is open.
  virtual void close();
  
private: // data
  
  /// internal file handle from boost library
  boost::filesystem::ifstream fin;
  
  /// internal file handle from boost library
  boost::filesystem::fstream finBin;
  
}; // class FileHandlerInputConcrete

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FileHandlerInputConcrete.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_FileHandlerInputConcrete_hh
