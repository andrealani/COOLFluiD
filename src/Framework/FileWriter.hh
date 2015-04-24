// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FileWriter_hh
#define COOLFluiD_Framework_FileWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include <boost/filesystem/path.hpp>

#include "Common/StringOps.hh"
#include "Common/FilesystemException.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides a basic interface for file writers
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API FileWriter {
public:

  /// Default constructor without arguments
  FileWriter();

  /// Default destructor
  virtual ~FileWriter();
  
  /// Opens and starts to write to the given file.
  /// @throw Common::FilesystemException
  virtual void writeToFile(const boost::filesystem::path& filepath);
  
  /// Gets the file extension to append to the file name
  virtual const std::string getWriterFileExtension() const = 0;

protected:

  /// Write the given file. This is an pure abstract method
  /// @throw Common::FilesystemException
  virtual void writeToFileStream(std::ofstream& fout) 
  {
    throw Common::NotImplementedException (FromHere(),"FileWriter::writeToFileStream()");
  }

  /// Get the name of the writer
  virtual const std::string getWriterName() const = 0;

}; // end of class FileWriter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FileWriter_hh
