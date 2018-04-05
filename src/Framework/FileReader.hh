// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FileReader_hh
#define COOLFluiD_Framework_FileReader_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <boost/filesystem/path.hpp>

#include "Common/FilesystemException.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Common {
    class FactoryRegistry;
  }
  
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides a basic interface for file readers
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API FileReader {
public: // functions

  /// Default constructor without arguments
  FileReader();

  /// Default destructor
  virtual ~FileReader();

  /// Read the given file. This is a template method
  /// @throw Common::FilesystemException
  virtual void readFromFile(const boost::filesystem::path& filepath);

  /// Gets the file extension to append to the file name
  virtual const std::string getReaderFileExtension() const = 0;
 
  /// Set the factory registry
  virtual void setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr);
  
  /// Get the factory registry
  Common::SafePtr<Common::FactoryRegistry> getFactoryRegistry();
  
protected: // functions
  
  /// Read a string
  virtual bool readString(std::ifstream& file) {return false;}
  
  /// Get the name of the reader
  virtual const std::string getReaderName() const = 0;

  /// Terminator character
  virtual const std::string getReaderTerminator() const = 0;

  /// Hook method to call when reading is finished.
  /// Possibly to deallocate temporary data.
  virtual void finish();

  /// Set the flag telling to read again
  void setReadAgain(bool readAgain)
  {
    m_readAgain = readAgain;
  }

  /// Get the number of times the file has been read on end (one after the other)
  CFuint getReadCount() const
  {
    return m_readCount;
  }

 private: // data
 
  /// factory registry to allow polymorphic creation of objects
  Common::SafePtr<Common::FactoryRegistry> m_fr;
  
  /// flag telling if another reading is needed
  bool m_readAgain;

  /// reading counter
  CFuint m_readCount;

}; // end of class FileReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FileReader_hh
