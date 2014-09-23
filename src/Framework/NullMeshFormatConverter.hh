// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullMeshFormatConverter_hh
#define COOLFluiD_Framework_NullMeshFormatConverter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"
#include "MeshFormatConverter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Provides an interface for a Null mesh format converter
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API NullMeshFormatConverter : public MeshFormatConverter {
public:

  /// Constructor
  NullMeshFormatConverter(const std::string& name) :
    MeshFormatConverter(name)
  {
  }

  /// Default destructor
  virtual ~NullMeshFormatConverter()
  {
  }

  /// Reads the all Dpl files and assembles the data in the converter.
  /// It is always in a convert or convertBack, but only exectuted once.
  /// @throw Common::FilesystemException if a file cannot be open
  /// @throw BadFormatException if a file is ill formated
  void readFiles(const boost::filesystem::path& filepath)
  {
    CFLog(VERBOSE,"NullMeshFormatConverter::readFiles() called!" << "\n");
  }

  /// Tries to check the file for conformity to the format.
  /// In a Null Mesfformat it doesn nothing other than
  /// issuing a wraning.
  void checkFormat(const boost::filesystem::path& filepath)
  {
    CFLog(VERBOSE,"NullMeshFormatConverter::checkFormat() called!" << "\n");
  }

  /// Writes the data read to the original format.
  /// Useful for debugging purposes.
  void convertBack(const boost::filesystem::path& filepath)
  {
    CFLog(VERBOSE,"NullMeshFormatConverter::convertBack() called!" << "\n");
  }

  /// Converts data from the file format to another format,
  /// taking into account the numerical method that will be
  /// used
  /// @param convertFromFileName name of the file to convert from
  /// @param fileName name of the file to write to
  void convert(const boost::filesystem::path& fromFilepath,
         const boost::filesystem::path& filepath)
  {
    CFLog(VERBOSE,"NullMeshFormatConverter::convert() called!" << "\n");
  }

protected:

  /// Adjust the node (state) numbering to make it stick to the
  /// COOLFluiD convention.
  void adjustToCFmeshNodeNumbering()
  {
    CFLog(VERBOSE,"NullMeshFormatConverter::adjustToCFmeshNodeNumbering() called!" << "\n");
  }

  /// Gets the target format.
  std::string getTargetFormat() const
  {
    return "Null";
  }

  /// Gets the origin format.
  std::string getOriginFormat() const
  {
    return "Null";
  }

  /// Checks if this object is a Null object.
  /// @return true since this is NullMeshFormatConverter
  bool isNull() const
  {
    return true;
  }

}; // end NullMeshFormatConverter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullMeshFormatConverter_hh
