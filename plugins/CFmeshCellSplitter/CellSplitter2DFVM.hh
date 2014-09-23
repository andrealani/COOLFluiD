// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshCellSplitter_CellSplitter2DFVM_hh
#define COOLFluiD_IO_CFmeshCellSplitter_CellSplitter2DFVM_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>


#include "Common/Table.hh"
#include "Framework/MeshFormatConverter.hh"
#include "Framework/CFmeshFileReader.hh"
#include "Framework/CFmeshFileWriter.hh"
#include "Framework/CFmeshReaderWriterSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace CFmeshCellSplitter {

//////////////////////////////////////////////////////////////////////////////

/**
 * A class that is capable of converting a 2D CFmesh file
 * to a 2D CFmesh file with split quads
 *
 * @author Thomas Wuilbaut
 *
 */
class CellSplitter2DFVM : public Framework::MeshFormatConverter {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::CFmeshFileReader
  <Framework::CFmeshReaderWriterSource> Reader;

  typedef Framework::CFmeshFileWriter
  <Framework::CFmeshReaderWriterSource> Writer;

  /**
   * Constructor
   */
  CellSplitter2DFVM(const std::string& name);

  /**
   * Destructor
   */
  virtual ~CellSplitter2DFVM();

  /**
   * Reads the all Dpl files and assembles the data in the converter.
   * It is always in a convert or convertBack, but only exectuted once.
   *
   * @throw Common::FilesystemException if a file cannot be open
   * @throw BadFormatException if a file is ill formated
   */
  void readFiles(const boost::filesystem::path& filepath) {}

  /**
   * Tries to check the file for conformity to the format.
   * Possibly not full proof.
   */
  void checkFormat(const boost::filesystem::path& filepath);

  /**
   * Writes the data read to the original format.
   * Useful for debugging purposes.
   */
  void convertBack(const boost::filesystem::path& filepath);

  /**
   * Converts data from the file format to another format,
   * taking into account the numerical method that will be
   * used
   * @param convertFromFileName name of the file to convert
   * @param fileName name of the file to convert
   */
  void convert(const boost::filesystem::path& fromFilepath,
         const boost::filesystem::path& filepath);

  /**
   * Configures this object.
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

protected:

  /**
   * Adjust the node (state) numbering to make it stick to the
   * COOLFluiD convention.
   */
  void adjustToCFmeshNodeNumbering() {}

  /**
   * Gets the target format.
   */
  std::string getTargetFormat() const
  {
    return "CFmesh";
  }

  /**
   * Gets the origin format.
   */
  std::string getOriginFormat() const
  {
    return "CFmesh";
  }

private:

  /**
   * Splits the elements
   */
  void split();

  /**
   * Converts current elements in 2D to 3D
   */
  void transformQuadElements();

  /**
   * Converts current quads to triangles
   */
  void convertElementsToTriag();

  /**
   * Change Element type information
   */
  void migrateElementTypes();

  /**
   * Change the face-state connectivity for the TRSs
   */
  void convertTRSConnectivity();

  /**
   * Split Quads into triangles
   */
  void splitQuads(std::vector<CFuint> quad);


private:

   /// the data to be read, extruded and rewriten to file
  /// this memory is owned here
  std::auto_ptr<Framework::CFmeshReaderWriterSource> _data;

  /// the file reader
  Reader _reader;

  /// the file writer
  Writer  _writer;

  /// connectivity pattern
  std::valarray<CFuint> _nodePattern;
  std::valarray<CFuint> _statePattern;

  /// id of the type TRIAG
  CFuint _triagTypeID;

  /// id of the type QUAD
  CFuint _quadTypeID;

  /// store the old element-node connectivity for
  /// constructing Top and Bottom TRS's
  Common::Table<CFuint>  _oldElemNode;

  /// store the old element-state connectivity for
  /// constructing Top and Bottom TRS's
  Common::Table<CFuint>  _oldElemState;

  /// temp Triangles
  std::vector<std::vector<CFuint> >  _newTriags;

}; // end class CellSplitter2DFVM

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshCellSplitter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshCellSplitter_CellSplitter2DFVM_hh
