// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshCellSplitter_CellSplitter3D_hh
#define COOLFluiD_IO_CFmeshCellSplitter_CellSplitter3D_hh

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
 * A class that is capable of subdividing
 *    a 3D CFmesh file (containing pyramids, prisms and hexahedras)
 * into 3D CFmesh file (containing only tetrahedras)
 *
 * @author Milan Zaloudek
 *
 * found at: Dompierre J., Labbe P., Vallet M.-G., Camarero R.: How to subdivide pyramids, prisms and hexahedra into tetrahedra, CERCA, Montreal
 *
 */
class CellSplitter3D : public Framework::MeshFormatConverter {
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
  CellSplitter3D(const std::string& name);

  /**
   * Destructor
   */
  virtual ~CellSplitter3D();

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
   * Converts current element to Tetrahedras
   */
  void convertElementsToTetra();

  /**
   * Updates the TRS data
   */
  void updateTRSData();

  /**
   * Change Element type information
   */
  void migrateElementTypes();

  /**
   * Converts current Pyramid to 2 Tetrahedras
   */
  void splitPyram(std::vector<CFuint> pyram);

  /**
   * Converts current Prism to 3 Tetrahedras
   */
  void splitPrism(std::vector<CFuint> prism);

  /**
   * Converts current Hexahedra to 4-5 Tetrahedras
   */
  void splitHexa(std::vector<CFuint> hexa);

  /**
   * Split Quads into triangles
   */
  void splitQuads(std::vector<CFuint> quad);

  /**
   * Searches for smallest node ID within current (arbitrary) face
   */
  void getSmallestID(CFuint nbNodes, std::vector<CFuint> OneElement);

  /**
   * Computes how many elements appear after splitting
   */
  void calculateNbNewElements();

  /**
   * Returns number of new tetras after splitting current hexa - either 5 or 6
   */
  CFuint nbNewTetras(CFuint iHexa);

private:

   /// the data to be read, extruded and rewriten to file
  /// this memory is owned here
  std::auto_ptr<Framework::CFmeshReaderWriterSource> _data;

  /// the file reader
  Reader _reader;

  /// the file writer
  Writer  _writer;

  /// smallest id of the (arbitrary) face
  CFint _SmallestFaceID;

  /// connectivity pattern
  std::valarray<CFuint> _pattern;

  /// id of the type TETRA
  CFuint _tetraTypeID;

  /// id of the type PYRAM
  CFuint _pyramTypeID;

  /// id of the type PRISM
  CFuint _prismTypeID;

  /// id of the type HEXA
  CFuint _hexaTypeID;

  /// how many new tetras are created
  CFuint _newTetrasSize;

  /// store the old element-node connectivity for
  /// constructing Top and Bottom TRS's
  Common::Table<CFuint>  _oldElemNode;

  /// store the old element-state connectivity for
  /// constructing Top and Bottom TRS's
  Common::Table<CFuint>  _oldElemState;

  /// temp Tetrahedras
  std::vector<std::vector<CFuint> >  _newTetras;
  /// temp Triangles
  std::vector<std::vector<CFuint> >  _newTriags;

}; // end class CellSplitter3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshCellSplitter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshCellSplitter_CellSplitter3D_hh
