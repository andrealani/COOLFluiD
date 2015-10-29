// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshExtruder_Extruder2DFVMMPI_hh
#define COOLFluiD_IO_CFmeshExtruder_Extruder2DFVMMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include <mpi.h>
#include <fstream>

#include "Common/Table.hh"
#include "Framework/MeshFormatConverter.hh"
#include "Framework/CFmeshFileReader.hh"
#include "Framework/CFmeshBinaryFileWriter.hh"
#include "Framework/CFmeshReaderWriterSource.hh"
#include "MathTools/FunctionParser.hh"
#include "Common/ParserException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace CFmeshExtruder {

//////////////////////////////////////////////////////////////////////////////

/**
 * A class that is capable of converting a 2D CFmesh file to 3D CFmesh file 
 * by extruding the mesh in z coordinate in parallel.
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class Extruder2DFVMMPI : public Framework::MeshFormatConverter {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  typedef Framework::CFmeshFileReader
  <Framework::CFmeshReaderWriterSource> Reader;
  
  typedef Framework::CFmeshBinaryFileWriter
  <Framework::CFmeshReaderWriterSource> Writer;
  
  /**
   * Constructor
   */
  Extruder2DFVMMPI(const std::string& name);

  /**
   * Destructor
   */
  virtual ~Extruder2DFVMMPI();

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
  
  /// Tell if this conveter can work in parallel.
  virtual bool isParallel() const {return true;}
  
protected:

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
   * Adjust the node (state) numbering to make it stick to the
   * COOLFluiD convention.
   */
  void adjustToCFmeshNodeNumbering() {}

  /**
   * Extrudes the data form 2D to 3D
   */
  void extrude();

  /**
   * Transforms the current nodes from 2D to 3D coordinates
   */
  void transformNodesTo3D();

  /**
   * Converts current elements in 2D to 3D
   */
  void convertCurrentElementsTo3D();

  /**
   * Change Element type information to 3D
   */
  void migrateElementTypes();

  /**
   * Creates the first layer of 3D elements
   * from the initial 2D information.
   * This is used prior to the following iterative procedure.
   * Due to algorithmic differences, the first layer cannot be done
   * in the same manner as the other ones.
   */
  void createFirstLayer();

  /**
   * Creates the next layers of 3D elements
   * from the first layer.
   * This is used in a iterative procedure.
   */
  void createAnotherLayer();

  /**
   * Updates the TRS data
   */
  void updateTRSData();

  /**
   * Extrudes the data in the current TRS's
   */
  void extrudeCurrentTRSs();

  /**
   * Creates a TRS in th top and bottom of the extrusion
   */
  void createTopBottomTRSs();

  /**
   * Creates a side TRS
   * @param name the name to give the TRS
   * @param layer the ID of the layers(1 or more) of nodes and states on which to base the
   *              GeometricEntity's of the TRS
   * @param bottom is a boolean that indicates to the function if the TRS to create
   *               is the bottom or the top one.
   */
  void createSideTRS(const std::string& name,
		     const std::vector<CFuint>& layer,
		     const bool& bottom);
  
  /**
   * Split 3D elements into tetrahedras
   */
  void splitPrism(const std::vector<CFuint>& prism);
  
  /**
   * Random z coordinate for nodes in inner layers
   */
  void randomNodes();

  /**
   * Split Quads into triangles
   */
  void splitQuads(const std::vector<CFuint>& quad);

  /**
   * Dirty hack if hybrid
   */
//   void reorderElementNodeState();

  /**
   * This struct is used for temporary definition of a prism
   */
  struct PrismData {
    /// Coordinates
    RealVector coord;

    /// ID
    RealVector prismID;
  }; // end PhysicalData
  
  /// get the deltaZ corresponding to the given x,y,l
  /// @param x  x-coordinate of the mesh point
  /// @param y  y-coordinate of the mesh point
  /// @param l  global ID of the current layer  
  CFreal getDeltaZ(const CFreal x, const CFreal y, const CFreal l)
  {
    if (_function != "0") {
      _eval[0] = x; _eval[1] = y; _eval[2] = l; 
      return _functionParser.Eval(&_eval[0]);
    }
    return _zDelta;
  }
  
private:
    
   /// the data to be read, extruded and rewriten to file
  /// this memory is owned here
  std::auto_ptr<Framework::CFmeshReaderWriterSource> _data;
 
  /// temp Tetras
  std::vector<std::vector<CFuint> >  _newTetras;
  
  /// temp Triangles
  std::vector<std::vector<CFuint> >  _newTriag;
  
  /// CFL function of iteration
  MathTools::FunctionParser _functionParser;
  
  /// a temporary vector to store CFL
  /// and pass the function parser
  RealVector _eval;
  
  /// a vector of string to hold the functions
  std::string _vars;
  
  /// a string holding the function definition
  std::string _function;
    
  /// the file reader
  Reader _reader;

  /// the file writer
  Writer  _writer;
  
  /// communicator
  MPI_Comm _comm;
  
  /// rank of this processor
  CFuint _myRank;
  
  /// number of processors
  CFuint _nbProc;
  
  /// the number of layers to grow in the z direction
  CFuint                               _nbLayers;
  
  /// the current layer
  CFuint                               _iLayer;
  
  /// number of layers local to this processor
  CFuint                               _nbLayersLocal;
  
  /// starting z coordinate local to this processor
  CFreal                               _zStartLocal;
  
  /// the size to grow in the z direction
  CFreal                                _zSize;
  
  /// split to tetrahedra?
  bool                                _split;
  
  /// random position of nodes in inner layers?
  bool                                _random;
  
  /// create one TRS called Periodic
  bool                                _periodic;
  
  /// the z diference per layer
  CFreal                                _zDelta;

  /// number of elements per layer
  CFuint                               _nbElemPerLayer;

  /// number of nodes per layer
  CFuint                               _nbNodesPerLayer;

  /// number of states per layer
  CFuint                               _nbStatesPerLayer;

  /// connectivity pattern per layer
  std::valarray<CFuint>                _nodePattern;

  /// connectivity pattern per layer
  std::valarray<CFuint>                _statePattern;

  /// id of the type TETRA
  CFuint _tetraTypeID;

  /// id of the type PRISM
  CFuint _prismTypeID;

  /// id of the type HEXA
  CFuint _hexaTypeID;

  /// store the old element-node connectivity for
  /// constructing Top and Bottom TRS's
  Common::Table<CFuint>  _oldElemNode;

  /// store the old element-state connectivity for
  /// constructing Top and Bottom TRS's
  Common::Table<CFuint>  _oldElemState;
  
}; // end class Extruder2DFVMMPI

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshExtruder

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshExtruder_Extruder2DFVMMPI_hh
