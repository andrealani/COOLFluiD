// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FAST2CFmesh_FAST2CFmeshConverter_hh
#define COOLFluiD_FAST2CFmesh_FAST2CFmeshConverter_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Framework/MeshFormatConverter.hh"
#include "Common/Table.hh"
#include "Common/CFMultiMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FAST2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * Provides an abstract interface for the format converters of the mesh files
 * in input.
 *
 * @author Andrea Lani
 *
 */
class FAST2CFmeshConverter : public Framework::MeshFormatConverter {

public: 

  /**
   * Constructor
   */
  FAST2CFmeshConverter(const std::string& name);

  /**
   * Destructor
   */
  ~FAST2CFmeshConverter();
  
  /**
   * Tries to check the file for conformity to the format.
   * Possibly not full proof.
   */
  void checkFormat(const boost::filesystem::path& filepath)
  {
  }

  /**
   * Writes the data read to the original format.
   * Useful for debugging purposes.
   */
  void convertBack(const boost::filesystem::path& filepath)
  {
  }

 /**
   * Adjust the node(= state) numbering to make it stick to the
   * COOLFluiD convention.
   */
  void adjustToCFmeshNodeNumbering() {}

protected: // functions

  /**
   * Gets the dimension.
   */
  CFuint getDimension() const {return m_dimension;}

  /**
   * Gets the number of variables
   */
  CFuint getNbVariables() const {return m_nbVars;}
  
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
    return "fgrid";
  }
  
  /**
   * Reads the all Dpl files and assembles the data in the converter.
   * It is always in a convert or convertBack, but only exectuted once.
   *
   * @throw Common::FilesystemException if a file cannot be open
   * @throw BadFormatException if a file is ill formated
   */
  void readFiles(const boost::filesystem::path& filepath);

private: // helper functions

  /**
   * Writes the FAST file when converting back to FAST
   * Useful for debugging purposes.
   * @see convertBack()
   */
  void writeFAST(const boost::filesystem::path& filepath);

  /**
   * Reads the FAST file
   *
   * @pre the extension of the mesh file is ".thor"
   * @param filepath the full file path
   * @throw Common::FilesystemException if the file cannot be open
  */
  void readFAST(const boost::filesystem::path& filepath);

  /**
   * Write in the COOLFluiD format the element list for a FEM mesh
   */
  void writeContinuousElements(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the state list for a FEM mesh
   */
  void writeContinuousStates(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the Topological Region Set data
   * considering to have FEM
   */
  void writeContinuousTrsData(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the mesh data considering
   * to have cell center FVM
   */
  void writeDiscontinuousElements(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the Topological Region Set data
   * considering to have FVM
   */
  void writeDiscontinuousTrsData(std::ofstream& fout);

  /**
   * Write in the COOLFluiD format the state list for a
   * cell centered FVM mesh
   */
  void writeDiscontinuousStates(std::ofstream& fout);
  
  /**
   * Write in the COOLFluiD format the node list
   */
  void writeNodes(std::ofstream& fout);

private: // data

  /// dimension of the mesh
  CFuint  m_dimension;
  
  /// nb of variables
  CFuint  m_nbVars;
  
  /// nb of cells in thor
  CFuint m_nbCells;
  
  /// nb of updatable nodes in thor
  CFuint m_nbUpdatableNodes;
  
  /// nb of faces in thor
  CFuint m_nbFaces;
  
  /// nb of patches in thor
  CFuint  m_nbPatches;
  
  // nodes coordinates
  Common::Table<CFreal>*            m_coordinate;
  
  // face connectivity
  Common::Table<CFuint>*            m_fcon;
    
  // boundary ID
  std::vector<int>              m_bcID;
  
  // map bcIDs to boundary face IDs 
  Common::CFMultiMap<int,CFuint>  m_mapBcToFaceIDs; 
  
  // element connectivity
  Common::Table<CFuint>*            m_econ;
  
}; // end class FAST2CFmeshConverter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FAST2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FAST2CFmesh_FAST2CFmeshConverter_hh
