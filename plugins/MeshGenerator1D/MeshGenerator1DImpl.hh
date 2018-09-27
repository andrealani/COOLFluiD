// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshGeneratorID_MeshGenerator1DImpl_hh
#define COOLFluiD_MeshGeneratorID_MeshGenerator1DImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <vector>

#include "Framework/MeshFormatConverter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshGenerator1D {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * Provides an abstract interface for the format converters of the mesh files
 * in input.
 *
 * @author Alessandro Munafo'
 *
 */
class MeshGenerator1DImpl : public Framework::MeshFormatConverter {
public:                 

  
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );


  /**
   * Constructor
   */
  MeshGenerator1DImpl(const std::string& name);
  
  /**
   * Destructor
   */
  ~MeshGenerator1DImpl();
  
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
  
protected:

  /**
   * Adjust the node (state) numbering to make it stick to the
   * COOLFluiD convention.
   */
  void adjustToCFmeshNodeNumbering() {}
  
  /**
   * Gets the dimension.
   */
  CFuint getDimension() const {return 1;}
  
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
    return "dat";
  }
  
  /**
   * Reads the all Dpl files and assembles the data in the converter.
   * It is always in a convert or convertBack, but only exectuted once.
   *
   * @throw Common::FilesystemException if a file cannot be open
   * @throw BadFormatException if a file is ill formated
   */
  void readFiles(const boost::filesystem::path& filepath);
  
private:

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
  
  /**
   * Read the radius info
   */
  void nozzleRadius();

  /**
   * Radius file for  unsteady shock tube and steady shock computations
   */ 
  void shockRadius();  

private:
  
  // Number of nodes 
  CFuint  _nbUpdatableNodes;

  // Number of cells
  CFuint  _nbCells; 

  // Vector of node coordinate                           
  std::vector <CFreal> _coordinate; 

  /// Configuration of boundary
  bool m_isPeriodicBoundary;
  
  /// Array of size 3 telling start position, end position, number of cells
  std::vector<CFreal> m_startEndN;
  
}; // end class MeshGenerator1DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshGenerator1D

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshGenerator1D_MeshGenerator1DImpl_hh
