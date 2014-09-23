// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshTools_ConvertStructMesh_hh
#define COOLFluiD_CFmeshTools_ConvertStructMesh_hh

//////////////////////////////////////////////////////////////////////////////

#include <set>

#include "Environment/ConcreteProvider.hh"
#include "Common/OwnedObject.hh"
#include "Common/StringOps.hh"
#include "NodeData.hh"
#include "FaceData.hh"
#include "Common/CFMultiMap.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class converts a structured mesh to unstructured Tecplot
 *
 * @author Andrea Lani
 *
 */
class ConvertStructMesh : public Common::OwnedObject{
public:
  
  typedef Environment::ConcreteProvider<ConvertStructMesh,1> PROVIDER;
  typedef const std::string& ARG1;
    
  /**
   * Constructor
   */
  ConvertStructMesh(const std::string& fileName);

  /**
   * Default destructor
   */
  virtual ~ConvertStructMesh();

  /**
   * Convert the mesh
   */
  void convert(const std::string& meshType);
  
  /**
   * Get the class name
   */
  static std::string getClassName() 
  {
    return "ConvertStructMesh";
  }
  
protected: // functions
  
  /**
   * Read topology file
   */
  virtual void readTopologyFile() = 0;
  
  /**
   * Read and build nodes with associated unique IDs
   */
  virtual void readAndBuildNodes() = 0;
  
  /**
   * Build elements
   */
  virtual void buildElements() = 0;
  
  /**
    * Write the Tecplot file
    */
  virtual void writeTecplotFile();
  
  /**
   * Write the Tecplot file for the boundary surfaces
   */
  virtual void writeBoundaryTecplotFile() = 0;
  
  /**
   * Write the CFmesh file
   */
  virtual void writeCFmeshFileFEM();
  
  /**
   * Write the CFmesh file split in triangles
   */
  virtual void writeCFmeshFileFEMSplit();
  
  /**
   * Write the CFmesh file
   */
  virtual void writeCFmeshFileFVM();
  
  /**
   * Write the CFmesh file split in triangles
   */
  virtual void writeCFmeshFileFVMSplit();
  
  /**
   * Get the global mesh file name
   */
  virtual std::string getGlobalMeshFile() const = 0;

protected: // data

  /// file name
  std::string _fileName;
  
  /// total number of blocks
  CFuint _nbBlocks;
  
  /// total number of elements
  CFuint _nbElements;
  
  /// max i in x direction
  CFuint _imax;
  
  /// max j in y direction
  CFuint _jmax; 
  
  /// max k in z direction
  CFuint _kmax;
  
  /// blocks topology information
  std::vector<Block> _blocks;
  
  /// node storage
  NodeData _nodes;
  
  /// element connectivity
  MathTools::CFMat<CFuint> _elements;
  
  /// face data storage
  FaceData _bface;
  
  /// map the BC index to the corresponding faces
  Common::CFMultiMap<int,FaceData::Itr> _mapBcIDToFace;
  
  /// array storing the unique BC IDs
  std::vector<int> _bcListIDs;
  
  /// temporary array for element IDs
  std::vector<int> _elemIDs;
  
}; // end of class ConvertStructMesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshTools_ConvertStructMesh_hh
