// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshTools_ConvertQuadMesh_hh
#define COOLFluiD_CFmeshTools_ConvertQuadMesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConvertStructMesh/ConvertStructMesh.hh"

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
class ConvertQuadMesh : public ConvertStructMesh {
public:

  /**
   * Constructor
   */
  ConvertQuadMesh(const std::string& fileName);

  /**
   * Default destructor
   */
  ~ConvertQuadMesh();

protected: // functions

  /**
   * Read topology file
   */
  void readTopologyFile() {}

  /**
   * Read and build nodes with associated unique IDs
   */
  void readAndBuildNodes();

  /**
   * Build elements
   */
  void buildElements();

  /**
   * Write the Tecplot file for the boundary surfaces
   */
  void writeBoundaryTecplotFile() {}

  /**
   * Add the boundary faces in @see FaceData
   */
  void addBoundaryFaces(CFuint iElem, CFuint i, CFuint j);

  /**
   * Get the global mesh file name
   */
  std::string getGlobalMeshFile() const
  {
    return _fileName + ".dat";
  }

}; // end of class ConvertQuadMesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshTools_ConvertQuadMesh_hh
