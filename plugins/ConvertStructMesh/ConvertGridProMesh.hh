// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshTools_ConvertGridProMesh_hh
#define COOLFluiD_CFmeshTools_ConvertGridProMesh_hh

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
class ConvertGridProMesh : public ConvertStructMesh {
public:

  /**
   * Constructor
   */
  ConvertGridProMesh(const std::string& fileName);

  /**
   * Default destructor
   */
  ~ConvertGridProMesh();

protected: // functions

  /**
   * Read topology file
   */
  void readTopologyFile();

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
  void writeBoundaryTecplotFile();

  /**
   * Get the number of faces in block face
   */
  CFuint getNbFacesInBlockFace(CFuint maxnx, CFuint maxny,
			       CFuint maxnz, int hexaFaceID) const;

  /**
   * Add the boundary face in @see FaceData
   */
  FaceData::Itr addBoundaryFace(CFuint iElem, CFuint i, CFuint j, CFuint k,
				int hexaFaceID, int bcID);

  /**
   * Get the global mesh file name
   */
  std::string getGlobalMeshFile() const
  {
    return _fileName + ".vol";
  }

protected: // data

}; // end of class ConvertGridProMesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshTools_ConvertGridProMesh_hh
