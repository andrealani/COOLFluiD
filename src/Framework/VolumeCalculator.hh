// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_VolumeCalculator_hh
#define COOLFluiD_Framework_VolumeCalculator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "NegativeVolumeException.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Node;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a calculator for volumes of elements
/// For 2D elements it calculates areas.
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
/// @author Andrea Lani

class Framework_API VolumeCalculator {
public:

  /// Default constructor without arguments.
  VolumeCalculator();

  /// Default destructor.
  virtual ~VolumeCalculator();

  /// Calculate the volume of a triangle
  /// @param coord coordinates of the nodes
  CFreal calculateTriagVolume(const std::vector<Node*>& nodes);

  /// Calculate the volume of a triangle in 3D
  /// @param coord coordinates of the nodes
  CFreal calculate3DTriagVolume(const std::vector<Node*>& nodes);

  /// Calculate the volume of a quadrilater
  /// @param coord coordinates of the nodes
  CFreal calculateQuadVolume(const std::vector<Node*>& nodes);

  /// Calculate the volume of a tetrahedra
  /// with coordinates given in a matrix nodes x dimension, (4x3).
  /// @param coord coordinates of the nodes
  CFreal calculateTetraVolume(const RealMatrix& coord);


  /// Calculate the volume of a tetrahedra
  /// with coordinates given in a vector of nodes*
  /// @param coord coordinates of the nodes
  CFreal calculateTetraVolume(const std::vector<Node*>& coord);

  /// Calculate the volume of a pyramid
  /// with coordinates given in a matrix nodes x dimension, (5x3).
  /// Equivalent to computation of 2 tetrahedra volumes.
  /// @param coord coordinates of the nodes
  CFreal calculatePyramVolume(const RealMatrix& coord);

  /// Calculate the volume of a pyramid
  /// with coordinates given in a vector of nodes*
  /// @param coord coordinates of the nodes
  CFreal calculatePyramVolume(const std::vector<Node*>& coord);

  /// Calculate the volume of a prism
  /// with coordinates given in a matrix nodes x dimension, (6x3).
  /// Equivalent to computation of 3 tetrahedra volumes.
  /// @param coord coordinates of the nodes
  CFreal calculatePrismVolume(const RealMatrix& coord);

  /// Calculate the volume of a prism
  /// with coordinates given in a vector of nodes*
  /// @param coord coordinates of the nodes
  CFreal calculatePrismVolume(const std::vector<Node*>& coord);

  /// Calculate the volume of a hexahedra
  /// with coordinates given in a matrix nodes x dimension, (8x3).
  /// Equivalent to computation of 2 prism volumes or
  /// 6 tetrahedra volumes.
  /// @param coord coordinates of the nodes
  CFreal calculateHexaVolume(const RealMatrix& coord);

  /// Calculate the volume of a hexa
  /// with coordinates given in a vector of nodes*
  /// @param coord coordinates of the nodes
  CFreal calculateHexaVolume(const std::vector<Node*>& coord);

  /// Check Node Numbering by computing the volume of a triangle
  /// @param coord coordinates of the nodes
  bool checkTriagNumbering(const std::vector<Node*>& nodes);

  /// Check Node Numbering by computing the volume of a quadrilater
  /// @param coord coordinates of the nodes
  bool checkQuadNumbering(const std::vector<Node*>& nodes);

  /// Check Node Numbering by computing the volume of a tetrahedra
  /// with coordinates given in a matrix nodes x dimension, (4x3).
  /// @param coord coordinates of the nodes
  bool checkTetraNumbering(const RealMatrix& coord);

  /// Check Node Numbering by computing the volume of a pyramid
  /// with coordinates given in a matrix nodes x dimension, (5x3).
  /// Equivalent to computation of 2 tetrahedra volumes.
  /// @param coord coordinates of the nodes
  bool checkPyramNumbering(const RealMatrix& coord);

  /// Check Node Numbering by computing the volume of a prism
  /// with coordinates given in a matrix nodes x dimension, (6x3).
  /// Equivalent to computation of 3 tetrahedra volumes.
  /// @param coord coordinates of the nodes
  bool checkPrismNumbering(const RealMatrix& coord);

  /// Check Node Numbering by computing the volume of a hexahedra
  /// with coordinates given in a matrix nodes x dimension, (8x3).
  /// Equivalent to computation of 2 prism volumes or
  /// 6 tetrahedra volumes.
  /// @param coord coordinates of the nodes
  bool checkHexaNumbering(const RealMatrix& coord);

protected:

  /// Helper function to compute the determinant of
  /// the tetra already set in the temporary matrix.
  CFreal calcVolume() const {
    return -(1./6.) * _detMat.determ4();
  } 

  /// Helper function to set the coordinates into
  /// the temporary matrix.
  /// This function is inline to ease compiler
  /// optimization.
  void setTetraIDs(const RealMatrix& coord, 
		   const CFint& i, const CFint& j,
                   const CFint& k, const CFint& l)
  {
    _detMat(0,0) = coord(i,0);
    _detMat(0,1) = coord(i,1);
    _detMat(0,2) = coord(i,2);

    _detMat(1,0) = coord(j,0);
    _detMat(1,1) = coord(j,1);
    _detMat(1,2) = coord(j,2);

    _detMat(2,0) = coord(k,0);
    _detMat(2,1) = coord(k,1);
    _detMat(2,2) = coord(k,2);

    _detMat(3,0) = coord(l,0);
    _detMat(3,1) = coord(l,1);
    _detMat(3,2) = coord(l,2);
  }

  /// Helper function to set the coordinates into
  /// the temporary matrix in case the input is not
  /// a matrix but a vector of Node*.
  /// This function is inline to ease compiler
  /// optimization.
  void setTetraIDs(const std::vector<Node*>& coord,
                   const CFint& i,
                   const CFint& j,
                   const CFint& k,
                   const CFint& l)
  {
    _detMat(0,0) = (*(coord[i]))[0];
    _detMat(0,1) = (*(coord[i]))[1];
    _detMat(0,2) = (*(coord[i]))[2];

    _detMat(1,0) = (*(coord[j]))[0];
    _detMat(1,1) = (*(coord[j]))[1];
    _detMat(1,2) = (*(coord[j]))[2];
    
    _detMat(2,0) = (*(coord[k]))[0];
    _detMat(2,1) = (*(coord[k]))[1];
    _detMat(2,2) = (*(coord[k]))[2];

    _detMat(3,0) = (*(coord[l]))[0];
    _detMat(3,1) = (*(coord[l]))[1];
    _detMat(3,2) = (*(coord[l]))[2];
  }
  
  /// hexa volume by subdivision into 6 tetra
  CFreal calculateHexaVolumeBy6Tetra(const std::vector<Node*>& coord);
  
  /// hexa volume by subdivision into 6 pyramids centered in the centroid
  CFreal calculateHexaVolumeBy6Pyram(const std::vector<Node*>& coord);
  
  /// hexa volume by long diagonal fast algorithm from J. Grandy (1997)
  CFreal calculateHexaVolumeByLD(const std::vector<Node*>& coord);
  
  /// hexa volume by tetrakis hexahedron fast and accurate algorithm from J. Grandy (1997)
  CFreal calculateHexaVolumeByTH(const std::vector<Node*>& coord);
  
private: // data

  /// temporary storage of matrix for determinant calculation
  /// hence this matrix has one more dimension
  RealMatrix _detMat;

  /// baricenter
  RealVector _center;

}; // end of class VolumeCalculator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_VolumeCalculator_hh
