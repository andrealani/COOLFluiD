// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ComputeShortestDistance_hh
#define COOLFluiD_MeshTools_ComputeShortestDistance_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntity.hh"
#include "Framework/NumericalJacobian.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the shortest distance from a point to a GeometricEntity
 *
 * @author Thomas Wuilbaut
 *
 */

class ComputeShortestDistance {
public:

  /**
   * Constructor.
   */
  ComputeShortestDistance();

  /**
   * Default destructor
   */
  ~ComputeShortestDistance();


  /**
   * Compute the shortest distance to the face
   */
  void compute(Framework::GeometricEntity* face, RealVector& coord, CFreal& minimumDistance);

  /**
   * Resize
   */
  void setDimension(const CFuint dim);

  /**
   * Set a relaxation factor
   */
  void setRelaxationFactor(const CFreal relax)
  {
    _relaxation = relax;
  }

  /**
   * Set the tolerance for stopping
   */
  void setTolerance(const CFreal tolerance)
  {
    _tolerance = tolerance;
  }

  /**
   * Set the maximum number of iterations
   */
  void setMaxIter(const CFuint maxIter)
  {
    _maxIter = maxIter;
  }

private: //function


  /**
   * Newton Loop: Setup and Initialize the vectors
   */
  void setupNewton();

  /**
   * Newton Loop: Computes the residual and the jacobian
   */
  void takeStep();

  /**
   * Newton Loop: Updates the Solution
   */
  void updateSolution();

  /**
   * Newton Loop: Computes the Residual
   */
  CFreal computeResidual();

private: //data

  /// Numerical jacobian calculator
  std::auto_ptr<Framework::NumericalJacobian> _numericalJacob;

  /// Pointer to the current face
  Framework::GeometricEntity * _currentFace;

  /// Mapped Coordinates
  RealVector _mappedCoord;

  /// Coordinates of the point to be projected
  RealVector _coord;

  /// Residual
  RealVector _residual;

  /// Perturbed Residual
  RealVector _otherResidual;

  /// Jacobian
  RealVector _jacobian;

  /// Relaxation Factor
  CFreal _relaxation;

  ///Tolerance for the stop condition
  CFreal _tolerance;

  ///Maximum number of iterationsfor the stop condition
  CFuint _maxIter;

  /// isResized
  bool _isSetup;

}; // end of class ComputeShortestDistance

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ComputeShortestDistance_hh
