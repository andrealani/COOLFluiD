// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ComputeWallDistanceNewton_hh
#define COOLFluiD_MeshTools_ComputeWallDistanceNewton_hh

//////////////////////////////////////////////////////////////////////////////

#include "MeshTools/ComputeWallDistance.hh"
#include "MeshTools/ComputeShortestDistance.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the distance from the states to the wall
 * and outputs to a file
 *
 * @author Thomas Wuilbaut
 *
 */
class ComputeWallDistanceNewton : public ComputeWallDistance {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  ComputeWallDistanceNewton(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeWallDistanceNewton();

  /**
   * Execute on a set of dofs
   */
  void execute();

protected: //data

  /// Newton method for computing the shortest distance to a face
  ComputeShortestDistance _shortestDistanceComputer;

  ///Tolerance for the stop condition of the computation of the shortest distance
  CFreal _tolerance;

  ///Relaxation actor for the newton procedure
  CFreal _relaxFactor;

  ///Maximum number of iterations for the newton procedure
  CFuint _maxIter;

}; // end of class ComputeWallDistanceNewton

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ComputeWallDistanceNewton_hh
