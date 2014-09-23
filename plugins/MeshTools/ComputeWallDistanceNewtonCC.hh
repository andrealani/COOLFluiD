// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ComputeWallDistanceNewtonCC_hh
#define COOLFluiD_MeshTools_ComputeWallDistanceNewtonCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "MeshTools/ComputeWallDistanceNewton.hh"

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
class ComputeWallDistanceNewtonCC : public ComputeWallDistanceNewton {
public:

  /**
   * Constructor.
   */
  ComputeWallDistanceNewtonCC(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeWallDistanceNewtonCC();

  /**
   * Execute on a set of dofs
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for State's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

}; // end of class ComputeWallDistanceNewtonCC

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ComputeWallDistanceNewtonCC_hh
