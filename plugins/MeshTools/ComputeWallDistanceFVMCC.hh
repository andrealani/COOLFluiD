// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ComputeWallDistanceFVMCC_hh
#define COOLFluiD_MeshTools_ComputeWallDistanceFVMCC_hh

//////////////////////////////////////////////////////////////////////

#include "ComputeWallDistance.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the distance from the states to the wall
 * and outputs to a file
 *
 * @author Thomas Wuilbaut
 *
 */
class ComputeWallDistanceFVMCC : public ComputeWallDistance {
public:

  /**
   * Constructor.
   */
  ComputeWallDistanceFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeWallDistanceFVMCC();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute on a set of dofs
   */
  void execute();


protected: //data

  /// socket for State's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool
		 <Framework::CellTrsGeoBuilder> > _cellBuilder;

}; // end of class ComputeWallDistanceFVMCC

//////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ComputeWallDistanceFVMCC_hh
