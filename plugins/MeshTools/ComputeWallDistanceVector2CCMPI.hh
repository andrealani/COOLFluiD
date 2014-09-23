// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ComputeWallDistanceVector2CCMPI_hh
#define COOLFluiD_MeshTools_ComputeWallDistanceVector2CCMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "MeshTools/ComputeWallDistance.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
      
  namespace Framework {
    class TRSDistributeData;
  }
  
  namespace MeshTools {
      
//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the distance from the states to the wall
 * and outputs to a file
 *
 * @author Andrea Lani
 * @author Thomas Wuilbaut
 *
 */
class ComputeWallDistanceVector2CCMPI : public ComputeWallDistance {
public:

  /**
   * Constructor.
   */
  ComputeWallDistanceVector2CCMPI(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeWallDistanceVector2CCMPI();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
  
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
  
  /**
   * Execute on a set of dofs
   */
  void execute2D();

  /**
   * Execute on a set of dofs
   */
  void execute3D();
  
  /**
   * Compute the wall distance (more complex but more accurate)
   */
  void computeWallDistance2D(Framework::TRSDistributeData& data);
  
  /**
   * Compute the wall distance
   */
  void computeWallDistance3D(std::vector<CFreal>& data);
  
private:
  
  /// temporary coefficient
  CFreal m_t;

  /// temporary distance between the inner state and the wall
  CFreal m_drXiXw;

  /// temporary distance between the inner state and the repositioned ghost state
  CFreal m_drXiXg;

  /// temporary internal node
  RealVector* m_innerNode;

  /// temporary node
  RealVector m_tempNode;

  /// temporary middle node
  RealVector m_midNode;

  /// temporary ghost node
  RealVector m_tempGhostNode;

  /// temporary face normal
  RealVector m_faceNormal;
  
#ifdef CF_HAVE_MPI
  /// communicator
  MPI_Comm m_comm;
#endif
  
  /// rank of this processor
  CFuint m_myRank;

  /// number for processors
  CFuint m_nbProc;
  
}; // end of class ComputeWallDistanceVector2CCMPI

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ComputeWallDistanceVector2CCMPI_hh
