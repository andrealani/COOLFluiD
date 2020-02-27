// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ComputeWallDistanceVectorFRMPI_hh
#define COOLFluiD_MeshTools_ComputeWallDistanceVectorFRMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/TRSDistributeData.hh"
#include "MeshTools/ComputeWallDistance.hh"
#include <vector>
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
    namespace FluxReconstructionMethod {
    class FluxReconstructionSolverData;
  }
  
  namespace MeshTools {
    
//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the distance from the states to the wall for FR
 * and outputs to a file 
 *
 * @author Ray Vandenhoeck
 * @author Andrea Lani
 * @author Thomas Wuilbaut
 *
 */
class ComputeWallDistanceVectorFRMPI : public ComputeWallDistance {
public:



  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor.
   */
  ComputeWallDistanceVectorFRMPI(const std::string& name);
  
  /**
   * Default destructor
   */
  ~ComputeWallDistanceVectorFRMPI();

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


  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > provideSockets();

    
  /**
   * This class stores face centers and normals to be distributed in parallel
   *
   * @author Andrea Lani
   */
  class TRSFaceDistributeData : public Framework::TRSDistributeData {
  public:
    std::vector<CFreal> faceCenters; /// face centers
    std::vector<CFreal> faceNormals; /// face normals
  };
  
  
  void execute3D();
  
  /**
   * Compute the wall distance
   */
  void computeWallDistance3D(std::vector<CFreal>& data);
  
  /**
   * Compute the wall distance (2D and 3D)
   */
  void computeWallDistance(TRSFaceDistributeData& data);
   
  void computeWallDistanceExtrapolate(TRSFaceDistributeData& data);   
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
  
  /// minimum state-face distance
  RealVector m_minStateFaceDistance;
  
#ifdef CF_HAVE_MPI
  /// communicator
  MPI_Comm m_comm;
#endif
  
  /// rank of this processor
  CFuint m_myRank;

  /// number for processors
  CFuint m_nbProc;

  /// Define the acceptable distance  
 
  CFreal m_acceptableDistance;
  
  // pointer to the data of the cell centered FVM method
  Common::SafePtr<FluxReconstructionMethod::FluxReconstructionSolverData> m_frData;
  
  /// socket for flux point normals
  Framework::DataSocketSink< CFreal > socket_flxPntNormals;
  
  }; // end of class ComputeWallDistanceVectorFRMPI

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ComputeWallDistanceVectorFRMPI_hh
