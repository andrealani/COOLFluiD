// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_StdSetup_hh
#define COOLFluiD_FluxReconstructionMethod_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MatrixInverterT.hh"

#include "Framework/Node.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to setup the FluxReconstruction method
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
class StdSetup : public FluxReconstructionSolverCom {

public: // functions

  /// Constructor
  explicit StdSetup(const std::string& name);

  /// Destructor
  ~StdSetup();

  /// Execute processing actions
  void execute();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();
   
private: // private functions    

  /**
   * Creates socket with the cell volume for each state
   */
  void computeStatesVolumes();
    
   /**
   * Gets the start and stop indexes of the range of faces with a certain orientation from the connectivity
   * where it is stored and puts it in a socket.
   */
  void createFaceOrientationStartIndexes();

  /**
   * Computes the unit normals and the the size of the projection vectors in face flux points.
   */
  void computeFaceJacobianVectorSizes();

  /**
   * Set the boundary condition type of the boundary faces
   */
  void setBndFacesBCType();
  
  /**
   * Set the coordinates of the states
   */
  void setStateCoords();

  /**
   * Sets the inverse jacobian matrix multiplied with the jacobian determinant at the given mapped coordinates
   */
  void setInvGeoJacobMatrixXJacobDet(std::vector< RealMatrix >& matr,
                                     Framework::GeometricEntity& cell,
                                     const std::vector< RealVector >& mappedCoord);
  
protected: // data
  
  /// storage for interpolated values in the mesh vertices
  Framework::DataSocketSource<RealVector> socket_nstates;

  /// socket with proxy to be able to use the data handle of nodal states
  /// uniformly independently from the actual storage type being RealVector or
  /// State*
  Framework::DataSocketSource< Framework::ProxyDofIterator< RealVector >* >
    socket_nstatesProxy;

  /// socket for state's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// socket with flags to check if a state has been updated or not
  Framework::DataSocketSource<bool> socket_isUpdated;
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for cell volumes
  Framework::DataSocketSource< CFreal > socket_volumes;
  
  /// socket for gradients
  Framework::DataSocketSource< std::vector< RealVector > > socket_gradients;
  
  /// socket for the gradients needed for artificial viscosity
  Framework::DataSocketSource< std::vector< RealVector > > socket_gradientsAV;
  
  /// socket for the values needed for the positivity preservation
  Framework::DataSocketSource< CFreal > socket_posPrev;
  
  /// socket for size of face normal jacobian in face flux points
  Framework::DataSocketSource< std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;
  
  /// socket for normals
  Framework::DataSocketSource< CFreal > socket_normals;

  /// array with IDs
  std::vector<CFuint> m_nodeIDToStateID;
  
  /// matrix inverter size 2
  MathTools::MatrixInverterT<2> inverter2;
  /// matrix inverter size 3
  MathTools::MatrixInverterT<3> inverter3;
  
  /// nodes to be set as coordinates for the states
  std::vector< Framework::Node > m_stateNodes;
  
};  // class StdSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdSetup_hh

