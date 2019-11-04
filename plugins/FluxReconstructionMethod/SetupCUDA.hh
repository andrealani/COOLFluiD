// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh
#define COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh

//////////////////////////////////////////////////////////////////////////////


#include "FluxReconstructionMethod/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a command to setup the FluxReconstruction method if porting to GPU is needed
 * 
 * @author Ray Vandenhoeck
 */
class SetupCUDA : public StdSetup {

public: // functions

  /// Constructor
  explicit SetupCUDA(const std::string& name);

  /// Destructor
  ~SetupCUDA();

  /// Execute processing actions
  void execute();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  
private: // private functions    

  /**
   * Creates sockets that store the needed normals
   */
  void createNormalSockets();
  
protected: // data
  
  /// storage for the normals in the solution points
  Framework::DataSocketSource< CFreal > socket_solPntNormals;
  
  /// storage for the normals in the solution points
  Framework::DataSocketSource< CFreal > socket_flxPntNormals;

  /// storage for the face directions
  Framework::DataSocketSource< CFint > socket_faceDir;
  
  /// variable for current face
  Framework::GeometricEntity* m_face;
  
  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;
  
  /// face local coordinates of the flux points on one face
  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoords;
  
  /// vector to store the face jacobians in
  std::vector< RealVector > m_faceJacobVecs;

  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< std::vector< CFint > > > m_faceMappedCoordDir;

  /// flx pnt - face connectivity per orient
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_faceFlxPntConnPerOrient;
  
};  // class SetupCUDA

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh


