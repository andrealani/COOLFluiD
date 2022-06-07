// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_SetupExtra_hh
#define COOLFluiD_FluxReconstructionMethod_SetupExtra_hh

//////////////////////////////////////////////////////////////////////////////


#include "FluxReconstructionMethod/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a command to setup the FluxReconstruction method with extra sockets
 * 
 * @author Ray Vandenhoeck
 */
class SetupExtra : public StdSetup {

public: // functions

  /// Constructor
  explicit SetupExtra(const std::string& name);

  /// Destructor
  ~SetupExtra();

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
  
  /// storage for the cell volumes
  Framework::DataSocketSource< CFreal > socket_cellVolumes;
  
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
  
};  // class SetupExtra

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_SetupExtra_hh


