// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_StdSolve_hh
#define COOLFluiD_FluxReconstructionMethod_StdSolve_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/ReconstructStatesFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the system using FluxReconstruction solver
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class StdSolve : public FluxReconstructionSolverCom {

public: // functions

  /// Constructor
  explicit StdSolve(const std::string& name);

  /// Destructor
  virtual ~StdSolve() {}

  /// Execute processing actions
  void execute();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: //data
  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// socket for normals
  Framework::DataSocketSink< CFreal > socket_normals;
  
  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;
  
  /// index of element type
  CFuint m_iElemType;
  
  /// variable for cell
  Framework::GeometricEntity* m_cell;
  
  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;
  
  /// vector containing pointers to the states in a cell
  std::vector< std::vector< Framework::State* >* > m_cellStatesFlxPnt;
  
  /// vector containing pointers to the internal fluxes in a cell (fr each sol pnt)
  std::vector< RealVector >* m_cellIntFlx;
  
  /// vector containing pointers to the fluxes in the flux points
  std::vector< Common::SafePtr< std::vector< RealVector > > > m_cellFlx;
  
  /// vector containing pointers to the face normals
  Common::SafePtr< std::vector< RealVector > > m_faceNormals;
  
  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;
  
  /// flux point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_flxPntsLocalCoords;
  
  /// flx pnt - face connectivity
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_faceFlxPntConn;
  
  /// face connectivity per orient
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceConnPerOrient;
  
  /// interface minus discontinuous flux in the flux points for each face
  std::vector< Common::SafePtr< std::vector< std::vector< RealVector > > > > m_corrFlxFactor;
  
  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;
  
  /// number of equations in the physical model
  CFuint m_nbrEqs;
  
  /// number of dimensions in the physical model
  CFuint m_dim;
  
  /// variable for current face orientation
  CFuint m_orient;
  
  /// variable for current face
  Framework::GeometricEntity* m_face;
  
  /// variable for current neighbouring cells
  std::vector< Framework::GeometricEntity* > m_cells;
  
  /// Riemann flux
  Common::SafePtr< RiemannFlux > m_riemannFluxComputer;
  
  /// variable for the states in the left and right cell
  std::vector< std::vector< Framework::State* >* > m_states;
  
  /// Interface fluxes at the flux points of a face
  std::vector< RealVector> m_flxPntRiemannFlux;
  
  /// Continuous flux at the solution points
  std::vector< RealMatrix> m_contFlx;
  
  /// Divergence of the continuous flux at the solution points
  std::vector< RealVector> m_divContFlx;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdSolve_hh

