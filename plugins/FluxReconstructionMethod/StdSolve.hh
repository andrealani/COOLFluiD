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
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/ReconstructStatesFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the system using a FluxReconstruction solver
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
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
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
    
protected: //functions

  /// compute the interface flux correction FI-FD
  void computeInterfaceFlxCorrection(CFuint faceID);
  
  /// compute the residual updates (-divFC)
  void computeResUpdates(CFuint elemIdx);
  
  /// add the residual updates to the RHS
  void updateRHS();
  
  /// add the updates to the wave speed
  void updateWaveSpeed();
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /**
   * Divides by jacobian determinant
   */
  void divideByJacobDet();

protected: //data
  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// socket for normals
  Framework::DataSocketSink< CFreal > socket_normals;
  
  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;
  
  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal > socket_updateCoeff;
  
  /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;
  
  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;
  
  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;
  
  /// index of element type
  CFuint m_iElemType;
  
  /// variable for cell
  Framework::GeometricEntity* m_cell;
  
  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;
  
  /// extrapolated states in the flux points of the cell
  std::vector< std::vector< Framework::State* > > m_cellStatesFlxPnt;
  
  /// vector containing pointers to the internal fluxes in a cell (fr each sol pnt)
  std::vector< RealVector >* m_cellIntFlx;
  
  /// vector containing pointers to the fluxes in the flux points
  std::vector< std::vector< RealVector > > m_cellFlx;
  
  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;
  
  /// flux point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_flxPntsLocalCoords;
  
  /// flx pnt - face connectivity per orient
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_faceFlxPntConnPerOrient;
  
  /// flx pnt - face connectivity
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;
  
  /// face connectivity per orient
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceConnPerOrient;
  
  /// interface minus discontinuous flux in the flux points for each face
  std::vector< std::vector< RealVector > > m_corrFlxFactor;
  
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
  
  /// Correction function computer
  Common::SafePtr< BaseCorrectionFunction > m_corrFctComputer;
  
  /// Correction function for current cell
  std::vector< std::vector< RealVector > > m_corrFct;
  
  /// Divergence of the correction function for current cell
  std::vector< std::vector< CFreal > > m_corrFctDiv;
  
  /// variable for the states in the left and right cell
  std::vector< std::vector< Framework::State* >* > m_states;
  
  /// Interface fluxes at the flux points of a face
  std::vector< RealVector> m_flxPntRiemannFlux;
  
  /// Continuous flux at the solution points
  std::vector< std::vector< RealVector> > m_contFlx;
  
  /// Divergence of the continuous flux at the solution points
  std::vector< RealVector> m_divContFlx;
  
  /// updates for the wave speed
  std::vector< CFreal > m_waveSpeedUpd;
  
  /// face Jacobian vector sizes (abs)
  std::vector< CFreal > m_faceJacobVecAbsSizeFlxPnts;
  
  /// coefficients for integration over a face
  Common::SafePtr< RealVector > m_faceIntegrationCoefs;
  
  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< std::vector< CFint > > > m_faceMappedCoordDir;
  
  /// local cell face - mapped coordinate direction
  Common::SafePtr< std::vector< CFint > > m_faceLocalDir;
  
  /// unit normal vector in flux points
  std::vector< RealVector > m_unitNormalFlxPnts;
  
  /// face Jacobian vector sizes
  std::vector< std::vector< CFreal > > m_faceJacobVecSizeFlxPnts;
  
  /// 1D flux points mapped coordinates
  Common::SafePtr< std::vector< CFreal > > m_flxPntsLocalCoords1D;
  
  /// flux point coordinates
  std::vector< RealVector > m_flxPntCoords;
  
  /// flux projection vectors in solution points for disc flux
  std::vector< std::vector< RealVector > > m_cellFluxProjVects;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdSolve_hh

