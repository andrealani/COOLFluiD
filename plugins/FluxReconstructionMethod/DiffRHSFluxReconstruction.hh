// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the (diffusive part of the) system using a FluxReconstruction solver
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class DiffRHSFluxReconstruction : public FluxReconstructionSolverCom {

public: // functions

  /// Constructor
  explicit DiffRHSFluxReconstruction(const std::string& name);

  /// Destructor
  virtual ~DiffRHSFluxReconstruction() {}

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

  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /// compute the divergence of the discontinuous flux (-divFD+divhFD)
  virtual void computeDivDiscontFlx(std::vector< RealVector >& residuals);
  
  /// add the residual updates to the RHS
  void updateRHS();
  
  /// add the updates to the wave speed
  void updateWaveSpeed();
  
  /// compute the correction -(FI)divh of a neighbouring cell
  void computeCorrection(CFuint side, std::vector< RealVector >& corrections);
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /**
   * Divides by jacobian determinant
   */
  void divideByJacobDet();
  
  /**
   * Set the data for the current face necessary to calculate FI
   */
  virtual void setFaceData(CFuint faceID);
  
  /**
   * Set the data for the current cell necessary to calculate the residual update
   */
  virtual void setCellData();
  
  /**
   * Compute the left and right states and gradients in the flx pnts
   */
  void computeFlxPntStatesAndGrads();
  
  virtual void computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux);
  
  /// prepare the computation of the diffusive flux
  virtual void prepareFluxComputation()
  {
  }

protected: //data
  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// socket for gradientsAV
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradientsAV;
  
  /// socket for positivity preservation values
  Framework::DataSocketSink< CFreal > socket_posPrev;
  
  /// storage of the rhs
  Framework::DataSocketSink< CFreal > socket_rhs;
  
  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink< CFreal > socket_updateCoeff;
  
  /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;
  
  /// diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;
  
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
  
  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;
  
  /// number of equations in the physical model
  CFuint m_nbrEqs;
  
  /// number of dimensions in the physical model
  CFuint m_dim;
  
  /// variable for current face orientation
  CFuint m_orient;
  
  /// number of solution pnts in the cell
  CFuint m_nbrSolPnts;;
  
  /// number of face flx pnts
  CFuint m_nbrFaceFlxPnts;
  
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
  
  /// flux point coordinates
  std::vector< RealVector > m_flxPntCoords;
  
  /// flux projection vectors in solution points for disc flux
  std::vector< std::vector< RealVector > > m_cellFluxProjVects;
  
  /// the gradients in the neighbouring cell
  std::vector< std::vector< std::vector< RealVector >* > > m_cellGrads;
  
  /// the corrected gradients in the flux points
  std::vector< std::vector< std::vector< RealVector* > > > m_cellGradFlxPnt;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  /// coefs to compute the derivative of the states in the sol pnts
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPnts;
  
  /// face inverse characteristic lengths
  std::vector< CFreal > m_faceInvCharLengths;
  
  /// cell volume
  std::vector< CFreal > m_cellVolume;
  
  /// ratio between convective and diffusive cfl limit
  CFreal m_cflConvDiffRatio;
  
  /// local cell face - flux point cell mapped coordinate per face connection orientation
  Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > m_faceFlxPntCellMappedCoords;
  
  /// flag telling whether to freeze the gradients
  bool m_freezeGrads;
  
  /// the discontinuous flux extrapolated to the flux points
  std::vector< RealVector > m_extrapolatedFluxes;
  
  /// dimensions on which to evaluate the flux in the flux points
  Common::SafePtr< std::vector< CFuint > >  m_flxPntFlxDim;
  
  /// average solution in a flux point
  RealVector m_avgSol;
  
  /// average gradients in a flux point
  std::vector< RealVector* > m_avgGrad;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstruction_hh

