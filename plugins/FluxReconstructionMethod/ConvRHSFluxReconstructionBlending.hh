// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to assemble the convective part of the system using a FluxReconstruction solver with order blending
 * 
 * @author Rayan Dhib
 * 
 * based on the class ConvRHSFluxReconstruction
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
class ConvRHSFluxReconstructionBlending : public FluxReconstructionSolverCom {

public: // functions

  /// Constructor
  explicit ConvRHSFluxReconstructionBlending(const std::string& name);

  /// Destructor
  virtual ~ConvRHSFluxReconstructionBlending() {}

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
  void computeInterfaceFlxCorrection();
  
  /// compute the divergence of the discontinuous flx (-divFD+divhFD)
  void computeDivDiscontFlx(std::vector< RealVector >& residuals);
  
  /// add the residual updates to the RHS
  void updateRHS();
  
  /**
   * add the residual updates to the RHS of both neighbor cells
   */
  void updateRHSBothSides();
  
  /// add the updates to the wave speed
  void updateWaveSpeed();
  
  /// compute the correction -(FI)divh of a neighbouring cell of the face
  void computeCorrection(CFuint side, std::vector< RealVector >& corrections);
  
  /**
   * compute the wave speed updates for this face
   */
  void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd, 
                               std::vector< CFreal >& waveSpeedUpdP0);
  
  /**
   * Set the data for the current face necessary to calculate FI-FD
   */
  void setFaceData(CFuint faceID);
  
  /**
   * Set the data for the current cell necessary to calculate the residual update
   */
  void setCellData();
  
  /**
   * Compute the left and right states in the flx pnts
   */
  void computeFlxPntStates();
  
  /// compute the volume term contribution to the gradients
  virtual void computeGradients();
  
  /// compute the face correction to the corrected gradients
  virtual void computeGradientFaceCorrections();

protected: //data
  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// socket for gradientsAV
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradientsAV;
  
  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;
  
  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal > socket_updateCoeff;
  
  /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;

  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPntsP0;
  
  /// storage of the volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// socket for output of the filtering
  Framework::DataSocketSink< CFreal > socket_alpha;
  
  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;
  
  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;
  
  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;
  
  /// index of element type
  CFuint m_iElemType;
  
  /// variable for cell
  Framework::GeometricEntity* m_cell;
  
  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;
  
  /// extrapolated states in the flux points of the cell
  std::vector< std::vector< Framework::State* > > m_cellStatesFlxPnt;

  std::vector< std::vector< Framework::State* > > m_cellStatesFlxPntP0;
  
  /// vector containing pointers to the fluxes in the flux points
  std::vector< std::vector< RealVector > > m_cellFlx;

  std::vector< std::vector< RealVector > > m_cellFlxP0;
  
  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoordsP0;
  
  /// flux point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_flxPntsLocalCoords;

  Common::SafePtr< std::vector< RealVector > > m_flxPntsLocalCoordsP0;
  
  /// flx pnt - face connectivity per orient
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_faceFlxPntConnPerOrient;

  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_faceFlxPntConnPerOrientP0;

  /// flx pnt - face connectivity
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;

  /// flx pnt - face connectivity for P0
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConnP0;   

  /// face connectivity per orient
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceConnPerOrient;
  
  /// number of equations in the physical model
  CFuint m_nbrEqs;
  
  /// number of dimensions in the physical model
  CFuint m_dim;
  
  /// variable for current face orientation
  CFuint m_orient;
  
  /// number of solution pnts in the cell
  CFuint m_nbrSolPnts;
  
  /// number of face flx pnts
  CFuint m_nbrFaceFlxPnts;

  CFuint m_nbrFaceFlxPntsP0;
  
  /// variable for current face
  Framework::GeometricEntity* m_face;
  
  /// variable for current neighbouring cells
  std::vector< Framework::GeometricEntity* > m_cells;
  
  /// Riemann flux
  Common::SafePtr< RiemannFlux > m_riemannFluxComputer;
  
  /// Correction function computer
  Common::SafePtr< BaseCorrectionFunction > m_corrFctComputer;
  
  /// Divergence of the correction function for current cell
  std::vector< std::vector< CFreal > > m_corrFctDiv;

  /// Divergence of the P0 correction function
  std::vector< std::vector< CFreal > > m_corrFctDivP0;
  
  std::vector< RealVector > m_P0State;
  std::vector< std::vector< RealVector > > m_statesP0;

  /// variable for the states in the left and right cell
  std::vector< std::vector< Framework::State* >* > m_states;

  /// vector of filtered interface fluxes
  std::vector< std::vector< RealVector > > m_Filtered_flxPntRiemannFlux;
  
  /// Interface fluxes at the flux points of a face
  std::vector< RealVector> m_flxPntRiemannFlux;

  std::vector< RealVector> m_flxPntRiemannFluxP0;
  
  /// Continuous flux at the solution points
  std::vector< std::vector< RealVector> > m_contFlx;

  std::vector< std::vector< RealVector> > m_contFlxP0;

  /// Filtered discontinuous flux at the solution points
  std::vector< std::vector< RealVector> > m_Filtered_DiscontFlux;
  
  /// Vandermonde Matrix
  RealMatrix  m_vdm;

  /// Inverse of the Vandermonde
  RealMatrix m_vdmInv;

  /// Face type index
  CFuint m_faceType;

  /// local cell flux point - face connectivity
  std::vector< CFuint > m_flxPntFaceConn;
  
  /// Divergence of the continuous flux at the solution points
  std::vector< RealVector> m_divContFlx;
  
  /// updates for the wave speed
  std::vector< CFreal > m_waveSpeedUpd;

  std::vector< CFreal > m_waveSpeedUpdP0;
  
  /// face Jacobian vector sizes (abs)
  std::vector< CFreal > m_faceJacobVecAbsSizeFlxPnts;

  std::vector< CFreal > m_faceJacobVecAbsSizeFlxPntsP0;
  
  /// coefficients for integration over a face
  Common::SafePtr< RealVector > m_faceIntegrationCoefs;
  
  Common::SafePtr< RealVector > m_faceIntegrationCoefsP0;

  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< std::vector< CFint > > > m_faceMappedCoordDir;
  
  /// local cell face - mapped coordinate direction
  Common::SafePtr< std::vector< CFint > > m_faceLocalDir;
  
  /// unit normal vector in flux points
  std::vector< RealVector > m_unitNormalFlxPnts;

  std::vector< RealVector > m_unitNormalFlxPntsP0;
  
  /// face Jacobian vector sizes
  std::vector< std::vector< CFreal > > m_faceJacobVecSizeFlxPnts;

  std::vector< std::vector< CFreal > > m_faceJacobVecSizeFlxPntsP0;
  
  /// flux point coordinates
  std::vector< RealVector > m_flxPntCoords;
  
  /// flux projection vectors in solution points for disc flux
  std::vector< std::vector< RealVector > > m_cellFluxProjVects;

  std::vector< std::vector< RealVector > > m_cellFluxProjVectsP0;
  
  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPntsP0;

  /// coefs to compute the derivative of the states in the sol pnts
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPnts;
  
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPntsP0;

  /// dimensions on which to evaluate the flux in the flux points
  Common::SafePtr< std::vector< CFuint > >  m_flxPntFlxDim;

  Common::SafePtr< std::vector< CFuint > >  m_flxPntFlxDimP0;
  
  /// the discontinuous flux extrapolated to the flux points
  std::vector< RealVector > m_extrapolatedFluxes;

  std::vector< RealVector > m_extrapolatedFluxesP0;
  
  std::vector< RealVector > m_Filtered_extrapolated;
  
  /// face local coordinates of the flux points on one face
  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoords;

  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoordsP0;

  /// local coordinates of the flux points on one face per face type
  Common::SafePtr<std::vector< std::vector< RealVector > > > m_faceFlxPntsLocalCoordsPerType;

  Common::SafePtr<std::vector< std::vector< RealVector > > > m_faceFlxPntsLocalCoordsPerTypeP0;

  std::vector< Framework::State* >* m_cellStatesP0;

  /// coefficients for integration over a face per face type
  Common::SafePtr<std::vector<  RealVector > > m_faceIntegrationCoefsPerType;

  Common::SafePtr<std::vector<  RealVector > > m_faceIntegrationCoefsPerTypeP0;

  /// dependencies of flx pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_flxSolDep;

  Common::SafePtr< std::vector< std::vector< CFuint > > > m_flxSolDepP0;

  /// dependencies of solution pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solSolDep;

  /// dependencies of flx pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solFlxDep;

  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solFlxDepP0;

  /// cell average solution coefficients for weighting P0 contributions
  Common::SafePtr< RealVector > m_cellAvgSolCoefs;

  /// nbr of sol pnts on which a flx pnt is dependent
  CFuint m_nbrSolDep;

  /// nbr of flx pnts a sol pnt influences
  CFuint m_nbrFlxDep;

  CFuint m_nbrFlxDepP0;

  /// nbr of sol pnts a sol pnt influences
  CFuint m_nbrSolSolDep;

  /// list of dimensions in which the flux will be evaluated in each sol pnt
  std::vector< std::vector< CFuint > > m_dimList;

  /// vector to store the face jacobians in
  std::vector< RealVector > m_faceJacobVecs;

  std::vector< RealVector > m_faceJacobVecsP0;
  
  /// left correction term projected on a normal
  RealVector m_projectedCorrL;
  
  /// right correction term projected on a normal
  RealVector m_projectedCorrR;
  
  /// array for jacobian determinants in sol pnts
  std::valarray<CFreal> m_jacobDet;

  /// Physical data temporary vector
  RealVector m_pData;
  
  /// Divergence of the continuous flux at the solution points of the left neighbour
  std::vector< RealVector> m_divContFlxL;
  
  /// Divergence of the continuous flux at the solution points of the right neighbour
  std::vector< RealVector> m_divContFlxR;
  
  /// FR order
  CFuint m_order;

  /// number of additionnal face normal directions per element for Triag (,terta and prism)
  CFuint m_ndimplus;

  /// Factor correcting Face normals direction (-1 factor needed for Tetra, due to the numbering convention the face normals are pointing inwards)
  CFreal m_mappedFaceNormalDir;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionBlending_hh

