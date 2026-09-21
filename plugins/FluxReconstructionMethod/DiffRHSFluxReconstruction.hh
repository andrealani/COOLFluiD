// Copyright (C) 2016 KU Leuven, Belgium
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

/// This is a standard command to assemble the (diffusive part of the) system using a FluxReconstruction solver.
///
/// The gradients follow the compact BR2 scheme, written like the fluxes. For
/// one cell, g = g(U) are the gradient variables of the diffusive variable
/// set, g^D the polynomial through their values at the solution points, g^D_f
/// its value extrapolated to the flux points of face f, g^I_f the interface
/// value of face f (the average of the two sides at an interior face, the
/// boundary value of the boundary condition at a boundary face) and h_f the
/// correction function of face f, as F^D, F^D_f, F^I_f and h_f for the flux.
/// The gradient the cell flux uses is corrected with every face of the cell,
///
///   q = grad g^D + sum_f (g^I_f - g^D_f) grad h_f
///
/// with grad h_f the unit normal of the face, scaled by the face Jacobian,
/// times div h_f. The interface flux uses, per side, the compact gradient of
/// that face alone, scaled by eta and taken at the flux points of the face,
///
///   q_f = grad g^D + eta (g^I_f - g^D_f) grad h_f
///
/// so the face flux depends on the two cells of the face only. The interface
/// diffusive flux of face f is the diffusive flux of the variable set at the
/// average of the two extrapolated states and at the average of the two
/// sides' q_f, with no penalty term; the volume flux of the cell uses q. eta is
/// the BR2 lifting multiplier (option BR2Eta, 5 by default): it sets how
/// strongly the jump of the face enters the gradient the interface flux sees,
/// and the usual BR2 condition wants it larger than the number of faces of the
/// element. Both gradients are built in mapped coordinates and divided by the
/// Jacobian determinant.
///
/// @author Alexander Papen
/// @author Ray Vandenhoeck
/// @author Rayan Dhib
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

  /**
   * Compute the common diffusive flux at the flux points of the current face:
   * the diffusive flux of the average extrapolated state and the average of
   * the two compact BR2 face gradients.
   */
  virtual void computeInterfaceFlxCorrection();

  /**
   * Fill the per-side cell metrics of the current face (solution point Jacobian
   * determinants and mapped coordinate plane normals) for both neighbouring
   * cells, from m_cells.
   */
  virtual void prepareFaceCellMetrics();

  /**
   * Set the gradient variables of the given states at the solution points of a
   * cell, one column per solution point.
   */
  virtual void computeCellGradVars(const std::vector< Framework::State* >& states, RealMatrix& gradVars);

  /**
   * Compute the compact BR2 gradient q_f of the current face for both cells and
   * extrapolate it to the flux points of the face. Per side it is the volume
   * term grad g^D of that cell plus eta times the correction
   * (g^I_f - g^D_f) grad h_f of this face, with g^I_f the average of the two
   * sides, divided by the Jacobian determinant and extrapolated to the flux
   * points. The correction of this face is the only one in it, so the result
   * depends on the two cells of the face alone. m_cellGrads is not changed.
   * @pre prepareFaceCellMetrics()
   * @param gradVarsSolPntsL gradient variables at the solution points of the left cell, nbrEqs x nbrSolPnts
   * @param gradVarsSolPntsR gradient variables at the solution points of the right cell, nbrEqs x nbrSolPnts
   * @param gradVarsFlxPntL gradient variables of the left cell extrapolated to the face flux points, nbrEqs x nbrFaceFlxPnts
   * @param gradVarsFlxPntR gradient variables of the right cell extrapolated to the face flux points, nbrEqs x nbrFaceFlxPnts
   * @param faceGrads output [side][iFlx][iEq], or CFNULL to write into m_cellGradFlxPnt
   */
  void computeCompactBR2FaceGradients(const RealMatrix& gradVarsSolPntsL, const RealMatrix& gradVarsSolPntsR,
                                      const RealMatrix& gradVarsFlxPntL, const RealMatrix& gradVarsFlxPntR,
                                      std::vector< std::vector< std::vector< RealVector* > > >* faceGrads = CFNULL);

  /**
   * Add the correction (g^I_f - g^D_f) grad h_f of the current interior face to
   * the gradients of both cells, with g^I_f the average of the two sides and g
   * the physical gradient variables or, for the artificial viscosity, the
   * conservative variables.
   * @param artificialViscosity true to add them to the artificial viscosity gradients
   */
  void addGradientFaceCorrections(const bool artificialViscosity = false);

  /**
   * Add the volume term grad g^D of the current cell m_cell to its gradients and
   * divide the result by the Jacobian determinant at the solution points.
   * @param artificialViscosity true for the artificial viscosity gradients
   */
  void addGradientVolumeTerm(const bool artificialViscosity = false);

  /**
   * Set the variables of the artificial viscosity gradients of the given
   * states, one column per state: the conservative variables.
   * @param nbrStates number of states to transform
   */
  void setAVGradientVars(const std::vector< Framework::State* >& states, const CFuint nbrStates, RealMatrix& values);

  /**
   * Compute the compact BR2 gradient of the artificial viscosity variables of
   * the current face for both cells.
   * @param faceGrads output [side][iFlx][iEq], or CFNULL to write into m_cellGradFlxPnt
   */
  void computeCompactBR2FaceGradientsAV(std::vector< std::vector< std::vector< RealVector* > > >* faceGrads = CFNULL);

  /**
   * Compute the compact BR2 gradient of the artificial viscosity variables at
   * the flux points of one boundary face of a cell: the volume term grad g^D
   * and eta times the correction of this face, with the interface value
   * g^I_f = 0.5*(g^D_f + g(U_ghost)), the average of the extrapolated values
   * and the values of the ghost states.
   * @param iFace local index of the face in the cell
   * @param unitNormals unit normals at the flux points of the face
   * @param flxPntCoords coordinates of the flux points of the face
   * @param faceGrads output [iFlx][iEq]
   */
  void computeCompactBR2BndFaceGradientAV(Framework::GeometricEntity& cell,
                                          Framework::GeometricEntity& face,
                                          const CFuint iFace,
                                          Common::SafePtr< BCStateComputer > bc,
                                          const std::vector< RealVector >& unitNormals,
                                          const std::vector< RealVector >& flxPntCoords,
                                          std::vector< std::vector< RealVector* > >& faceGrads);

  /**
   * Compute the artificial viscosity residual of the boundary faces of a cell:
   * at every boundary flux point the common flux eps q_f.n, with eps the
   * artificial viscosity extrapolated to the flux point, q_f the compact BR2
   * gradient of the face and n the unit normal, lifted into the cell with the
   * correction functions.
   * @param isFaceOnBoundary tells for every face of the cell whether it is a boundary face
   * @param faceBCIdx index of the boundary condition of every face of the cell
   * @param epsilons artificial viscosity at the solution points of the cell
   * @param residual output, size nbrSolPnts*nbrEqs
   */
  void computeBndFacesAVResidual(Framework::GeometricEntity& cell,
                                 const std::vector< bool >& isFaceOnBoundary,
                                 const std::vector< CFuint >& faceBCIdx,
                                 const std::vector< CFreal >& epsilons,
                                 RealVector& residual);

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
  virtual void computeFlxPntStatesAndGrads();
  
  virtual void computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux);
  
  /// prepare the computation of the diffusive flux
  virtual void prepareFluxComputation()
  {
  }

  /**
   * Prepares the evaluation of the diffusive flux at a solution point, for the
   * physics that need point data beyond the state and the gradient (the wall
   * distance of the turbulence models).
   * @param stateID local ID of the state of the solution point
   * Default: prepareFluxComputation().
   */
  virtual void prepareSolPntFluxComputation(const CFuint stateID)
  {
    prepareFluxComputation();
  }

  /**
   * Prepares the evaluation of the diffusive flux at a flux point of the current
   * face, see prepareSolPntFluxComputation().
   * @param iFlx index of the flux point on the face
   * Default: prepareFluxComputation().
   */
  virtual void prepareFlxPntFluxComputation(const CFuint iFlx)
  {
    prepareFluxComputation();
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

  /// Max number of face flx pnts (relevant only for prisms)
  CFuint m_nbFaceFlxPntsMax;

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
  
  /// coefficients for integration over a face per face type
  Common::SafePtr<std::vector<  RealVector > > m_faceIntegrationCoefsPerType;

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
  
  /// the discontinuous flux extrapolated to the flux points
  std::vector< RealVector > m_extrapolatedFluxes;
  
  /// dimensions on which to evaluate the flux in the flux points
  Common::SafePtr< std::vector< CFuint > >  m_flxPntFlxDim;
  
  /// average solution in a flux point
  RealVector m_avgSol;
  
  /// average gradients in a flux point
  std::vector< RealVector* > m_avgGrad;

  /// face local coordinates of the flux points on one face
  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoords;

  /// local coordinates of the flux points on one face per face type
  Common::SafePtr<std::vector< std::vector< RealVector > > > m_faceFlxPntsLocalCoordsPerType;

  /// dependencies of flx pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_flxSolDep;

  /// dependencies of solution pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solSolDep;

  /// dependencies of flx pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solFlxDep;

  /// nbr of sol pnts on which a flx pnt is dependent
  CFuint m_nbrSolDep;

  /// nbr of flx pnts a sol pnt influences
  CFuint m_nbrFlxDep;

  /// nbr of sol pnts a sol pnt influences
  CFuint m_nbrSolSolDep;

  /// list of dimensions in which the flux will be evaluated in each sol pnt
  std::vector< std::vector< CFuint > > m_dimList;

  /// vector to temporarily store the gradients
  std::vector< RealVector* > m_tempGrad;

  /// number of flx pnts in a cell
  CFreal m_nbrTotalFlxPnts;
  
  /// Factor correcting Face normals direction (-1 factor needed for Tetra, due to the numbering convention the face normals are pointing inwards)
  CFreal m_mappedFaceNormalDir;

  /// vector with the face Jacobian vectors
  std::vector< RealVector > m_faceJacobVecs;
  
  /// vector of arrays with jacobian determinants for the sol pnts
  std::vector< std::valarray<CFreal> > m_jacobDets;
  
  /// FR order
  CFuint m_order;

  /// Element shape/type
  CFGeoShape::Type elemShape;
  
  /// number of additionnal face normal directions for Triag (,terta and prism)
  CFuint m_ndimplus;

  /// Coeff (-1 or 1) indicating the direction of the computed mapped normals for Triag =-1 (normals pointing outwards)
  CFreal m_mappedNormalDir;
  
  bool m_addRiemannToGradJacob;
  
  bool m_addRiemannToGradCrossCellJacob;
  
  bool m_addFluxToGradCrossCellJacob;

  /// multiplier of the face lifting in the compact BR2 face gradient
  CFreal m_br2Eta;

  /// solution point Jacobian determinants of both cells of the current face
  std::vector< std::valarray<CFreal> > m_solJacobDet;

  /// mapped coordinate plane normals at the solution points of both cells of the current face
  std::vector< std::vector< std::vector< RealVector > > > m_neighbCellFluxProjVects;

  /// tells whether prepareFaceCellMetrics() filled the per-side cell metrics of the current face
  bool m_faceCellMetricsPrepared;

  /// compact BR2 face gradient at the solution points of both cells of the current face, [side][iSol][iEq]
  std::vector< std::vector< std::vector< RealVector > > > m_compactGradsSolPnts;

  /// correction projected on a normal, size dim
  RealVector m_projectedCorr;

  /// gradient variables at the solution points of both cells of the current face, [side](iEq,iSol)
  std::vector< RealMatrix > m_gradVarsSolPnts;

  /// gradient variables extrapolated to the flux points of the current face, [side](iEq,iFlx)
  std::vector< RealMatrix > m_gradVarsFlxPnt;

  /// gradient variables extrapolated to one flux point
  RealVector m_flxPntGradVars;

  /// solution point state data passed to the diffusive variable set
  std::vector< RealVector* > m_gradVarStatePtrs;

  /// states extrapolated to the flux points of a boundary face, for the artificial viscosity boundary gradient
  std::vector< Framework::State* > m_bndIntStates;

  /// ghost states at the flux points of that boundary face
  std::vector< Framework::State* > m_bndGhostStates;

  /// the first nbrFaceFlxPnts states of m_bndIntStates, for every number of face flux points
  std::vector< std::vector< Framework::State* > > m_bndIntStatesFlxPnt;

  /// the first nbrFaceFlxPnts states of m_bndGhostStates, for every number of face flux points
  std::vector< std::vector< Framework::State* > > m_bndGhostStatesFlxPnt;

  /// unit normals at the flux points of a boundary face, for every number of face flux points
  std::vector< std::vector< RealVector > > m_bndUnitNormalFlxPnts;

  /// coordinates of the flux points of a boundary face, for every number of face flux points
  std::vector< std::vector< RealVector > > m_bndFlxPntCoords;

  /// gradient variables of the ghost states, (iEq,iFlx)
  RealMatrix m_bndGhostGradVars;

  /// mapped coordinate plane normals at the solution points of the cell of a boundary face, [iDim][iSol]
  std::vector< std::vector< RealVector > > m_bndCellFluxProjVects;

  /// Jacobian determinants at the solution points of the cell of a boundary face
  std::valarray< CFreal > m_bndSolJacobDet;

  /// artificial viscosity compact gradient at the flux points of a boundary face, [iFlx][iEq]
  std::vector< std::vector< RealVector > > m_bndFaceGradsAV;

  /// pointers to m_bndFaceGradsAV, [iFlx][iEq]
  std::vector< std::vector< RealVector* > > m_bndFaceGradPtrsAV;

  /// artificial viscosity common flux at a boundary flux point
  RealVector m_bndFlxPntFluxAV;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstruction_hh

