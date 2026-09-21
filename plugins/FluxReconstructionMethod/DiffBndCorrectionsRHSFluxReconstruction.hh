#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/BndFaceDiffData.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   * @author Rayan Dhib
   *
   */
class DiffBndCorrectionsRHSFluxReconstruction : public FluxReconstructionSolverCom {

public:
  typedef Framework::BaseMethodCommandProvider<
      FluxReconstructionSolverData,DiffBndCorrectionsRHSFluxReconstruction > PROVIDER;

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSFluxReconstruction(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSFluxReconstruction();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * set the BC state computer
   */
  void setBcStateComputer(Common::SafePtr< BCStateComputer > bcStateComputer)
  {
    m_bcStateComputer = bcStateComputer;
  }

  /**
   * Computes the data of the diffusive boundary flux at the flux points of a
   * boundary face built by the caller, with the boundary condition and the
   * diffusive flux of this command. The RHS is not updated.
   * @param face    the boundary face
   * @param orient  orientation of the face
   * @param data    the data of the diffusive boundary flux
   */
  void computeBndFaceDiffData(Framework::GeometricEntity* face,
                              const CFuint orient,
                              BndFaceDiffData& data);

protected: // functions

  /**
   * Execute on the current TRS
   */
  virtual void executeOnTrs();

  /// compute the states and ghost states in the flx pnts
  virtual void computeFlxPntStates();

  /**
   * Computes the compact BR2 gradient of the current boundary face at its flux
   * points, in m_cellGradFlxPnt: the volume term of the gradient variables of
   * the interior cell plus eta, the BR2 multiplier, times the correction of
   * this face alone, from the gradient variables extrapolated to the flux
   * points to the boundary values of the boundary condition. m_cellGrads is
   * not changed.
   * @pre computeFlxPntStates() and setIntCellMetrics()
   */
  void computeCompactBR2BndFaceGradient();

  /// set the cell metrics of the interior cell needed by the compact face gradient
  void setIntCellMetrics();
  
  /**
   * Compute the diffusive flux at the flux points of the current boundary face,
   * at the boundary state U_b and with the boundary gradients q_b that the
   * boundary condition returns from the compact face gradient.
   */
  virtual void computeInterfaceFlxCorrection();
  
  /// compute the total correction
  void computeCorrection(std::vector< RealVector >& corrections);

  /// add the residual updates to the RHS
  void updateRHS();

  /// add the updates to the wave speed
  void updateWaveSpeed();
  
  /// set the bnd face data necessary to compute FI-FD
  virtual void setBndFaceData(CFuint faceID);
  
  /// prepare the computation of the diffusive flux
  virtual void prepareFluxComputation()
  {
  }

  /**
   * Prepares the evaluation of the diffusive flux at a flux point of the current
   * boundary face, for the physics that need point data beyond the state and the
   * gradient (the wall distance of the turbulence models).
   * @param iFlx index of the flux point on the face
   * Default: prepareFluxComputation().
   */
  virtual void prepareFlxPntFluxComputation(const CFuint iFlx)
  {
    prepareFluxComputation();
  }
  
  virtual void computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux);
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeWaveSpeedUpdates(CFreal& waveSpeedUpd);

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// socket for gradientsAV
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradientsAV;
  
  /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// the BCStateComputer for this BC
  Common::SafePtr< BCStateComputer > m_bcStateComputer;

  /// variable for current face
  Framework::GeometricEntity* m_face;

  /// variable for current cell
  Framework::GeometricEntity* m_intCell;

  /// variable for current face orientation
  CFuint m_orient;
  
  /// number of dimensions in the physical model
  CFuint m_dim;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;

  /// update for the wave speed in the neighbouring cell
  CFreal m_waveSpeedUpd;
  
  /// the gradients in the neighbouring cell
  std::vector< std::vector< RealVector >* > m_cellGrads;
  
  /// the corrected gradients in the flux points
  std::vector< std::vector< RealVector* > > m_cellGradFlxPnt;
  
  /// the ghost gradients in the flux points
  std::vector< std::vector< RealVector* > > m_flxPntGhostGrads;
  
  /// cell volume
  CFreal m_cellVolume;
  
  /// ratio between convective and diffusive cfl limit
  CFreal m_cflConvDiffRatio;

  /// number of equations in the physical model
  CFuint m_nbrEqs;
  
  /// number of solution pnts in the cell
  CFuint m_nbrSolPnts;
  
  /// number of flux pnts on a face
  CFuint m_nbrFaceFlxPnts;
  
  /// Factor correcting Face normals direction (-1 factor needed for Tetra, due to the numbering convention the face normals are pointing inwards)
  CFreal m_mappedFaceNormalDir;
  
  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;
  
  // All flux points of a cell
  Common::SafePtr<std::vector< RealVector > > m_allCellFlxPnts;
  
  /// flux points mapped coordinates
  std::vector< RealVector > m_flxPntsLocalCoords;
  
  /// flux point coordinates
  std::vector< RealVector > m_flxPntCoords;
  
  /// flx pnt - face connectivity
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;
  
  /// face connectivity per orient
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceConnPerOrient;
  
  /// coefficients for integration over a face
  Common::SafePtr< RealVector > m_faceIntegrationCoefs;

  /// coefficients for integration over a face per face type
  Common::SafePtr<std::vector<  RealVector > > m_faceIntegrationCoefsPerType;
  
  /// face Jacobian vector sizes (abs)
  std::vector< CFreal > m_faceJacobVecAbsSizeFlxPnts;
  
  /// extrapolated states in the flux points of the cell
  std::vector< Framework::State* > m_cellStatesFlxPnt;
  
  /// vector containing pointers to the fluxes in the flux points
  std::vector< RealVector > m_cellFlx;
  
  /// ghost flux point solutions
  std::vector< Framework::State* > m_flxPntGhostSol;
  
  /// Riemann flux
  Common::SafePtr< RiemannFlux > m_riemannFluxComputer;
  
  /// Interface fluxes at the flux points of a face
  std::vector< RealVector> m_flxPntRiemannFlux;
  
  /// Correction function computer
  Common::SafePtr< BaseCorrectionFunction > m_corrFctComputer;
  
  /// Divergence of the correction function
  std::vector< std::vector< CFreal > > m_corrFctDiv;
  
  /// corrections to be added to RHS because of the boundary faces
  std::vector< RealVector> m_corrections;
  
  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< CFint > > m_faceMappedCoordDir;
  
  /// unit normal vector in flux points
  std::vector< RealVector > m_unitNormalFlxPnts;
  
  /// face Jacobian vector sizes
  std::vector< CFreal > m_faceJacobVecSizeFlxPnts;
  
  /// diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  /// average solution in a flux point
  RealVector m_avgSol;

  /// storage of m_bndStateFlux
  std::vector< RealVector > m_bndStateFluxStore;

  /// boundary states the diffusive flux is evaluated at
  std::vector< RealVector* > m_bndStateFlux;
  
  /// average gradients in a flux point
  std::vector< RealVector* > m_avgGrad;

  /// local coordinates of the flux points on one face
  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoords;

  /// local coordinates of the flux points on one face per face type
  Common::SafePtr<std::vector< std::vector< RealVector > > > m_faceFlxPntsLocalCoordsPerType;

  /// dependencies of flx pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_flxSolDep;

  /// nbr of sol pnts on which a flx pnt is dependent
  CFuint m_nbrSolDep;
  
  /// vector with face Jacobian vectors
  std::vector< RealVector > m_faceJacobVecs;

  /// Physical data temporary vector
  RealVector m_pData;
  
  /// FR order
  CFuint m_order;

  /// multiplier of the face lifting in the compact BR2 face gradient
  CFreal m_br2Eta;

  /// number of additional face normal directions for Triag (, tetra and prism)
  CFuint m_ndimplus;

  /// list of dimensions in which the flux will be evaluated in each sol pnt
  std::vector< std::vector< CFuint > > m_dimList;

  /// mapped coordinate plane normals at the solution points of the interior cell
  std::vector< std::vector< RealVector > > m_cellFluxProjVects;

  /// solution point Jacobian determinants of the interior cell
  std::valarray< CFreal > m_solJacobDet;

  /// coefs to compute the derivative of the states in the sol pnts
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPnts;

  /// dependencies of solution pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solSolDep;

  /// nbr of sol pnts a sol pnt influences
  CFuint m_nbrSolSolDep;

  /// correction projected on a normal
  RealVector m_projectedCorr;

  /// gradient variables at the solution points of the interior cell
  RealMatrix m_gradVarsSolPnts;

  /// solution point state data passed to the diffusive variable set
  std::vector< RealVector* > m_gradVarStatePtrs;

  /// compact BR2 gradient of the current boundary face at the solution points [iSol][iEq]
  std::vector< std::vector< RealVector > > m_compactGradsSolPnts;

  /// storage of m_flxPntGradVars
  std::vector< RealVector > m_flxPntGradVarsStore;

  /// gradient variables extrapolated to the flux points of the current boundary face
  std::vector< RealVector* > m_flxPntGradVars;

  /// storage of m_bndGradVars
  std::vector< RealVector > m_bndGradVarsStore;

  /// boundary values of the gradient variables at the flux points
  std::vector< RealVector* > m_bndGradVars;

  /// boundary gradients of the diffusive boundary flux at the flux points [iFlx][iEq]
  std::vector< std::vector< RealVector* > > m_bndGradFlxPnt;

private: // functions

  /**
   * Resizes the members of data to the current number of flux points.
   * @param data  the data of the diffusive boundary flux
   */
  void resizeBndFaceDiffData(BndFaceDiffData& data);

}; // end of class DiffBndCorrectionsRHSFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstruction_hh
