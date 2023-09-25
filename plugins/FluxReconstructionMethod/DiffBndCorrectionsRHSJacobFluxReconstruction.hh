#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for an implicit scheme
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   *
   */
class DiffBndCorrectionsRHSJacobFluxReconstruction : public DiffBndCorrectionsRHSFluxReconstruction {

public:
  typedef Framework::BaseMethodCommandProvider<
      FluxReconstructionSolverData,DiffBndCorrectionsRHSJacobFluxReconstruction > PROVIDER;

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSJacobFluxReconstruction(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSJacobFluxReconstruction();

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

protected: // functions

  /**
   * Execute on the current TRS
   */
  virtual void executeOnTrs();
  
  /**
   * set the local indexes of the other faces (not the current boundary face)
   * @pre m_faces is set
   */
  void setOtherFacesLocalIdxs();

  /**
   * set the face neighbour states
   * @pre setOtherFacesLocalIdxs()
   */
  void setFaceNeighbourStates();
  
  /**
   * compute the contribution of the diffusive boundary face term to the Jacobian
   */
  void computeJacobDiffBndContribution();
  
  /**
   * recompute the cell gradients from the current cell and the neighbouring cells,
   * after perturbation
   */
  void computePerturbedGradients();
  
  /**
   * recompute analytically the cell gradients from the current cell and the neighbouring cells,
   * after perturbation
   */
  void computePerturbedGradientsAnalytical();
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);
  
  /**
   * compute the term for the gradient computation for the cell
   */
  virtual void computeCellGradTerm(RealMatrix& gradTerm);
  
  /**
   * compute the terms for the gradient computation for a face
   */
  virtual void computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR);
  
  /**
   * compute the contribution of the convective boundary flux correction to the Jacobian
   */
  void extrapolatePerturbedState();
  
  /**
   * store backups of values before perturbing the states
   */
  void storeBackups();
  
  /**
   * restore values after perturbing a state
   */
  void restoreFromBackups();

protected: // data
  
  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to the numerical Jacobian computer
  Common::SafePtr<Framework::NumericalJacobian> m_numJacob;

  /// accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> m_acc;

  /// variable for faces
  const std::vector< Framework::GeometricEntity* >* m_faces;

  /// vector containing pointers to the left and right states with respect to a face
  std::vector< std::vector< std::vector< Framework::State* >* > > m_faceNghbrStates;

  /// perturbed updates to the residuals
  RealVector m_pertResUpdates;

  /// derivative of update to one CV-residual
  RealVector m_derivResUpdates;

  /// updates to the gradients
  std::vector< std::vector< std::vector< RealVector > > > m_gradUpdates;

  /// perturbed boundary cell gradients
  std::vector< std::vector< RealVector >* > m_pertGrads;

  /// Jacobian determinants
  std::valarray< CFreal > m_solJacobDet;

  /// cell local indexes of the other faces (not the boundary face itself)
  std::vector< CFuint > m_otherFaceLocalIdxs;

  /// pointer to booleans telling whether a face is on the boundary
  Common::SafePtr< std::vector< bool > > m_isFaceOnBoundary;

  /// pointer to neighbouring cell side vector
  Common::SafePtr< std::vector< CFuint > > m_nghbrCellSide;

  /// pointer to current cell side vector
  Common::SafePtr< std::vector< CFuint > > m_currCellSide;

  /// pointer to orientation vector
  Common::SafePtr< std::vector< CFuint > > m_faceOrients;

  /// pointer to BC index vector
  Common::SafePtr< std::vector< CFuint > > m_faceBCIdx;

  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;
  
  /// perturbed corrections due to the boundary faces for the Jacobian
  std::vector< RealVector> m_pertCorrections;
  
  /// unperturbed updates to the residuals
  RealVector m_resUpdates;
  
  /// backup of the gradients in the neighbouring cell
  std::vector< std::vector< RealVector > > m_cellGradsBackUp;
  
  /// coefs to compute the derivative of the states in the sol pnts
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPnts;
  
  /// flx pnt - face connectivity per orient
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_faceFlxPntConnPerOrient;
  
  /// extrapolated perturbed states in the flux points of the cell
  std::vector< std::vector< Framework::State* > > m_pertCellStatesFlxPnt;
  
  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< std::vector< CFint > > > m_faceMappedCoordDirPO;
  
  /// index of the perturbed solution point
  CFuint m_pertSol;
  
  /// index of the perturbed variable
  CFuint m_pertVar;
  
  /// dependencies of sol pnts on flx pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solFlxDep;

  /// nbr of flx pnts on which a sol pnt is dependent
  CFuint m_nbrFlxDep;
  
  /// backup of extrapolated states in the flux points of the cell
  std::vector< RealVector > m_cellStatesFlxPntBackup;
  
  /// influenced flx pnt idx (by perturbation)
  CFuint m_influencedFlxPnt;

  /// influenced flx pnts idx (by perturbation)
  std::vector< CFuint> m_influencedFlxPnts;

  /// Number of influenced flx pnts (by perturbation)
  CFuint m_NbInfluencedFlxPnts;
  
  /// Element shape
  CFGeoShape::Type elemShape;
  
  /// backup of interface fluxes at the flux points of a face
  std::vector< RealVector> m_flxPntRiemannFluxBackup;
  
  /// the corrected gradients in the flux points backup
  std::vector< std::vector< RealVector* > > m_cellGradFlxPntBackup;
  
  /// list of dimensions in which the flux will be evaluated in each sol pnt
  std::vector< std::vector< CFuint > > m_dimList;
  
  /// flux projection vectors in solution points for disc flux
  std::vector< std::vector< RealVector > > m_cellFluxProjVects;
  
  /// term for the gradient computation
  RealMatrix m_gradTerm;
  
  /// term for the gradient computation before perturbation
  RealMatrix m_gradTermBefore;
  
  /// dependencies of solution pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solSolDep;

  /// nbr of sol pnts a sol pnt influences
  CFuint m_nbrSolSolDep;
  
  /// correction projected on a normal
  RealVector m_projectedCorr;
  
  /// face local coordinates of the flux points on one face
  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoords;

  /// local coordinates of the flux points on one face per face type
  Common::SafePtr<std::vector< std::vector< RealVector > > > m_faceFlxPntsLocalCoordsPerType;
  
  /// vector to store the face jacobians in
  std::vector< RealVector > m_faceJacobVecs;
  
  /// term for the gradient computation on a face
  RealMatrix m_gradTermFace;
  
  /// term for the ghost gradient computation on a face
  RealMatrix m_ghostGradTerm;
  
  /// temp term for the gradient computation on a face
  RealMatrix m_gradTermTemp;
  
  /// array for jacobian determinants in sol pnts
  std::valarray<CFreal> m_jacobDet;
  
  /// number of flx pnts in one element
  CFuint m_nbrTotalFlxPnts;
  
  /// extrapolated states in the flux points of the cell
  std::vector< Framework::State* > m_cellStatesFlxPnt2;
  
  /// perturbations
  RealVector m_eps;

  /// number of additionnal face normal directions for Triag (,terta and prism)
  CFuint m_ndimplus;

}; // end of class DiffBndCorrectionsRHSJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstruction_hh
