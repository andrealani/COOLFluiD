#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BaseBndFaceTermComputer_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BaseBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "MathTools/MatrixInverter.hh"

#include "FluxReconstructionMethod/BCStateComputer.hh"
//#include "FluxReconstructionMethod/FaceDiffusiveFlux.hh"
#include "FluxReconstructionMethod/ReconstructStatesFluxReconstruction.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms at the boundaries.
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 */
class BaseBndFaceTermComputer : public FluxReconstructionSolverStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      FluxReconstructionSolverData,BaseBndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  BaseBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~BaseBndFaceTermComputer();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BaseBndFaceTermComputer";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  
  /**
   * set the BCStateComputer
   */
  void setBcStateComputer(Common::SafePtr< BCStateComputer > bcStateComputer)
  {
    m_bcStateComputer = bcStateComputer;
  }

  /**
   * set the orientation of the faces
   */
  void setFaceOrientation(CFuint orient)
  {
    m_orient = orient;
  }

  /**
   * set the current face
   */
  void setCurrentFace(Framework::GeometricEntity* face)
  {
    m_face = face;
  }

  /**
   * returns m_unitNormalPntSet
   */
  Common::SafePtr< std::vector< RealVector > > getUnitNormalPntSet()
  {
    return &m_unitNormalPntSet;
  }

  /**
   * returns m_faceJacobPntSet
   */
  Common::SafePtr< std::vector< RealVector > > getFaceJacobPntSet()
  {
    return &m_faceJacobPntSet;
  }

  /**
   * Set up private data and data
   */
  virtual void setup();

  /**
   * Unset up private data and data
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * set face term data
   */
  virtual void setFaceTermData();

  /**
   * compute face data
   * @pre setCurrentFace and setFaceTermData()
   */
  void computeFaceData();

  /**
   * compute neighbour cell data
   * @pre setCurrentFace
   */
  virtual void computeNeighbourCellData();

  /**
   * reconstruct solution in the flux points
   * @pre setFaceTermData()
   */
  void reconstructFluxPntsStates(const std::vector< Framework::State* >& cellIntStates, bool onlyExtraVars = false);

  /**
   * reconstruct gradients in the flux points
   * @pre setFaceTermData()
   */
  void reconstructFluxPntsGradients(const std::vector< std::vector< RealVector >* >& cellIntGrads);

  /**
   * backup and reconstruct physical variable in the boundary face in the required points
   * and reconstruct the ghost states
   * @pre setVolumeTermData()
   */
  void backupAndReconstructPhysVar(const CFuint iVar, const std::vector< Framework::State* >& cellStates);

  /**
   * backup physical variable in the boundary face in the required points
   */
  void backupPhysVar(const CFuint iVar);

  /**
   * restore physical variable in the boundary face in the required points
   */
  void restorePhysVar(const CFuint iVar);

  /**
   * compute the convective face term for this face
   * @pre reconstructFluxPntsStates(), setFaceTermData() and set the geometrical data of the face
   */
  void computeConvFaceTerm(RealVector& resUpdates);

  /**
   * compute the convective face term and the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeConvFaceTermAndWaveSpeedUpdates(RealVector& resUpdates,CFreal& waveSpeedUpd);

  /**
   * face term contribution to the gradients
   * @pre reconstructFluxPntsStates(), setFaceTermData() and set the geometrical data of the face
   */
  void computeGradientFaceTerm(std::vector< std::vector< RealVector > >& gradUpdates);

  /**
   * face term contribution to the gradients of the extra variables
   * @pre reconstructFluxPntsStates(), setFaceTermData() and set the geometrical data of the face
   */
  void computeGradientExtraVarsFaceTerm(std::vector< std::vector< RealVector > >& gradUpdates);

  /**
   * compute the diffusive face term for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeDiffFaceTerm(RealVector& resUpdates);

  /**
   * compute the diffusive face term and the contribution to the update coefficient for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,CFreal& updateCoefContr);

  /**
   * set cell mapped coordinates of point set in which the states should be reconstructed
   * @pre setFaceOrientation()
   */
  virtual void setPointSet(const std::vector< RealVector >& faceMappedCoords);

  /**
   * compute face data for given set of points
   * @pre setCurrentFace and setOutputData()
   */
  virtual void computeFacePntSetData();

  /**
   * reconstruct solution in a given set of points
   * @pre setNodeSet()
   */
  std::vector< RealVector >& reconstructGivenPntsStates(const std::vector< Framework::State* >& cellIntStates);

  /**
   * reconstruct gradients in a given set of points
   * @pre setNodeSet()
   */
  std::vector< std::vector< RealVector > >&
      reconstructGivenPntsGrads(const std::vector< std::vector< RealVector >* >& cellIntGrads);

protected: // functions

  /**
   * compute the ghost states from the internal states, unit normals and coordinates
   * @pre internal solutions have to be reconstructed first
   */
  void computeGhostStates();

  /**
   * compute the perturbed ghost states
   * from the perturbed internal states, unit normals and coordinates
   * @pre internal solution has to be perturbed first
   * @pre computeGhostStates()
   */
  void computePerturbedGhostStates();

  /**
   * compute the ghost gradients from the internal gradients, unit normals and coordinates
   * @pre internal gradients have to be reconstructed first
   */
  virtual void computeGhostGradients();

  /**
   * compute the flux point coordinates
   * @pre setFaceTermData() and set the geometrical face data
   */
  void computeFlxPntCoords();

  /**
   * compute the flux integrals from the fluxes in the flux points
   * @pre fluxes in flux points are computed
   */
  void computeFaceTermFromFlxPntFluxes(RealVector& resUpdates);

  /**
   * compute face term contribution to the gradients from the solution in the flux points
   * @pre diffusiveVarSet->setGradientsVars()
   */
  void computeGradFaceTermFromFlxPntSol(std::vector< std::vector< RealVector > >& gradUpdates);

protected: // data

  /// the socket containing the extra variables
  Framework::DataSocketSink< RealVector > socket_extraVars;

  /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;

  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;

  /// Diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;

  /// Strategy that reconstructs the states in a given number of nodes
  Common::SafePtr< ReconstructStatesFluxReconstruction > m_statesReconstr;

  /// Strategy that sets the ghost states corresponding to the boundary condition
  Common::SafePtr< BCStateComputer> m_bcStateComputer;

  /// Riemann flux
  Common::SafePtr< RiemannFlux > m_riemannFluxComputer;

//   /// Face diffusive flux
//   Common::SafePtr< FaceDiffusiveFlux > m_faceDiffFluxComputer;

  /// ratio between convective and diffusive cfl limit
  CFreal m_cflConvDiffRatio;

  /// reconstruction coefficients for the flux points
  Common::SafePtr< RealMatrix > m_flxPntsRecCoefs;

  /// derivation coefficients for the solution points
  Common::SafePtr< RealMatrix > m_solPntsDerivCoefs;

  /// local cell face - flux point connectivity per face connection orientation
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;

  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< CFint > > m_faceMappedCoordDir;

  /// flux point index (in the matrix flxPntRecCoefs) for reconstruction
  Common::SafePtr< std::vector< CFuint > > m_flxPntMatrixIdxForReconstruction;

  /// solution point index (in the cell) for reconstruction
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solPntIdxsForReconstruction;

  /// flux point index (in the matrix m_solPntsDerivCoefs) for derivation
  Common::SafePtr< std::vector< CFuint > > m_flxPntMatrixIdxForDerivation;

  /// solution point index (in the cell) for derivation
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solPntIdxsForDerivation;

  /// coefficients for integration over a face
  Common::SafePtr< RealVector > m_faceIntegrationCoefs;

  /// face local coordinates of face flux points
  Common::SafePtr< std::vector< RealVector > > m_faceFlxPntsFaceLocalCoords;

  /// local cell face - flux point cell mapped coordinate per face connection orientation
  Common::SafePtr< std::vector< std::vector< RealVector > > > m_faceFlxPntCellMappedCoords;

  /// current face
  Framework::GeometricEntity* m_face;

  /// face orientation
  CFuint m_orient;

  /// cell volume
  CFreal m_cellVolume;

  /// unit normal vector in flux points
  std::vector< RealVector > m_unitNormalFlxPnts;

  /// face Jacobian vector sizes
  std::vector< CFreal > m_faceJacobVecSizeFlxPnts;

  /// flux projection vector sizes (abs)
  std::vector< CFreal > m_faceJacobVecAbsSizeFlxPnts;

  /// face inverse characteristic lengths
  std::vector< CFreal > m_faceInvCharLengths;

  /// internal cell extra variables
  std::vector< RealVector* > m_cellExtraVars;

  /// internal flux point solutions
  std::vector< Framework::State* > m_flxPntIntSol;

  /// ghost flux point solutions
  std::vector< Framework::State* > m_flxPntGhostSol;

  /// flux point extra variables
  std::vector< RealVector* > m_flxPntIntExtraVars;

  /// flux point extra variables
  std::vector< RealVector* > m_flxPntGhostExtraVars;

  /// flux point coordinates
  std::vector< RealVector > m_flxPntCoords;

  /// vector with pointers to all the states in flux points
  std::vector< Framework::State* > m_allSol;

  /// vector with pointers to all the extra variables in flux points
  std::vector< RealVector* > m_allExtraVars;

  /// RealVector variable for internal flux point solutions
  std::vector< RealVector* > m_flxPntIntRVSol;

  /// RealVector variable for ghost flux point solutions
  std::vector< RealVector* > m_flxPntGhostRVSol;

  /// Riemann fluxes in flux points
  std::vector< RealVector > m_flxPntRiemannFlux;

  /// gradients in the internal flux points
  std::vector< std::vector< RealVector* > > m_flxPntIntGrads;

  /// gradients in the ghost flux points
  std::vector< std::vector< RealVector* > > m_flxPntGhostGrads;

  /// pointers to gradients in the internal flux points
  std::vector< std::vector< RealVector* >* > m_flxPntIntGradPtrs;

  /// pointers to gradients in the ghost flux points
  std::vector< std::vector< RealVector* >* > m_flxPntGhostGradPtrs;

  /// backup for physical variable in flux points
  RealVector m_backupPhysVar;

  /// backup for ghost states in flux points
  std::vector< RealVector > m_backupGhostStates;

  /// number of flux points the Riemann flux has to be evaluated in
  std::vector< CFuint > m_nbrFlxPnts;

  /// helper term for gradient computation
  RealVector m_gradTerm;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// number of dimensions in the physical model
  CFuint m_dim;

  /// number of extra variables in the physical model
  CFuint m_nbrExtraVars;

  /// face mapped coordinates of a given set of points
  std::vector< RealVector > m_faceMappedCoordsPntSet;

  /// cell mapped coordinates of a given set of points
  std::vector< RealVector > m_cellMappedCoordsPntSet;

  /// solution reconstruction coefficients of a given set of points
  std::vector< std::vector< CFreal > > m_solRecCoefsPntSet;

  /// states in given point set
  std::vector< RealVector > m_statesPntSet;

  /// Internal solutions in given set of points
  std::vector< Framework::State* > m_intSolPntSet;

  /// Ghost solutions in given set of points
  std::vector< Framework::State* > m_ghostSolPntSet;

  /// Pointers to internal solutions in given set of points
  std::vector< RealVector* > m_intSolPntSetPtrs;

  /// Pointers to ghost solutions in given set of points
  std::vector< RealVector* > m_ghostSolPntSetPtrs;

  /// unit normal vector in points set
  std::vector< RealVector > m_unitNormalPntSet;

  /// face Jacobian vector sizes in points set
  std::vector< RealVector > m_faceJacobPntSet;

  /// absolute coordinates in points set
  std::vector< RealVector > m_coordsPntSet;

  /// gradients in given point set
  std::vector< std::vector< RealVector > > m_gradsPntSet;

  /// internal gradients in given point set
  std::vector< std::vector< RealVector* > > m_intGradsPntSet;

  /// ghost gradients in given point set
  std::vector< std::vector< RealVector* > > m_ghostGradsPntSet;
  
private:

  /// Physical data temporary vector
  RealVector m_pData;

}; // class BaseBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BaseBndFaceTermComputer_hh

