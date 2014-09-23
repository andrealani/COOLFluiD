#ifndef COOLFluiD_Numerics_SpectralFD_BaseFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_BaseFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFD/FaceDiffusiveFlux.hh"
#include "SpectralFD/ReconstructStatesSpectralFD.hh"
#include "SpectralFD/RiemannFlux.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms.
 *
 * @author Kris Van den Abeele
 */
class BaseFaceTermComputer : public SpectralFDMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,BaseFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  BaseFaceTermComputer(const std::string& name);

  /// Destructor
  ~BaseFaceTermComputer();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BaseFaceTermComputer";
  }

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * set the orientation of the faces
   */
  void setFaceOrientation(const CFuint orient)
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
   * Set up private data and data
   */
  virtual void setup();

  /**
   * Unset up private data and data
   */
  virtual void unsetup();

  /**
   * set face term data
   */
  virtual void setFaceTermData();

  /**
   * compute face data
   * @pre setCurrentFace
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
  void reconstructFluxPntsStates(const std::vector< std::vector< Framework::State* >* >& cellStates, bool onlyExtraVars = false);

  /**
   * reconstruct gradients in the flux points
   * @pre setFaceTermData()
   */
  void reconstructFluxPntsGradients(const std::vector< std::vector< std::vector< RealVector >* > >& cellGrads);

  /**
   * reconstruct gradients in the flux points, on one side
   * @pre setFaceTermData()
   */
  void reconstructFluxPntsGradients(const CFuint side,
                                    const std::vector< std::vector< RealVector >* >& cellGrads);

  /**
   * backup and reconstruct physical variable in the given cell in the required points
   * @pre setVolumeTermData()
   */
  void backupAndReconstructPhysVar(const CFuint side, const CFuint iVar, const std::vector< Framework::State* >& cellStates);

  /**
   * backup physical variable in one cell in the required points
   */
  void backupPhysVar(const CFuint side, const CFuint iVar);

  /**
   * restore physical variable in one cell in the required points
   */
  void restorePhysVar(const CFuint side, const CFuint iVar);

  /**
   * compute the convective face term for this face
   * @pre reconstructFluxPntsStates(), setFaceTermData() and set the geometrical data of the face
   */
  void computeConvFaceTerm(std::vector< RealVector >& resUpdates);

  /**
   * compute the convective face term and the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeConvFaceTermAndWaveSpeedUpdates(std::vector< RealVector >& resUpdates,
                                              std::vector< CFreal >& waveSpeedUpd);

  /**
   * face term contribution to the gradients
   * @pre reconstructFluxPntsStates(), setFaceTermData() and set the geometrical data of the face
   */
  void computeGradientFaceTerm(std::vector< std::vector< std::vector< RealVector > > >& gradUpdates);

  /**
   * face term contribution to the gradients of the extra variables
   * @pre reconstructFluxPntsStates(), setFaceTermData() and set the geometrical data of the face
   */
  void computeGradientExtraVarsFaceTerm(std::vector< std::vector< std::vector< RealVector > > >& gradUpdates);

  /**
   * compute the diffusive face term for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeDiffFaceTerm(std::vector< RealVector >& resUpdates);

  /**
   * compute the diffusive face term and the contribution to the update coefficient for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeDiffFaceTermAndUpdateCoefContributions(std::vector< RealVector >& resUpdates,
                                                            std::vector< CFreal >& updateCoefContr);

protected: // functions

  /**
   * compute the residual updates from the fluxes in the flux points
   * @pre fluxes in flux points are computed
   */
  void computeFaceTermFromFlxPntFluxes(std::vector< RealVector >& resUpdates);

  /**
   * compute face term contribution to the gradients from the solution in the flux points
   */
  void computeGradFaceTermFromFlxPntSol(std::vector< std::vector< std::vector< RealVector > > >& gradUpdates);

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
  Common::SafePtr< ReconstructStatesSpectralFD > m_statesReconstr;

  /// Riemann flux
  Common::SafePtr< RiemannFlux > m_riemannFluxComputer;

  /// Face diffusive flux
  Common::SafePtr< FaceDiffusiveFlux > m_faceDiffFluxComputer;

  /// ratio between convective and diffusive cfl limit
  CFreal m_cflConvDiffRatio;

  /// reconstruction coefficients for the flux points
  Common::SafePtr< RealMatrix > m_flxPntsRecCoefs;

  /// derivation coefficients for the solution points
  Common::SafePtr< RealMatrix > m_solPntsDerivCoefs;

  /// local cell face - flux point connectivity per face connection orientation
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_faceFlxPntConn;

  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< std::vector< CFint > > > m_faceMappedCoordDir;

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

  /// face flux points face local coordinates
  Common::SafePtr< std::vector< RealVector > > m_faceFlxPntsFaceLocalCoords;

  /// local cell face - flux point cell mapped coordinate per face connection orientation
  Common::SafePtr< std::vector< std::vector< std::vector< RealVector > > > > m_faceFlxPntCellMappedCoords;

  /// current face
  Framework::GeometricEntity* m_face;

  /// face orientation
  CFuint m_orient;

  /// cell volumes
  std::vector< CFreal > m_cellVolumes;

  /// unit normal vector in flux points
  std::vector< RealVector > m_unitNormalFlxPnts;

  /// face Jacobian vector sizes
  std::vector< std::vector< CFreal > > m_faceJacobVecSizeFlxPnts;

  /// face Jacobian vector sizes (abs)
  std::vector< CFreal > m_faceJacobVecAbsSizeFlxPnts;

  /// face inverse characteristic lengths
  std::vector< CFreal > m_faceInvCharLengths;

  /// left and right cell extra variables
  std::vector< std::vector< RealVector* > > m_cellExtraVars;

  /// flux point solutions in left and right cells
  std::vector< std::vector< Framework::State* > > m_flxPntSol;

  /// vector with pointers to all the states in flux points
  std::vector< Framework::State* > m_allSol;

  /// flux point extra variables in left and right cells
  std::vector< std::vector< RealVector* > > m_flxPntExtraVars;

  /// vector with pointers to all the extra variables in flux points
  std::vector< RealVector* > m_allExtraVars;

  /// RealVector variable for flux point solutions in left and right cells
  std::vector< std::vector< RealVector* > > m_flxPntRVSol;

  /// Riemann flux numerical damping in flux points
  std::vector< RealVector > m_flxPntRiemannFlux;

  /// gradients in the flux points in left and right cells
  std::vector< std::vector< std::vector< RealVector* > > > m_flxPntGrads;

  /// pointers to gradients in the flux points in left and right cells
  std::vector< std::vector< std::vector< RealVector* >* > > m_flxPntGradPtrs;

  /// backup for physical variable in flux points
  std::vector< CFreal > m_backupPhysVar;

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

private:

  /// Physical data temporary vector
  RealVector m_pData;

}; // class BaseFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BaseFaceTermComputer_hh

