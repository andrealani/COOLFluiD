#ifndef COOLFluiD_Numerics_SpectralFV_BaseFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_BaseFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFV/FaceDiffusiveFlux.hh"
#include "SpectralFV/ReconstructStatesSpectralFV.hh"
#include "SpectralFV/RiemannFlux.hh"
#include "SpectralFV/SpectralFVElementData.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms
 *
 * @author Kris Van den Abeele
 */
class BaseFaceTermComputer : public SpectralFVMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFVMethodData,BaseFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  BaseFaceTermComputer(const std::string& name);

  /// Destructor
  ~BaseFaceTermComputer();

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
  void computeNeighbourCellData();

  /**
   * reconstruct solution in the flux points
   * @pre setFaceTermData()
   */
  void reconstructFluxPntsStates(const std::vector< Framework::State* >& cellLStates,
                                 const std::vector< Framework::State* >& cellRStates);

  /**
   * reconstruct face averaged solutions
   * @pre setFaceTermData() and reconstructFluxPntsStates()
   */
  virtual void reconstructFaceAvgState(const std::vector< Framework::State* >& cellLStates,
                                       const std::vector< Framework::State* >& cellRStates) = 0;

  /**
   * reconstruct gradients in the flux points
   * @pre setFaceTermData()
   */
  void reconstructFluxPntsGradients(const std::vector< std::vector< RealVector >* >& cellLGrads,
                                    const std::vector< std::vector< RealVector >* >& cellRGrads,
                                    const CFuint nbrCellLGrads, const CFuint nbrCellRGrads);

  /**
   * reconstruct gradients in the flux points, on one side
   * @pre setFaceTermData()
   */
  void reconstructFluxPntsGradients(const CFuint side,
                                    const std::vector< std::vector< RealVector >* >& cellGrads,
                                    const CFuint nbrCellGrads);

  /**
   * backup and reconstruct physical variable in the left cell in the required points
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
  void computeConvFaceTerm(RealVector& resUpdates);

  /**
   * compute the convective face term and the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeConvFaceTermAndWaveSpeedUpdates(RealVector& resUpdates,
                                              CFreal& waveSpeedUpdL,CFreal& waveSpeedUpdR);

  /**
   * face term contribution to the gradients
   * @pre reconstructFluxPntsStates(), setFaceTermData() and set the geometrical data of the face
   */
  void computeGradientFaceTerm(std::vector< std::vector< RealVector > >& gradUpdates);

  /**
   * compute the diffusive face term for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeDiffFaceTerm(RealVector& resUpdates);

  /**
   * compute the diffusive face term and the contribution to the update coefficient for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeDiffFaceTermAndUpdateCoefContributions
                                                      (RealVector& resUpdates,
                                                       CFreal& updateCoefContrL, CFreal& updateCoefContrR);

protected: // functions

  /**
   * compute the flux integrals from the fluxes in the flux points
   * @pre fluxes in flux points are computed
   */
  virtual void computeFaceFluxIntegralFromFlxPntFluxes(RealVector& resUpdates) = 0;

  /**
   * compute face term contribution to the gradients from the solution in the flux points
   * @pre diffusiveVarSet->setGradientsVars()
   */
  virtual void computeGradFaceTermFromFlxPntSol(std::vector< std::vector< RealVector > >& gradUpdates) = 0;

protected: // data

  /// socket for face surfaces
  Framework::DataSocketSink<CFreal > socket_faceSurf;

  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;

  /// Diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;

  /// Strategy that reconstructs the states in a given number of nodes
  Common::SafePtr< ReconstructStatesSpectralFV > m_statesReconstr;

  /// Riemann flux
  Common::SafePtr< RiemannFlux > m_riemannFluxComputer;

  /// Face diffusive flux
  Common::SafePtr< FaceDiffusiveFlux > m_faceDiffFluxComputer;

  /// ratio between convective and diffusive cfl limit
  CFreal m_cflConvDiffRatio;

  /// flux point on SV faces reconstruction coefficients
  Common::SafePtr< std::vector< std::vector< std::vector< std::vector< CFreal > > > > >
      m_svFaceflxPntsRecCoefs;

  /// current face
  Framework::GeometricEntity* m_face;

  /// face orientation
  CFuint m_orient;

  /// face surface
  CFreal m_surf;

  /// cell volumes
  std::vector< CFreal > m_cellVolumes;

  /// average normal vector
  RealVector m_avgNormal;

  /// unit normal vector
  RealVector m_avgUnitNormal;

  /// unit normal vector in flux points
  std::vector< RealVector > m_unitNormalFlxPnt;

  /// flux point solutions in left and right cells
  std::vector< std::vector< Framework::State* > > m_flxPntSol;

  /// flux point averaged variables in left and right cells
  std::vector< std::vector< RealVector* > > m_flxPntExtraVars;

  /// solutions averaged over face
  std::vector< Framework::State* > m_faceAvgSol;

  /// extra variables averaged over face
  std::vector< RealVector* > m_faceAvgExtraVars;

  /// vector with pointers to all the states in flux points
  std::vector< Framework::State* > m_allSol;

  /// vector with pointers to all the extra variables in flux points
  std::vector< RealVector* > m_allExtraVars;

  /// RealVector variable for flux point solutions in left and right cells
  std::vector< std::vector< RealVector* > > m_flxPntRVSol;

  /// Riemann fluxes in flux points
  std::vector< RealVector > m_flxPntRiemannFlx;

  /// gradient variables in the flux points in left and right cells
  std::vector< std::vector< RealVector* > > m_flxPntGradVars;

  /// gradients in the flux points in left and right cells
  std::vector< std::vector< std::vector< RealVector* > > > m_flxPntGrads;

  /// pointers to gradients in the flux points in left and right cells
  std::vector< std::vector< std::vector< RealVector* >* > > m_flxPntGradPtrs;

  /// backup for physical variable in flux points
  std::vector< CFreal > m_backupPhysVar;

  /// number of flux points the Riemann flux has to be evaluated in
  std::vector< CFuint > m_nbrFlxPnts;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// number of extra variables in the physical model
  CFuint m_nbrExtraVars;

}; // class BaseFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_BaseFaceTermComputer_hh

