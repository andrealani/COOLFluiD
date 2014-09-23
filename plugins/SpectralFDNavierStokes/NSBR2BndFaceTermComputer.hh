#ifndef COOLFluiD_Numerics_SpectralFD_NSBR2BndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_NSBR2BndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/BR2BndFaceTermComputer.hh"

// #include "SpectralFDNavierStokes/NSFaceDiffusiveFluxBR2Approach.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms at the boundaries
 * for the Navier-Stokes equations using the Bassi-Rebay II scheme.
 *
 * @author Kris Van den Abeele
 */
class NSBR2BndFaceTermComputer : public BR2BndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NSBR2BndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  NSBR2BndFaceTermComputer(const std::string& name);

  /// Destructor
  ~NSBR2BndFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSBR2BndFaceTermComputer";
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
   * compute the ghost gradients from the internal gradients, unit normals and coordinates
   * @pre internal gradients have to be reconstructed first
   */
  virtual void computeGhostGradients();

  /**
   * compute the diffusive face term for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(),
   *      setFaceTermData() and set the geometrical data of the face
   */
//   virtual void computeDiffFaceTerm(RealVector& resUpdates);

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
   * reconstruct gradients in point set
   * @pre setOutputData()
   */
  virtual std::vector< std::vector< RealVector > >& reconstructGivenPntsGrads(const std::vector< Framework::State* >& cellIntStates);

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

  /// NSFaceDiffusiveFluxIPApproach object
//   Common::SafePtr< NSFaceDiffusiveFluxBR2Approach > m_nsFaceDiffusiveFluxBR2;

  /// gradients in the internal flux points
  std::vector< std::vector< RealVector* > > m_flxPntIntGradVarGrads;

  /// gradients in the ghost flux points
  std::vector< std::vector< RealVector* > > m_flxPntGhostGradVarGrads;

  /// internal gradient variable gradients in a given point set
  std::vector< std::vector< RealVector* > > m_intGradVarGradsPntSet;

}; // class NSBR2BndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NSBR2BndFaceTermComputer_hh
