#ifndef COOLFluiD_Numerics_SpectralFD_NSIPBndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_NSIPBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/IPBndFaceTermComputer.hh"

// #include "SpectralFDNavierStokes/NSFaceDiffusiveFluxIPApproach.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms at the boundaries
 * for the Navier-Stokes equations if the scheme is fully compact.
 *
 * @author Kris Van den Abeele
 */
class NSIPBndFaceTermComputer : public IPBndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NSIPBndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  NSIPBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~NSIPBndFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSIPBndFaceTermComputer";
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

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

  /// NSFaceDiffusiveFluxIPApproach object
//   Common::SafePtr< NSFaceDiffusiveFluxIPApproach > m_nsFaceDiffusiveFluxIP;

  /// gradients in the internal flux points
  std::vector< std::vector< RealVector* > > m_flxPntIntGradVarGrads;

  /// gradients in the ghost flux points
  std::vector< std::vector< RealVector* > > m_flxPntGhostGradVarGrads;

}; // class NSIPBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NSIPBndFaceTermComputer_hh
