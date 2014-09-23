#ifndef COOLFluiD_Numerics_SpectralFD_NavierStokesBndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_NavierStokesBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms at the boundaries for the Navier-Stokes equations.
 *
 * @author Kris Van den Abeele
 */
class NavierStokesBndFaceTermComputer : public BaseBndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NavierStokesBndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  NavierStokesBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~NavierStokesBndFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NavierStokesBndFaceTermComputer";
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
   * compute the diffusive face term and the contribution to the update coefficient for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,CFreal& updateCoefContr);

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

}; // class NavierStokesBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NavierStokesBndFaceTermComputer_hh

