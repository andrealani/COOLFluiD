#ifndef COOLFluiD_Numerics_SpectralFD_NavierStokesFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_NavierStokesFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/BaseFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms for the Navier-Stokes equations.
 *
 * @author Kris Van den Abeele
 */
class NavierStokesFaceTermComputer : public BaseFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NavierStokesFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  NavierStokesFaceTermComputer(const std::string& name);

  /// Destructor
  ~NavierStokesFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NavierStokesFaceTermComputer";
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
  virtual void computeDiffFaceTermAndUpdateCoefContributions(std::vector< RealVector >& resUpdates,
                                                             std::vector< CFreal >& updateCoefContr);

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

}; // class NavierStokesFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NavierStokesFaceTermComputer_hh
