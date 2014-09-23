#ifndef COOLFluiD_Numerics_SpectralFD_NSIPFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_NSIPFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/IPFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms for the
 * Navier-Stokes equations with a fully compact approach.
 *
 * @author Kris Van den Abeele
 */
class NSIPFaceTermComputer : public IPFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NSIPFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  NSIPFaceTermComputer(const std::string& name);

  /// Destructor
  ~NSIPFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSIPFaceTermComputer";
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

}; // class NSIPFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NSIPFaceTermComputer_hh
