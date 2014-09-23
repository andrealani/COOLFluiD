#ifndef COOLFluiD_Numerics_SpectralFD_NSBR2FaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_NSBR2FaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/BR2FaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms for the
 * Navier-Stokes equations with the Bassi-Rebay II approach.
 *
 * @author Kris Van den Abeele
 */
class NSBR2FaceTermComputer : public BR2FaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NSBR2FaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  NSBR2FaceTermComputer(const std::string& name);

  /// Destructor
  ~NSBR2FaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSBR2FaceTermComputer";
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

}; // class NSBR2FaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NSBR2FaceTermComputer_hh
