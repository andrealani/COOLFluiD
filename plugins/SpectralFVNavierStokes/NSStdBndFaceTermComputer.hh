#ifndef COOLFluiD_Numerics_SpectralFV_NSStdBndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_NSStdBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFV/StdBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the boundary face terms
 * in a standard way using Gaussian quadrature. This computer is specifically for Navier-Stokes.
 *
 * @author Kris Van den Abeele
 */
class NSStdBndFaceTermComputer : public StdBndFaceTermComputer {

public:  // methods

  /// Constructor
  NSStdBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~NSStdBndFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSStdBndFaceTermComputer";
  }

  /// Set up private data and data
  void setup();

  /// Unset up private data and data
  void unsetup();

  /**
   * compute the diffusive face term and the contribution to the update coefficient for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeDiffFaceTermAndUpdateCoefContributions(RealVector& resUpdates,CFreal& updateCoefContr);

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

}; // class NSStdBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_NSStdBndFaceTermComputer_hh
