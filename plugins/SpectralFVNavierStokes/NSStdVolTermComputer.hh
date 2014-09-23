#ifndef COOLFluiD_Numerics_SpectralFV_NSStdVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_NSStdVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFV/StdVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the volume terms in a cell
 * in a standard way using Gaussian quadrature. This computer is specifically
 * for Navier-Stokes
 *
 * @author Kris Van den Abeele
 */
class NSStdVolTermComputer : public StdVolTermComputer {

public:  // methods

  /// Constructor
  NSStdVolTermComputer(const std::string& name);

  /// Destructor
  ~NSStdVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSStdVolTermComputer";
  }

  /// Set up private data and data
  void setup();

  /**
   * Unsetup private data
   */
  void unsetup();

  /**
   * volume term contribution to the gradients
   * @pre reconstructStates(), setVolumeTermData() and setCellFaceNormTransfM()
   */
  void computeGradientVolumeTerm(std::vector< std::vector< RealVector > >& gradUpdates);

protected: // data

  /// Navier-Stokes variable set
  /// @note it is assumed that the diffusive term is a Navier-Stokes term here,
  /// but it should not be a big problem to make it more general.
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

}; // class NSStdVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_NSStdVolTermComputer_hh

