#ifndef COOLFluiD_Numerics_SpectralFV_NSQuadFreeVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_NSQuadFreeVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFV/QuadFreeVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the volume terms in a cell
 * in a quadrature-free way. This computer is specifically
 * for Navier-Stokes
 *
 * @author Kris Van den Abeele
 */
class NSQuadFreeVolTermComputer : public QuadFreeVolTermComputer {

public:  // methods

  /// Constructor
  NSQuadFreeVolTermComputer(const std::string& name);

  /// Destructor
  ~NSQuadFreeVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSQuadFreeVolTermComputer";
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

}; // class NSQuadFreeVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_NSQuadFreeVolTermComputer_hh

