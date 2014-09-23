#ifndef COOLFluiD_Numerics_SpectralFD_NSCompactVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_NSCompactVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/CompactVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the volume terms in a cell for the Navier-Stokes equations,
 * if the scheme is fully compact.
 *
 * @author Kris Van den Abeele
 */
class NSCompactVolTermComputer : public CompactVolTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NSCompactVolTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  NSCompactVolTermComputer(const std::string& name);

  /// Destructor
  ~NSCompactVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSCompactVolTermComputer";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /**
   * compute the diffusive volume term for this cell
   * @pre reconstructFaceStates(), setVolumeTermData(), setCellFaceNormTransfM() and
   *      setCellFaceTermData()
   */
  void computeCellDiffVolumeTerm(RealVector& resUpdates);

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

}; // class NSCompactVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NSCompactVolTermComputer_hh

