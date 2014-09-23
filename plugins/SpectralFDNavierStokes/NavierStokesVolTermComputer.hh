#ifndef COOLFluiD_Numerics_SpectralFD_NavierStokesVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_NavierStokesVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/BaseVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the volume terms in a cell for the Navier-Stokes equations
 *
 * @author Kris Van den Abeele
 */
class NavierStokesVolTermComputer : public BaseVolTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NavierStokesVolTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  NavierStokesVolTermComputer(const std::string& name);

  /// Destructor
  ~NavierStokesVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NavierStokesVolTermComputer";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /**
   * volume term contribution to the gradients
   * @pre reconstructStates(), setVolumeTermData() and setCellFaceNormTransfM()
   */
  virtual void computeGradientVolumeTerm(std::vector< std::vector< RealVector > >& gradUpdates);

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

}; // class NavierStokesVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NavierStokesVolTermComputer_hh

