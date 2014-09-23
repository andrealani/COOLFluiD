#ifndef COOLFluiD_SpectralFD_TVBLimiterEuler3DSpectralFD_hh
#define COOLFluiD_SpectralFD_TVBLimiterEuler3DSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/ReconstructStatesSpectralFD.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

#include "SpectralFDNavierStokes/TVBLimiterEulerSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
    }
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that applies an elementwise TVB limiter to the solution for 3D Euler,
 * taking into account the requirement of pressure positivty (pressure is limited instead of rhoE)
 *
 * @author Kris Van den Abeele
 *
 */
class TVBLimiterEuler3DSpectralFD : public TVBLimiterEulerSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit TVBLimiterEuler3DSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~TVBLimiterEuler3DSpectralFD();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

private: // functions

  /// compute the pressure from the conserved variables
  void computePressFromConsVar(const RealVector& consVar, CFreal& press);

  /// compute the rhoE from pressure and the other conserved variables
  void computeRhoEFromPressAndOtherConsVar(RealVector& consVar, const CFreal& press);

protected: // data

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> m_eulerVarSet;

  /// heat capacity ratio minus one
  CFreal m_gammaMinusOne;

}; // class TVBLimiterEuler3DSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_TVBLimiterEuler3DSpectralFD_hh
