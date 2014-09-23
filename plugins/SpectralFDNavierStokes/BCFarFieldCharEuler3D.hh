#ifndef COOLFluiD_Numerics_SpectralFD_BCFarFieldCharEuler3D_hh
#define COOLFluiD_Numerics_SpectralFD_BCFarFieldCharEuler3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

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
 * This class represents a characteristic far field boundary condition for the 3D Euler/Navier-Stokes equations
 *
 * @author Kris Van den Abeele
 */
class BCFarFieldCharEuler3D : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCFarFieldCharEuler3D(const std::string& name);

  /// Destructor
  ~BCFarFieldCharEuler3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCFarFieldCharEuler3D";
  }

  /// Set up private data and data
  void setup();

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                             std::vector< std::vector< RealVector* > >& ghostGrads,
                             const std::vector< RealVector >& normals,
                             const std::vector< RealVector >& coords);

protected: // data

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// static temperature
  CFreal m_tempInf;

  /// static pressure
  CFreal m_presInf;

  /// x velocity
  CFreal m_uInf;

  /// y velocity
  CFreal m_vInf;

  /// z velocity
  CFreal m_wInf;

  /// density
  CFreal m_rhoInf;

  /// speed of sound
  CFreal m_cInf;

}; // class BCFarFieldCharEuler3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BCFarFieldCharEuler3D_hh

