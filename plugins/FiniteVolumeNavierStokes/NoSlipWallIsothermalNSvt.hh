#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSvt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/MirrorVelocity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class MapGeoToTrsAndIdx;
  }

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the Isothermal wall
   * BC in case of update variables including velocity components and temperature
   * [p v T y_i] or [rho_i v T] etc.
   *
   * @author Andrea Lani
   *
   */
template <class MODEL>
class NoSlipWallIsothermalNSvt : public MirrorVelocity {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallIsothermalNSvt(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NoSlipWallIsothermalNSvt();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Apply boundary condition on the given ghost state
   * @param innerState state vector associated to the internal cell
   * @param ghostState ghost state vector
   * @param ghostT  temperature in the ghost state
   */
  virtual void setGhostStateImpl(const Framework::State& innerState,
				 Framework::State& ghostState);

  /**
   * This function tells if there any negative ghost value
   */
  bool isAnyNegativeGhostValue(const RealVector& values) const
  {
    for (CFuint i = 0; i < values.size(); ++i) {
      if (values[i] < 0.0) return true;
    }
    return false;
  }

  /**
   * Set the wall temperature
   */
  virtual void setWallTemperature(Framework::GeometricEntity *const face);

  /**
   * Get the convective heat flux corresponding to the given wall temperature
   */
  CFreal convectiveHeatFlux(const Framework::State& innerState,
			    Framework::State& ghostState,
			    const CFreal Twall);

protected:

  /**
   * Impose adiabatic conditions
   */
  void imposeAdiabatic(const Framework::State& innerState,
		       Framework::State& ghostState);

  /**
   * Compute the wall temperature from radiative equilibrium
   */
  void computeTwallFromRadEq(Framework::GeometricEntity *const face);

  /**
   * Compute the wall temperature from file equilibrium
   */
  void computeTwallFromFile(Framework::GeometricEntity *const face);

protected:

  /// pointer to the mapping face - TRS
  Common::SafePtr<Framework::MapGeoToTrsAndIdx> m_mapGeoToTrs;
  
  /// diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffVar;

  /// map each wall faceID to the corresponding wall temperature
  Common::CFMap<Framework::TopologicalRegionSet*, RealVector*> m_mapTrs2Twall;

  /// array of gradients
  std::vector<RealVector*> m_gradients;

  // array of average values (p, u, v, T, ...)
  RealVector m_avValues;

  /// temperature at the wall
  CFreal m_wallTempIn;

  /// temperature variable ID
  CFuint m_tempID;
  
  /// roto-translational temperature in the ghost state
  CFreal m_ghostT;
  
  // array of states involved in the flux computation
  std::vector<RealVector*> m_gradStates;

  // array of values (p, u, v, T, ...)
  RealMatrix m_gradValues;
  
  /// velocity at the wall
  std::vector<CFreal> m_wallVelocity;
  
  /// temperature at the wall
  CFreal m_wallTemp;

  /// minimum temperature at the ghost
  CFreal m_ghostTempMin;
  
  /// flag to activate the radiative equilibrium
  bool m_radiativeEq;
  
  /// flag to fix the wall temperature directly in the ghost cell
  bool m_fixTWallInGhost;
  
  /// wall emissivity
  CFreal m_emissivity;
  
  /// max variation of wall temperature between two time steps
  CFreal m_maxRadEqDT;
  
  /// distant body temperature
  CFreal m_dbTemp;
  
  /// interactive flag to run adiabatic
  bool m_adiabatic;
  
  /// name of the file where the temperature distribution is provided
  std::string m_fileNameTw;
  
}; // end of class NoSlipWallIsothermalNSvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NoSlipWallIsothermalNSvt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSvt_hh
