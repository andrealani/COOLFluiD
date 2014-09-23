#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSGReKOPvt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSGReKOPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the no slip wall bc
   * for K-Omega and SA turbulence models
   *
   * @author Khalil Bensassi
   *
   */
template <class CVARSET, class DVARSET>
class NoSlipWallIsothermalNSGReKOPvt : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallIsothermalNSGReKOPvt(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalNSGReKOPvt();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);


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
   * Reposition the node if any of the given ghostState values is < 0
   */
  void repositionNode(const CFreal& innerValue, CFreal& ghostValue);

  /**
   * This function makes a linear interpolation between the values in the
   * inner state and the ghost state ones
   */
  void linearInterpolate(const CFreal& innerValue, const CFreal& wallValue, CFreal& ghostValue)
  {
    ghostValue = innerValue - (innerValue - wallValue)*(m_drXiXg/m_drXiXw);
  }

 private:

  /// physical model convective variable set
  Common::SafePtr<CVARSET> _varSetTurb;

  /// physical model diffusive variable set
  Common::SafePtr<DVARSET> _diffVarTurb;

  /// X-component of a velocity vector of the wall
  CFreal _xWallVelocity;

  /// Y-component of a velocity vector of the wall
  CFreal _yWallVelocity;

  /// Z-component of a velocity vector of the wall
  CFreal _zWallVelocity;
  
  /// turb intensity variable ID
  CFuint m_kID;
  
  /// turb intensity in the ghost state
  CFreal m_ghostK;

  /// turb intensity at the wall
  CFreal m_wallK;
  
  /// temperature at the wall
  CFreal m_wallT;
  
  /// minimum turb intensity at the ghost
  CFreal m_ghostKMin;
  
  /// turb intensity in the ghost state
  CFreal m_ghostGa;
  
  /// turb intensity in the ghost state
  CFreal m_ghostRe;

}; // end of class NoSlipWallIsothermalNSGReKOPvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NoSlipWallIsothermalNSGReKOPvt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSGReKOPvt_hh
