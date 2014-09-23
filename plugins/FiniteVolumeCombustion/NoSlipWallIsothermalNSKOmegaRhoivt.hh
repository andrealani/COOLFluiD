#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSKOmegaRhoivt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSKOmegaRhoivt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

	namespace Framework {
		class PhysicalChemicalLibrary;
	}

		namespace Numerics {

			namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the no slip wall bc
   * for K-Omega and SA turbulence models using Rho_i, u, v, T as variables
   *
   * @author Alessandro Mazzetti
   * 
   *
   */
template <class CVARSET, class DVARSET>
class NoSlipWallIsothermalNSKOmegaRhoivt : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallIsothermalNSKOmegaRhoivt(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalNSKOmegaRhoivt();

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
  
  //data
  /// pointer to the physical chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// interactive flag to run adiabatic
  bool m_adiabatic;

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
  
  
  // ghost and inner massTtmTt, massTtmS and RoverMtot are
  // divide "inner" and "ghost" just to be sure even if they 
  // should be equal
  /// sum of Moles of species over mass of mixture - inner
  CFreal m_massTtmTt_In;
  
  /// Rgas over total molar mass of mixture - inner
  CFreal m_RoverMtot_In;
  
  /// sum of Moles of species over mass of mixture - ghost
  CFreal m_massTtmTt_G;
  
  /// Rgas over total molar mass of mixture - ghost
  CFreal m_RoverMtot_G;
  
  /// mass species fractions in the internal cell
  RealVector m_ysIn;
  
  /// molar masses
  RealVector m_mmasses;
  
  ///  Moles of species over mass of mixture
  RealVector m_massTtmS_In;
  
  ///  Moles of species over mass of mixture
  RealVector m_massTtmS_G;
  
  ///  wall mass fractions
  RealVector m_ysW;
  
  ///  ghost mass fractions
  RealVector m_ysG;
  
  ///  inner densities
  RealVector innerRho;
  
}; // end of class NoSlipWallIsothermalNSKOmegaRhoivt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NoSlipWallIsothermalNSKOmegaRhoivt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSKOmegaRhoivt_hh
