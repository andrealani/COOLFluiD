#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNSTurb2D_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNSTurb2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class NavierStokes2DVarSet;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the no slip wall bc
   * for K-Omega and SA turbulence models
   *
   * @author Thomas Wuilbaut
   *
   */
class NoSlipWallAdiabaticNSTurb2D : public FVMCC_BC {

public:
  
  typedef Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DPuvt> ConvTurb2DVarSet;
  typedef Physics::NavierStokes::NavierStokesTurbVarSet<Physics::NavierStokes::NavierStokes2DVarSet, 0> DiffTurb2DVarSet;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallAdiabaticNSTurb2D(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallAdiabaticNSTurb2D();

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
  
private:

  /// physical model convective variable set
  Common::SafePtr<ConvTurb2DVarSet> _varSetTurb;

  /// physical model diffusive variable set
  Common::SafePtr<DiffTurb2DVarSet> _diffVarTurb;

  /// X-component of a velocity vector of the wall
  CFreal _xWallVelocity;

  /// Y-component of a velocity vector of the wall
  CFreal _yWallVelocity;

  /// turb intensity variable ID
  CFuint m_kID;
  
  /// turb intensity in the ghost state
  CFreal m_ghostK;
  
  /// turb intensity at the wall
  CFreal m_wallK;

  /// minimum turb intensity at the ghost
  CFreal m_ghostKMin;
  
}; // end of class NoSlipWallAdiabaticNSTurb2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNSTurb2D_hh
