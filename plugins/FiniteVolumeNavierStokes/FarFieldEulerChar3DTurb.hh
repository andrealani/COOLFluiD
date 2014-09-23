#ifndef COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar3DTurb_hh
#define COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar3DTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/FarFieldEulerChar3D.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class NavierStokes3DVarSet;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to implement turbulent farfield boundary condition
 * using the characteristics
 *
 * @author Milan Zaloudek
 *
 */
class FarFieldEulerChar3DTurb : public FarFieldEulerChar3D {

public:

  typedef Physics::NavierStokes::MultiScalarVarSet
  <Physics::NavierStokes::Euler3DVarSet> ConvTurb3DVarSet;

  typedef Physics::NavierStokes::NavierStokesTurbVarSet<
    Physics::NavierStokes::NavierStokes3DVarSet, 0> DiffTurb3DVarSet;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FarFieldEulerChar3DTurb(const std::string& name);

  /**
   * Default destructor
   */
  ~FarFieldEulerChar3DTurb();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

private:

  // physical model var set
  Common::SafePtr<ConvTurb3DVarSet> _varSetTurb;

 /// physical model diffusive variable set
  Common::SafePtr<DiffTurb3DVarSet> _diffVarTurb;

  // value of turbulent vars at farfield
  std::vector<CFreal> _turbVars;

}; // end of class FarFieldEulerChar3DTurb

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar3DTurb_hh
