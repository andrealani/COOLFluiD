#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorVelocity3DTurb_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorVelocity3DTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
      class NavierStokes3DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a 3D mirror BC for turbulent computations
   * making sure that turbulent variables are mirrored too - not true at MirrorEuler3D
   *
   * @author Milan Zaloudek
   *
   */
class MirrorVelocity3DTurb : public FVMCC_BC {
public:

  typedef Physics::NavierStokes::MultiScalarVarSet
  <Physics::NavierStokes::Euler3DVarSet> ConvTurb3DVarSet;
 
  typedef Physics::NavierStokes::NavierStokesTurbVarSet<
    Physics::NavierStokes::NavierStokes3DVarSet, 0> DiffTurb3DVarSet;
  
  /**
   * Constructor
   */
  MirrorVelocity3DTurb(const std::string& name);

  /**
   * Default destructor
   */
  ~MirrorVelocity3DTurb();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private: // data
  /// physical model var set
  Common::SafePtr<ConvTurb3DVarSet> _varSetTurb;

  /// physical model diffusive variable set
  Common::SafePtr<DiffTurb3DVarSet> _diffVarTurb;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// values of turbulent vars
  std::vector<CFreal> _turbVars;

}; // end of class MirrorVelocity3DTurb

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorVelocity3DTurb_hh
