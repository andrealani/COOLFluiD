#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb3DTtPtRadialAlpha_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb3DTtPtRadialAlpha_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
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
  *
  * @author Milan Zaloudek
  *
  */

class SubInletNSTurb3DTtPtRadialAlpha : public FVMCC_BC {
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
  SubInletNSTurb3DTtPtRadialAlpha(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletNSTurb3DTtPtRadialAlpha();

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

  /// total temperature
  CFreal     _tTotal;

  /// total pressure
  CFreal     _pTotal;

  /// switch of swirling velocity component
  bool       _isSwirling;

  /// swirling intensity 0->zero u-velocity, 1->only u-velocity
  CFreal     _swirl;

  /// values of turbulent vars
  std::vector<CFreal> _turbVars;
}; // end of class SubInletNSTurb3DTtPtRadialAlpha

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb3DTtPtRadialAlpha_hh
