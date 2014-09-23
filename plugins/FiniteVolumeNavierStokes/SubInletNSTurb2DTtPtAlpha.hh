#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb2DTtPtAlpha_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb2DTtPtAlpha_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
      class NavierStokes2DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command with the initial conditions given for tTotal, pTotal and alpha
   *
   * @author Mehmet Sarp Yalim
   * @author Andrea Lani
   *
   */
class SubInletNSTurb2DTtPtAlpha : public FVMCC_BC {
public:

  typedef Physics::NavierStokes::MultiScalarVarSet
  <Physics::NavierStokes::Euler2DVarSet> ConvTurb2DVarSet;

  typedef Physics::NavierStokes::NavierStokesTurbVarSet<
    Physics::NavierStokes::NavierStokes2DVarSet, 0> DiffTurb2DVarSet;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletNSTurb2DTtPtAlpha(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletNSTurb2DTtPtAlpha();

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
  Common::SafePtr<ConvTurb2DVarSet> _varSetTurb;

  /// physical model diffusive variable set
  Common::SafePtr<DiffTurb2DVarSet> _diffVarTurb;

//  /// physical model var set
//  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSetTurb;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// total temperature
  CFreal     _tTotal;

  /// total pressure
  CFreal     _pTotal;

  /// alpha
  CFreal     _alpha;

  /// values of turbulent vars
  std::vector<CFreal> _turbVars;

}; // end of class SubInletNSTurb2DTtPtAlpha

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb2DTtPtAlpha_hh
