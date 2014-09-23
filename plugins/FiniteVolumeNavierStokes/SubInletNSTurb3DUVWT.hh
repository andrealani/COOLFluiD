#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb3DUVWT_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb3DUVWT_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a subsonic inlet command with the initial conditions given for U,V,T and
  * the turbulent variables
  *
  * @author Thomas Wuilbaut
  * @author Joao Pinto
  * @author Milan Zaloudek
  *
  */

class SubInletNSTurb3DUVWT : public FVMCC_BC {
public:
  
  typedef Physics::NavierStokes::MultiScalarVarSet<
  Physics::NavierStokes::Euler3DPvt<Physics::NavierStokes::Euler3DVarSet> > ConvTurb3DVarSet;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletNSTurb3DUVWT(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletNSTurb3DUVWT();

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

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// x velocity
  CFreal                                   _uinf;

  /// y velocity
  CFreal                                   _vinf;

  /// z velocity
  CFreal                                   _winf;

  /// static temperature
  CFreal                                   _temperature;

  /// Turbulent Variables K, Omega...
  std::vector<CFreal>                      _turbVars;

}; // end of class SubInletNSTurb3DUVWT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletNSTurb3DUVWT_hh
