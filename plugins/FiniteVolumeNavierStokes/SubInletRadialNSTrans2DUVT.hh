#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletRadialNSTrans2DUVT_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletRadialNSTrans2DUVT_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a radial subsonic inlet command with the initial conditions given for U,V,T and
  * the turbulent variables
  *
  * @author Milan Zaloudek
  *
  */

class SubInletRadialNSTrans2DUVT : public FVMCC_BC {
public:

  typedef Physics::NavierStokes::MultiScalarVarSet
  <Physics::NavierStokes::Euler2DVarSet> ConvTrans2DVarSet;

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletRadialNSTrans2DUVT(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletRadialNSTrans2DUVT();

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
  Common::SafePtr<ConvTrans2DVarSet> _varSetTurb;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// velocity magnitude
  CFreal _Vinf;

  /// static temperature
  CFreal _temperature;

  /// Turbulent Variables K, Omega...
  std::vector<CFreal> _turbVars;

}; // end of class SubInletRadialNSTrans2DUVT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletRadialNSTrans2DUVT_hh
