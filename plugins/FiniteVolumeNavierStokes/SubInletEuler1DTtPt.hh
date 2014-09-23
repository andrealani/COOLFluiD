#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEuler1DTtPt_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEuler1DTtPt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler1DVarSet;
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
   * @author Alessandro Munafò
   *
   */
class SubInletEuler1DTtPt : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEuler1DTtPt(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletEuler1DTtPt();

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
  Common::SafePtr<Physics::NavierStokes::Euler1DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// total temperature
  CFreal     _tTotal;

  /// total pressure
  CFreal     _pTotal;

}; // end of class SubInletEuler1DTtPt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEuler1DTtPt_hh
