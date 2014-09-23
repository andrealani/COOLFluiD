#ifndef COOLFluiD_Numerics_FiniteVolume_TransOutletNSTurb2D_hh
#define COOLFluiD_Numerics_FiniteVolume_TransOutletNSTurb2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//   namespace Physics {
//     namespace NavierStokes {
//       class ConvTurb2DVarSet;
//     }
//   }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an outlet command for partly sub and partly supersonic flow
 *
 * @author Milan Zaloudek
 *
 */
class TransOutletNSTurb2D : public FVMCC_BC {
public:

  typedef Physics::NavierStokes::MultiScalarVarSet
  <Physics::NavierStokes::Euler2DVarSet> ConvTurb2DVarSet;

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  TransOutletNSTurb2D(const std::string& name);

  /**
   * Default destructor
   */
  ~TransOutletNSTurb2D();

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

  /// physical model (in conservative variables)
  Common::SafePtr<ConvTurb2DVarSet> _varSetTurb;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// static pressure
  CFreal _pressure;

}; // end of class TransOutletEuler2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_TransOutletNSTurb2D_hh
