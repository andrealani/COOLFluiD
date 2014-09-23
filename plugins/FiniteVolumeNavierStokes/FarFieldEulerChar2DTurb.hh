#ifndef COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar2DTurb_hh
#define COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar2DTurb_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/FarFieldEulerChar2D.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokes2DVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to implement farfield boundary condition in turbulent computations
 * using th characteristics
 *
 * @author Thomas Wuilbaut
 * @author Milan Zaloudek
 *
 */


class FarFieldEulerChar2DTurb : public FarFieldEulerChar2D {

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
  FarFieldEulerChar2DTurb(const std::string& name);

  /**
   * Default destructor
   */
  ~FarFieldEulerChar2DTurb();

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
  Common::SafePtr<ConvTurb2DVarSet> _varSetTurb;
  
  /// physical model diffusive variable set
  Common::SafePtr<DiffTurb2DVarSet> _diffVarTurb;
  
  // value of turbulent vars at farfield
  std::vector<CFreal> _turbVars;

}; // end of class FarFieldEulerChar2DTurb

//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar2DTurb_hh
