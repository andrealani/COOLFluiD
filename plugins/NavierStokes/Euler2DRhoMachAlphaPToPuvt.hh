#ifndef COOLFluiD_Physics_NavierStokes_Euler2DRhoMachAlphaPToPuvt_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DRhoMachAlphaPToPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

      class EulerTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from RhoMachAlphaP
 * to conservative variables
 *
 * @author Thomas Wuilbaut
 *
 */
class Euler2DRhoMachAlphaPToPuvt : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DRhoMachAlphaPToPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DRhoMachAlphaPToPuvt();
  
  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
  /**
   * Set up
   */
  void setup(const CFuint maxNbTransStates);

private:

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler2DRhoMachAlphaPToPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DRhoMachAlphaPToPuvt_hh
