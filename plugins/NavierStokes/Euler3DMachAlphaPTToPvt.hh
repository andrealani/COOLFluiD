#ifndef COOLFluiD_Physics_NavierStokes_Euler3DMachAlphaPTToPvt_hh
#define COOLFluiD_Physics_NavierStokes_Euler3DMachAlphaPTToPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

      class EulerTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from MachAlphaPT
 * to conservative variables
 *
 * @author Andrea Lani
 *
 */
class Euler3DMachAlphaPTToPvt : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler3DMachAlphaPTToPvt(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DMachAlphaPTToPvt();

  /**
   * Set up
   */
  void setup(const CFuint maxNbTransStates);

  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
private:

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler3DMachAlphaPTToPvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler3DMachAlphaPTToPvt_hh
