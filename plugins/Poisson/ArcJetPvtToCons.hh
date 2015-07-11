#ifndef COOLFluiD_Physics_ArcJet_ArcJetPvtToCons_hh
#define COOLFluiD_Physics_ArcJet_ArcJetPvtToCons_hh

//////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace ArcJet {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from [p v T Bx By Bz Phi] to conservative variables
 *
 * @author Andrea Lani
 *
 */
template <typename BASE>
class ArcJetPvtToCons : public BASE {
public:
  
  /**
   * Default constructor without arguments
   */
  ArcJetPvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~ArcJetPvtToCons();

  /**
   * Transform a set of state vectors into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
  
}; // end of class ArcJetPvtToCons

//////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#include "ArcJetPvtToCons.ci"

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ArcJet_ArcJetPvtToCons_hh
