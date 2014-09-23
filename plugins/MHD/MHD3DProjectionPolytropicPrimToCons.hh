#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicPrimToCons_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicPrimToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {
      class MHDProjectionPolytropicTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive 
 * to conservative variables
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 */
class MHD3DProjectionPolytropicPrimToCons : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  MHD3DProjectionPolytropicPrimToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~MHD3DProjectionPolytropicPrimToCons();
 
  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<MHDProjectionPolytropicTerm> _model;

}; // end of class MHD3DProjectionPolytropicPrimToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicPrimToCons_hh
