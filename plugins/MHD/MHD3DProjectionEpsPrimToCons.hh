#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionEpsPrimToCons_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionEpsPrimToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {
      class MHDProjectionEpsTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive 
 * to conservative variables
 *
 * @author Andrea Lani
 */
class MHD3DProjectionEpsPrimToCons : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  MHD3DProjectionEpsPrimToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~MHD3DProjectionEpsPrimToCons();
 
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
  Common::SafePtr<MHDProjectionEpsTerm> _model;

}; // end of class MHD3DProjectionEpsPrimToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionEpsPrimToCons_hh
