#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionEpsConsToPrim_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionEpsConsToPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {
      class MHDProjectionEpsTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to primitive variables
 *
 * @author Andrea Lani
 */
class MHD3DProjectionEpsConsToPrim : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  MHD3DProjectionEpsConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~MHD3DProjectionEpsConsToPrim();
 
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

}; // end of class MHD3DProjectionEpsConsToPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionEpsConsToPrim_hh
