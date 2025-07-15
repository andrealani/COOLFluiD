#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionConsToPrimE_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionConsToPrimE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {
      class MHDProjectionTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to primitive variables
 *
 * @author Andrea Lani
 * @author Radka Keslerova
 */
class MHD3DProjectionConsToPrimE : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  MHD3DProjectionConsToPrimE(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~MHD3DProjectionConsToPrimE();
 
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
  Common::SafePtr<MHDProjectionTerm> _model;

}; // end of class MHD3DProjectionConsToPrimE

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionConsToPrimE_hh
