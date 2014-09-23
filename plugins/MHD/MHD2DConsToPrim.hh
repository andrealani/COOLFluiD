#ifndef COOLFluiD_Physics_MHD_MHD2DConsToPrim_hh
#define COOLFluiD_Physics_MHD_MHD2DConsToPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {
      class MHDTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to Roe to primitive variables
 *
 * @author Andrea Lani
 */
class MHD2DConsToPrim : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  MHD2DConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~MHD2DConsToPrim();
  
  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& pdata, Framework::State& result);
  
private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<MHDTerm> _model;

}; // end of class MHD2DConsToPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD2DConsToPrim_hh
