#ifndef COOLFluiD_Physics_EntropyMHD_EntropyMHD3DPrimToCons_hh
#define COOLFluiD_Physics_EntropyMHD_EntropyMHD3DPrimToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * Transformer from primitive to conservative variables
 * for EntropyMHD 3D MHD with GLM projection.
 *
 * @author Rayan Dhib
 */
class EntropyMHD3DPrimToCons : public Framework::VarSetTransformer {
public:

  EntropyMHD3DPrimToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);
  ~EntropyMHD3DPrimToCons();

  void transform(const Framework::State& state, Framework::State& result);
  void transformFromRef(const RealVector& data, Framework::State& result);

}; // end of class EntropyMHD3DPrimToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_EntropyMHD_EntropyMHD3DPrimToCons_hh
