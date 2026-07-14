#ifndef COOLFluiD_Physics_EntropyMHD_EntropyMHD3DConsToPrim_hh
#define COOLFluiD_Physics_EntropyMHD_EntropyMHD3DConsToPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * Transformer from conservative to primitive variables
 * for EntropyMHD 3D MHD with GLM projection.
 *
 * @author Rayan Dhib
 */
class EntropyMHD3DConsToPrim : public Framework::VarSetTransformer {
public:

  EntropyMHD3DConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model);
  ~EntropyMHD3DConsToPrim();

  void transform(const Framework::State& state, Framework::State& result);
  void transformFromRef(const RealVector& data, Framework::State& result);

}; // end of class EntropyMHD3DConsToPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_EntropyMHD_EntropyMHD3DConsToPrim_hh
