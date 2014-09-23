#ifndef COOLFluiD_DiscontGalerkin_DGElemTypeData_hh
#define COOLFluiD_DiscontGalerkin_DGElemTypeData_hh

//////////////////////////////////////////////////////////////////////////////


#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Common/Quartet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

    namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

typedef Common::Quartet<Framework::BlockAccumulator*,
                       RealMatrix*,
                       RealVector*,
                       std::vector<RealVector>*> DGElemTypeData;

//////////////////////////////////////////////////////////////////////////////

    } // namespace DiscontGalerkin

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_DiscontGalerkin_DGElemTypeData_hh




