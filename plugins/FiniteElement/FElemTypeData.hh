#ifndef COOLFluiD_Numerics_FiniteElement_FElemTypeData_hh
#define COOLFluiD_Numerics_FiniteElement_FElemTypeData_hh

//////////////////////////////////////////////////////////////////////////////


#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Common/Quartet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class BlockAccumulator; }

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

typedef Common::Quartet<Framework::BlockAccumulator*,
                       RealMatrix*,
                       RealVector*,
                       std::vector<RealVector>*> FElemTypeData;

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_FElemTypeData_hh




