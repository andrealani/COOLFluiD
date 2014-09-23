#include "Framework/NullEquationFilter.hh"
#include "CellCenterFVMData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NullEquationFilter<CellCenterFVMData>,
                       CellCenterFVMData,
                       EquationFilter<CellCenterFVMData>,
                       FiniteVolumeModule>
nullEquationFilterProv("Null");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
