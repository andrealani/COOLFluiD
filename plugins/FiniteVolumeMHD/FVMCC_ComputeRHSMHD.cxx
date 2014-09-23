#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/FVMCC_ComputeRHSMHD.hh"
#include "FiniteVolume/FVMCC_ComputeRHS.hh"
#include "FiniteVolume/FVMCC_ComputeRhsJacob.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ComputeRHSMHD<FVMCC_ComputeRHS>,
		      CellCenterFVMData,
		      FiniteVolumeMHDModule>
fvmcc_computeRHSMHDProv("FVMCCMHD");

MethodCommandProvider<FVMCC_ComputeRHSMHD<FVMCC_ComputeRhsJacob>,
		      CellCenterFVMData,
		      FiniteVolumeMHDModule>
fvmcc_computeRhsJacobMHDProv("NumJacobMHD");

//////////////////////////////////////////////////////////////////////////////
 
} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
