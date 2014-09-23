#include "FiniteVolumeCUDA/FVMCC_ComputeRHSCell.hh"
#include "FiniteVolumeCUDA/FVMCC_ComputeRhsJacobCell.hh"
#include "FiniteVolumeMHD/FVMCC_ComputeRHSMHD.hh"

#include "FiniteVolumeCUDA/FiniteVolumeCUDA.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/VarSetListT.hh"
#include "FiniteVolume/LaxFriedFlux.hh"
#include "FiniteVolume/LeastSquareP1PolyRec2D.hh"
#include "FiniteVolume/LeastSquareP1PolyRec3D.hh"
#include "FiniteVolume/BarthJesp.hh"
#include "MHD/MHD2DProjectionConsT.hh"
#include "MHD/MHD3DProjectionConsT.hh"
#include "MHD/MHD2DProjectionPrimT.hh"
#include "MHD/MHD3DProjectionPrimT.hh"
#include "MHD/MHDProjectionPrimToConsT.hh"
#include "FiniteVolumeMHD/LaxFriedFluxTanaka.hh"
#include "MHD/MHD2DProjectionVarSet.hh"
#include "MHD/MHD3DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

#define FVMCC_COMPUTE_MHD_RHS_PROV(__dim__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<FVMCC_ComputeRHSMHD<FVMCC_ComputeRHSCell<LaxFriedFlux, \
							       VarSetListT<MHD##__dim__##__svars__##T, MHD##__dim__##__uvars__##T>, \
							       LeastSquareP1PolyRec##__dim__ , BarthJesp, __nbBThreads__> >, \
		      CellCenterFVMData, FiniteVolumeCUDAModule>	\
fvmcc_computeRhsMHD##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);

// 48 block threads (default)
FVMCC_COMPUTE_MHD_RHS_PROV(2D, ProjectionCons, ProjectionCons, 48, "LaxFriedMHD2DCons")
FVMCC_COMPUTE_MHD_RHS_PROV(3D, ProjectionCons, ProjectionCons, 48, "LaxFriedMHD3DCons")
FVMCC_COMPUTE_MHD_RHS_PROV(2D, ProjectionCons, ProjectionPrim, 48, "LaxFriedMHD2DPrim")
FVMCC_COMPUTE_MHD_RHS_PROV(3D, ProjectionCons, ProjectionPrim, 48, "LaxFriedMHD3DPrim")
#undef FVMCC_COMPUTE_MHD_RHS_PROV

#define FVMCC_COMPUTE_MHD_RHS_PROV_TANAKA(__dim__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<FVMCC_ComputeRHSMHD<FVMCC_ComputeRHSCell<LaxFriedFluxTanaka<MHD##__dim__##ProjectionVarSet>, \
							       VarSetListT<MHD##__dim__##__svars__##T, MHD##__dim__##__uvars__##T>, \
							       LeastSquareP1PolyRec##__dim__ , BarthJesp, __nbBThreads__> >, \
		      CellCenterFVMData, FiniteVolumeCUDAModule>	\
fvmcc_computeRhsMHDTanaka##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);

// 48 block threads (default)
FVMCC_COMPUTE_MHD_RHS_PROV_TANAKA(2D, ProjectionCons, ProjectionCons, 48, "LaxFriedTanakaMHD2DCons")
FVMCC_COMPUTE_MHD_RHS_PROV_TANAKA(3D, ProjectionCons, ProjectionCons, 48, "LaxFriedTanakaMHD3DCons")
FVMCC_COMPUTE_MHD_RHS_PROV_TANAKA(2D, ProjectionCons, ProjectionPrim, 48, "LaxFriedTanakaMHD2DPrim")
FVMCC_COMPUTE_MHD_RHS_PROV_TANAKA(3D, ProjectionCons, ProjectionPrim, 48, "LaxFriedTanakaMHD3DPrim")
#undef FVMCC_COMPUTE_MHD_RHS_PROV_TANAKA

//////////////////////////////////////////////////////////////////////////////

#define FVMCC_COMPUTE_MHD_RHS_JACOB_PROV(__dim__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<FVMCC_ComputeRHSMHD<FVMCC_ComputeRhsJacobCell<LaxFriedFlux, \
								    VarSetListT<MHD##__dim__##__svars__##T, MHD##__dim__##__uvars__##T>, \
								    LeastSquareP1PolyRec##__dim__ , BarthJesp, __nbBThreads__> >, \
		      CellCenterFVMData, FiniteVolumeCUDAModule>	\
fvmcc_computeRhsJacobMHD##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);

// 48 block threads (default)
FVMCC_COMPUTE_MHD_RHS_JACOB_PROV(2D, ProjectionCons, ProjectionCons, 48, "NumJacobLaxFriedMHD2DCons")
FVMCC_COMPUTE_MHD_RHS_JACOB_PROV(3D, ProjectionCons, ProjectionCons, 48, "NumJacobLaxFriedMHD3DCons")
FVMCC_COMPUTE_MHD_RHS_JACOB_PROV(2D, ProjectionCons, ProjectionPrim, 48, "NumJacobLaxFriedMHD2DPrim")
FVMCC_COMPUTE_MHD_RHS_JACOB_PROV(3D, ProjectionCons, ProjectionPrim, 48, "NumJacobLaxFriedMHD3DPrim")
#undef FVMCC_COMPUTE_MHD_RHS_JACOB_PROV

#define FVMCC_COMPUTE_MHD_RHS_JACOB_PROV_TANAKA(__dim__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<FVMCC_ComputeRHSMHD<FVMCC_ComputeRhsJacobCell<LaxFriedFluxTanaka<MHD##__dim__##ProjectionVarSet>, \
								    VarSetListT<MHD##__dim__##__svars__##T, MHD##__dim__##__uvars__##T>, \
								    LeastSquareP1PolyRec##__dim__ , BarthJesp, __nbBThreads__> >, \
		      CellCenterFVMData, FiniteVolumeCUDAModule>	\
fvmcc_computeRhsJacobMHDTanaka##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);

// 48 block threads (default)
FVMCC_COMPUTE_MHD_RHS_JACOB_PROV_TANAKA(2D, ProjectionCons, ProjectionCons, 48, "NumJacobLaxFriedTanakaMHD2DCons")
FVMCC_COMPUTE_MHD_RHS_JACOB_PROV_TANAKA(3D, ProjectionCons, ProjectionCons, 48, "NumJacobLaxFriedTanakaMHD3DCons")
FVMCC_COMPUTE_MHD_RHS_JACOB_PROV_TANAKA(2D, ProjectionCons, ProjectionPrim, 48, "NumJacobLaxFriedTanakaMHD2DPrim")
FVMCC_COMPUTE_MHD_RHS_JACOB_PROV_TANAKA(3D, ProjectionCons, ProjectionPrim, 48, "NumJacobLaxFriedTanakaMHD3DPrim")
#undef FVMCC_COMPUTE_MHD_RHS_JACOB_PROV_TANAKA
      
//////////////////////////////////////////////////////////////////////////////

   } // namespace FiniteVolume
    
  } // namespace Numerics

} // namespace COOLFluiD
