#include "FiniteVolumeCUDA/FVMCC_ComputeRhsJacobCell.hh"
#include "Framework/MeshData.hh"
#include "Framework/BlockAccumulatorBaseCUDA.hh"
#include "Framework/CellConn.hh"
#include "Config/ConfigOptionPtr.hh"
#include "Common/CUDA/CudaDeviceManager.hh"
#include "Common/CUDA/CFVec.hh"
#include "Common/CUDA/CudaTimer.hh"
#include "FiniteVolume/CellData.hh"

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

#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfConsT.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfRhoiViTiT.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfRhoiViTiToConsT.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfConsToRhoiViTiT.hh"
#include "FiniteVolumeMultiFluidMHD/AUSMPlusUpFluxMultiFluid.hh"
#include "FiniteVolumeMultiFluidMHD/AUSMFluxMultiFluid.hh"

#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell2DProjectionConsT.hh"
#include "FiniteVolumeMaxwell/StegerWarmingMaxwellProjection2D.hh"



//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Physics::MHD;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

#define FVMCC_MHD_RHS_JACOB_PROV(__dim__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<FVMCC_ComputeRhsJacobCell<LaxFriedFlux, \
						VarSetListT<MHD##__dim__##__svars__##T, MHD##__dim__##__uvars__##T>, \
						LeastSquareP1PolyRec##__dim__ , BarthJesp, __nbBThreads__>, \
		      CellCenterFVMData, FiniteVolumeCUDAModule>	\
fvmcc_RhsJacobMHD##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);

// 48 block threads (default)
FVMCC_MHD_RHS_JACOB_PROV(2D, ProjectionCons, ProjectionCons, 48, "CellNumJacobLaxFriedMHD2DCons")
FVMCC_MHD_RHS_JACOB_PROV(3D, ProjectionCons, ProjectionCons, 48, "CellNumJacobLaxFriedMHD3DCons")
FVMCC_MHD_RHS_JACOB_PROV(2D, ProjectionCons, ProjectionPrim, 48, "CellNumJacobLaxFriedMHD2DPrim")
FVMCC_MHD_RHS_JACOB_PROV(3D, ProjectionCons, ProjectionPrim, 48, "CellNumJacobLaxFriedMHD3DPrim")
#undef FVMCC_MHD_RHS_JACOB_PROV

#define FVMCC_MHD_RHS_JACOB_PROV_TANAKA(__dim__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<FVMCC_ComputeRhsJacobCell<LaxFriedFluxTanaka<MHD##__dim__##ProjectionVarSet>, \
						VarSetListT<MHD##__dim__##__svars__##T, MHD##__dim__##__uvars__##T>, \
						LeastSquareP1PolyRec##__dim__ , BarthJesp, __nbBThreads__>, \
		      CellCenterFVMData, FiniteVolumeCUDAModule>	\
fvmcc_RhsJacobMHDTanaka##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);

// 48 block threads (default)
FVMCC_MHD_RHS_JACOB_PROV_TANAKA(2D, ProjectionCons, ProjectionCons, 48, "CellNumJacobLaxFriedTanakaMHD2DCons")
FVMCC_MHD_RHS_JACOB_PROV_TANAKA(3D, ProjectionCons, ProjectionCons, 48, "CellNumJacobLaxFriedTanakaMHD3DCons")
FVMCC_MHD_RHS_JACOB_PROV_TANAKA(2D, ProjectionCons, ProjectionPrim, 48, "CellNumJacobLaxFriedTanakaMHD2DPrim")
FVMCC_MHD_RHS_JACOB_PROV_TANAKA(3D, ProjectionCons, ProjectionPrim, 48, "CellNumJacobLaxFriedTanakaMHD3DPrim")
#undef FVMCC_MHD_RHS_JACOB_PROV_TANAKA




//Provider for StegerWarmingProjectionMaxwell2D
#define FVMCC_MAXWELL_RHS_JACOB_PROV_STEGER(__dim__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<FVMCC_ComputeRhsJacobCell<StegerWarmingMaxwellProjection2D<Maxwell##__dim__##ProjectionVarSet>, \
						VarSetListT<Maxwell##__dim__##__svars__##T, Maxwell##__dim__##__uvars__##T>, \
						LeastSquareP1PolyRec##__dim__ , BarthJesp, __nbBThreads__>, \
		      CellCenterFVMData, FiniteVolumeCUDAModule>	\
fvmcc_RhsJacobMaxwellSteger##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);

// 48 block threads (default)
FVMCC_MAXWELL_RHS_JACOB_PROV_STEGER(2D, ProjectionCons, ProjectionCons, 48, "CellNumJacobStegerWarmingMaxwell2DCons")

#undef FVMCC_MAXWELL_RHS_JACOB_PROV_STEGER


// Provider for AUSMPlusUpFlux 
#define FVMCC_MULTIFLUIDMHD_RHS_JACOB_PROV_AUSMPLUSUP(__dim__,__half__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<FVMCC_ComputeRhsJacobCell<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell##__dim__##ProjectionVarSet> >, \
			              VarSetListT<EulerMFMHD##__dim__##__half__##__svars__##T, EulerMFMHD##__dim__##__half__##__uvars__##T>, \
				      LeastSquareP1PolyRec##__dim__ , BarthJesp, __nbBThreads__>, \
		      CellCenterFVMData, FiniteVolumeCUDAModule>	\
fvmcc_RhsJacobMultiFluidMHDAUSMPlusUp##__dim__##__half__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);

// 48 block threads (default)
FVMCC_MULTIFLUIDMHD_RHS_JACOB_PROV_AUSMPLUSUP(2D, Half,Cons, RhoiViTi, 48, "CellNumJacobAUSMPlusUpEulerMFMHD2DHalfRhoiViTi")

#undef FVMCC_MULTIFLUIDMHD_RHS_JACOB_PROV_AUSMPLUSUP




//////////////////////////////////////////////////////////////////////////////

template <typename PHYS>
HOST_DEVICE inline void setState(CFreal* state, CFreal* statePtr, 
				 CFreal* node, CFreal* nodePtr)
{
  // copy the state node data to shared memory
  for (CFuint i = 0; i < PHYS::DIM; ++i) {node[i] = nodePtr[i];}
  // copy the state data to shared memory
  for (CFuint i = 0; i < PHYS::NBEQS; ++i) {state[i] = statePtr[i];} 
}
      
//////////////////////////////////////////////////////////////////////////////
      
template <typename PHYS>
HOST_DEVICE inline void setFaceNormal(FluxData<PHYS>* fd, CFreal* normal)
{  
  CudaEnv::CFVecSlice<CFreal,PHYS::DIM> n(normal);
  const CFreal area = n.norm2();
  fd->setFaceArea(area);
  const CFreal ovArea = 1./area;
  CudaEnv::CFVecSlice<CFreal,PHYS::DIM> un(fd->getUnitNormal());
  for (CFuint i = 0; i < PHYS::DIM; ++i) {
    un[i] = n[i]*ovArea;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

template <typename PHYS, typename PTR>
HOST_DEVICE void setFluxData(const CFuint f, const CFint stype, 
			     const CFuint stateID, const CFuint cellID, 
			     KernelData<CFreal>* kd, FluxData<PHYS>* fd,
			     PTR cellFaces)
{  
  fd->setStateID(RIGHT, stateID);
  CFreal* statePtrR = (stype > 0) ? &kd->states[stateID*PHYS::NBEQS] : &kd->ghostStates[stateID*PHYS::NBEQS];  
  CFreal* nodePtrR = (stype > 0) ? &kd->centerNodes[stateID*PHYS::DIM] : &kd->ghostNodes[stateID*PHYS::DIM];  
  setState<PHYS>(fd->getState(RIGHT), statePtrR, fd->getNode(RIGHT), nodePtrR);
  
  fd->setIsBFace(stype < 0);
  fd->setStateID(LEFT, cellID);
  const CFuint faceID = cellFaces[f*kd->nbCells + cellID];
  fd->setIsOutward(kd->isOutward[faceID] == cellID);
  
  CFreal* statePtrL = &kd->states[cellID*PHYS::NBEQS];
  CFreal* nodePtrL = &kd->centerNodes[cellID*PHYS::DIM];
  setState<PHYS>(fd->getState(LEFT), statePtrL, fd->getNode(LEFT), nodePtrL);
  setFaceNormal<PHYS>(fd, &kd->normals[faceID*PHYS::DIM]);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, CFuint SIZE>
void print(const std::string& name, T* array) 
{
  std::cout << name << " = \t";
  for (CFuint i = 0; i < SIZE; ++i) {
    std::cout.precision(10); std::cout << array[i] << " ";
  }
  std::cout << "\n";
}
      
//////////////////////////////////////////////////////////////////////////////

template <typename T, CFuint SIZE>
void printArray(T* array) 
{
  for (CFuint i = 0; i < SIZE; ++i) {
    std::cout << array[i] << " ";
  }
  std::cout << "\n";
}

//////////////////////////////////////////////////////////////////////////////

template <typename MODEL>
HOST_DEVICE void computeFaceCentroid(const CellData::Itr* cell, const CFuint faceIdx, 
				     const CFreal* nodes, CFreal* midFaceCoord)
{  
  CudaEnv::CFVecSlice<CFreal, MODEL::DIM> coord(midFaceCoord);
  coord = 0.;
  const CFuint nbFaceNodes = cell->getNbFaceNodes(faceIdx);
  const CFreal ovNbFaceNodes = 1./(static_cast<CFreal>(nbFaceNodes));
  for (CFuint n = 0; n < nbFaceNodes; ++n) {
    const CFuint cellNodeID = cell->getNodeID(faceIdx, n);
    const CFuint nodeID = cell->getNodeID(faceIdx,n);
    const CFreal* faceNode = &nodes[nodeID*MODEL::DIM];
    for (CFuint d = 0; d < MODEL::DIM; ++d) {
      coord[d] += faceNode[d];
    }
  }
  coord *= ovNbFaceNodes;
}

//////////////////////////////////////////////////////////////////////////////

template <typename PHYS, typename POLYREC>
__global__ void computeGradientsKernel(typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor,
				       const CFuint nbCells,
				       CFreal* states, 
				       CFreal* nodes,
				       CFreal* centerNodes,
				       CFreal* ghostStates,
				       CFreal* ghostNodes,
				       CFreal* uX,
				       CFreal* uY,
				       CFreal* uZ,
				       CFreal* limiter,
				       CFreal* updateCoeff, 
				       CFreal* rhs,
				       CFreal* normals,
				       CFint* isOutward,
				       const CFuint* cellInfo,
				       const CFuint* cellStencil,
				       const CFuint* cellFaces,
				       const CFuint* cellNodes,
				       const CFint*  neighborTypes,
				       const Framework::CellConn* cellConn)
{    
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (cellID < nbCells) {    
    KernelData<CFreal> kd (nbCells, states, nodes, centerNodes, ghostStates, ghostNodes, updateCoeff, 
			   rhs, normals, uX, uY, uZ, isOutward);
    
    // compute and store cell gradients at once 
    POLYREC polyRec(dcor);
    CellData cells(nbCells, cellInfo, cellStencil, cellFaces, cellNodes, neighborTypes, cellConn);
    CellData::Itr cell = cells.getItr(cellID);
    polyRec.computeGradients(&states[cellID*PHYS::NBEQS], &centerNodes[cellID*PHYS::DIM], &kd, &cell);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

template <typename PHYS, typename POLYREC, typename LIMITER>
__global__ void computeLimiterKernel(typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE>* dcol,
				     typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor,
				     const CFuint nbCells,
				     CFreal* states, 
				     CFreal* nodes,
				     CFreal* centerNodes,
				     CFreal* ghostStates,
				     CFreal* ghostNodes,
				     CFreal* uX,
				     CFreal* uY,
				     CFreal* uZ,
				     CFreal* limiter,
				     CFreal* updateCoeff, 
				     CFreal* rhs,
				     CFreal* normals,
				     CFint* isOutward,
				     const CFuint* cellInfo,
				     const CFuint* cellStencil,
				     const CFuint* cellFaces,
				     const CFuint* cellNodes,
				     const CFint*  neighborTypes,
				     const Framework::CellConn* cellConn)
{    
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
 
  if (cellID < nbCells) {    
    // compute all cell quadrature points at once (size of this array is overestimated)
    CFreal midFaceCoord[PHYS::DIM*PHYS::DIM*2];
    
    CellData cells(nbCells, cellInfo, cellStencil, cellFaces, cellNodes, neighborTypes, cellConn);
    CellData::Itr cell = cells.getItr(cellID);
    const CFuint nbFacesInCell = cell.getNbFacesInCell();
    for (CFuint f = 0; f < nbFacesInCell; ++f) { 
      computeFaceCentroid<PHYS>(&cell, f, nodes, &midFaceCoord[f*PHYS::DIM]);
    }
    
    // compute cell-based limiter at once
    KernelData<CFreal> kd (nbCells, states, nodes, centerNodes, ghostStates, ghostNodes, updateCoeff, 
			   rhs, normals, uX, uY, uZ, isOutward);
    LIMITER limt(dcol);
    
    if (dcor->currRes > dcor->limitRes && (dcor->limitIter > 0 && dcor->currIter < dcor->limitIter)) {	
      limt.limit(&kd, &cell, &midFaceCoord[0], &limiter[cellID*PHYS::NBEQS]);
    }
    else {
      if (!dcor->freezeLimiter) {
	// historical modification of the limiter
	CudaEnv::CFVec<CFreal,PHYS::NBEQS> tmpLimiter;
	limt.limit(&kd, &cell, &midFaceCoord[0], &tmpLimiter[0]);
	CFuint currID = cellID*PHYS::NBEQS;
	for (CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar, ++currID) {
	  limiter[currID] = min(tmpLimiter[iVar],limiter[currID]);
	}
      }
    }
  }
}
  
//////////////////////////////////////////////////////////////////////////////

template <typename SCHEME, typename POLYREC, typename LIMITER>
__global__ void computeFluxJacobianKernel(typename SCHEME::BASE::template DeviceConfigOptions<NOTYPE>* dcof,
					  typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor,
					  typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE>* dcol,
					  typename NumericalJacobian::DeviceConfigOptions<typename SCHEME::MODEL>* dcon,
					  typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE>* dcop,
					  const CFuint nbCells,
					  const CFuint startCellID,
					  CFreal* states, 
					  CFreal* nodes,
					  CFreal* centerNodes,
					  CFreal* ghostStates,
					  CFreal* ghostNodes,
					  CFreal* blockJacob,
					  CFuint* blockStart,
					  CFreal* uX,
					  CFreal* uY,
					  CFreal* uZ,
					  CFreal* limiter,
					  CFreal* updateCoeff, 
					  CFreal* rhs,
					  CFreal* normals,
					  CFint* isOutward,
					  const CFuint* cellInfo,
					  const CFuint* cellStencil,
					  const CFuint* cellFaces,
					  const CFuint* cellNodes,
					  const CFint* neighborTypes,
					  const Framework::CellConn* cellConn)
{  
  typedef typename SCHEME::MODEL PHYS;
  
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x + startCellID;
  
  if (cellID < nbCells) {
    // compute and store cell gradients at once 
    POLYREC polyRec(dcor);
    SCHEME  fluxScheme(dcof);
    LIMITER limt(dcol);
    NumericalJacobian::DeviceFunc<PHYS> numJacob(dcon);
    KernelData<CFreal> kd (nbCells, states, nodes, centerNodes, ghostStates, 
			   ghostNodes, updateCoeff, rhs, normals, uX, uY, uZ, isOutward);
    
    // compute all cell quadrature points at once (array size can be overestimated in 3D)
    const CFuint MAX_NB_FACES = PHYS::DIM*2;
    CFreal midFaceCoord[PHYS::DIM*MAX_NB_FACES];
    CudaEnv::CFVec<CFreal,PHYS::NBEQS> fluxDiff;
    CudaEnv::CFVec<CFreal,PHYS::NBEQS> resBkp;
    FluxData<PHYS> currFd; currFd.initialize();
    typename SCHEME::MODEL pmodel(dcop);
    
    // reset the rhs and update coefficients to 0
    CudaEnv::CFVecSlice<CFreal,PHYS::NBEQS> res(&rhs[cellID*PHYS::NBEQS]);
    res = 0.;
    updateCoeff[cellID] = 0.;
    
    CellData cells(nbCells, cellInfo, cellStencil, cellFaces, cellNodes, neighborTypes, cellConn);
    CellData::Itr cell = cells.getItr(cellID);
    const CFuint nbFacesInCell = cell.getNbActiveFacesInCell();
    const CFuint nbRows = nbFacesInCell + 1;
    const CFuint bStartCellID = blockStart[cellID];
    
    // this block accumulator represents a column block (nbFaces+1 x 1)
    BlockAccumulatorBaseCUDA acc(nbRows, 1, PHYS::NBEQS, &blockJacob[bStartCellID]);
    acc.reset();
    
    // compute the face flux and flux numerical jacobian within the same loop
    for (CFuint f = 0; f < nbFacesInCell; ++f) { 
      const CFint stype = cell.getNeighborType(f);
      
      if (stype != 0) { // skip all partition faces
	const CFuint stateID = cell.getNeighborID(f);
	setFluxData(f, stype, stateID, cellID, &kd, &currFd, cellFaces);
	
	// compute face quadrature points (face centroids)
	CFreal* faceCenters = &midFaceCoord[f*PHYS::DIM];
	computeFaceCentroid<PHYS>(&cell, f, nodes, faceCenters);
	
	// extrapolate solution on quadrature points on both sides of the face
	polyRec.extrapolateOnFace(&currFd, faceCenters, uX, uY, uZ, limiter);

	fluxScheme(&currFd, &pmodel);

        for (CFuint iEq = 0; iEq < PHYS::NBEQS; ++iEq) {
	  CFreal value = currFd.getResidual()[iEq];
	  res[iEq]   -= value;  // update the residual 
	  resBkp[iEq] = value;  // backup the current face-based residual
	}


	// update the update coefficient
	updateCoeff[cellID] += currFd.getUpdateCoeff();
		
	// only contribution from internal faces is computed here  
	if (stype > 0) { 	  
	  currFd.setIsPerturb(true);
	  // flux jacobian computation
	  for (CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar) {
	    // here we perturb the current variable for the left cell state
	    numJacob.perturb(iVar, &currFd.getState(LEFT)[iVar]);
	    
	    // extrapolate solution on quadrature points on both sides of the face
	    const CFreal rstateBkpL = currFd.getRstate(LEFT)[iVar];
	    polyRec.extrapolateOnFace(iVar, &currFd, faceCenters, uX, uY, uZ, limiter);
	    fluxScheme(&currFd, &pmodel); // compute the convective flux across the face
	    
	    // compute the numerical jacobian of the flux
	    CudaEnv::CFVecSlice<CFreal,PHYS::NBEQS> resPert(currFd.getResidual());
	    numJacob.computeDerivative(&resBkp, &resPert, &fluxDiff);
	    
	    // contribution to the row corresponding of the current cell
	    // this subblock gets all contributions from all face cells
	    acc.addValues(0, 0, iVar, &fluxDiff[0]);
	    
	    // contribution to row corresponding to the f+1 cell: 
	    // this is the flux jacobian contribution for the neighbor cells
	    // due to the currently perturbed cell state and is opposite in sign
	    // because the outward normal for neighbors is inward for the current cell
	    fluxDiff *= -1.0;
	    acc.addValues(f+1, 0, iVar, &fluxDiff[0]); 
	    
	    // restore perturbed states
	    currFd.getRstate(LEFT)[iVar] = rstateBkpL;
	    numJacob.restore(&currFd.getState(LEFT)[iVar]);
	  }
	  
	  currFd.setIsPerturb(false);
	}
      }
    }
    
    
  }
}

//////////////////////////////////////////////////////////////////////////////
  
template <typename SCHEME, typename POLYREC, typename LIMITER>
void computeFluxJacobianCPU(typename SCHEME::BASE::template DeviceConfigOptions<NOTYPE>* dcof,
			    typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor,
			    typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE>* dcol,
			    typename NumericalJacobian::DeviceConfigOptions<typename SCHEME::MODEL>* dcon,
			    typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE>* dcop,
			    const CFuint nbCells,
			    CFreal* states, 
			    CFreal* nodes,
			    CFreal* centerNodes,
			    CFreal* ghostStates,
			    CFreal* ghostNodes, 
			    CFreal* blockJacob,
			    CFuint* blockStart,
			    CFreal* uX,
			    CFreal* uY,
			    CFreal* uZ,
			    CFreal* limiter,
			    CFreal* updateCoeff, 
			    CFreal* rhs,
			    CFreal* normals,
			    CFint* isOutward,
			    const CFuint* cellInfo,
			    const CFuint* cellStencil,
			    const CFuint* cellFaces,
			    const CFuint* cellNodes,
			    const CFint* neighborTypes,
			    const Framework::CellConn* cellConn)
{ 
  using namespace std;
  
  typedef typename SCHEME::MODEL PHYS;
  
  CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
  timer.start();
  
  FluxData<PHYS> fd; fd.initialize();
  FluxData<PHYS>* currFd = &fd;
  cf_assert(currFd != CFNULL);
  SCHEME fluxScheme(dcof);
  POLYREC polyRec(dcor);
  LIMITER limt(dcol);
  
  CellData cells(nbCells, cellInfo, cellStencil, cellFaces, cellNodes, neighborTypes, cellConn);
  KernelData<CFreal> kd(nbCells, states, nodes, centerNodes, ghostStates, ghostNodes, updateCoeff, 
			rhs, normals, uX, uY, uZ, isOutward);
  
  const CFuint MAX_NB_FACES = PHYS::DIM*2;
  CFreal midFaceCoord[PHYS::DIM*MAX_NB_FACES];
  CudaEnv::CFVec<CFreal,PHYS::NBEQS> fluxDiff;
  CudaEnv::CFVec<CFreal,PHYS::NBEQS> resBkp;
  CudaEnv::CFVec<CFreal,PHYS::NBEQS> tmpLimiter;
  NumericalJacobian::DeviceFunc<PHYS> numJacob(dcon);
  PHYS pmodel(dcop);
  
  // compute the cell-based gradients
  for (CFuint cellID = 0; cellID < nbCells; ++cellID) {
    CellData::Itr cell = cells.getItr(cellID);
    polyRec.computeGradients(&states[cellID*PHYS::NBEQS], &centerNodes[cellID*PHYS::DIM], &kd, &cell);
  }
  
  // printGradients<PHYS::NBEQS>(uX, uY, uZ, nbCells);
  CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::computeFluxJacobianCPU() => computing gradients took " << timer.elapsed() << " s\n");
  timer.start();
  
  // compute the cell based limiter
  for (CellData::Itr cell = cells.begin(); cell <= cells.end(); ++cell) {
    // compute all cell quadrature points at once (size of this array is overestimated)
    const CFuint nbFacesInCell = cell.getNbFacesInCell();
    for (CFuint f = 0; f < nbFacesInCell; ++f) { 
      computeFaceCentroid<PHYS>(&cell, f, nodes, &midFaceCoord[f*PHYS::DIM]);
    }
    
    const CFuint cellID = cell.getCellID();
    if (dcor->currRes > dcor->limitRes && (dcor->limitIter > 0 && dcor->currIter < dcor->limitIter)) {	
      // compute cell-based limiter
      limt.limit(&kd, &cell, &midFaceCoord[0], &limiter[cellID*PHYS::NBEQS]);
    }
    else {
      if (!dcor->freezeLimiter) {
	// historical modification of the limiter
	limt.limit(&kd, &cell, &midFaceCoord[0], &tmpLimiter[0]);
	CFuint currID = cellID*PHYS::NBEQS;
	for (CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar, ++currID) {
	  limiter[currID] = min(tmpLimiter[iVar],limiter[currID]);
	}
      }
    }
  }
  
  // printLimiter<PHYS::NBEQS>(limiter, nbCells);
  
  CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::computeFluxJacobianCPU() => computing limiter took " << timer.elapsed() << " s\n");
  timer.start();
  
  // compute the fluxes and the jacobian
  for (CellData::Itr cell = cells.begin(); cell <= cells.end(); ++cell) {
    // reset the rhs and update coefficients to 0
    const CFuint cellID = cell.getCellID();
    CudaEnv::CFVecSlice<CFreal,PHYS::NBEQS> res(&rhs[cellID*PHYS::NBEQS]);
    res = 0.;
    updateCoeff[cellID] = 0.;
    
    const CFuint nbFacesInCell = cell.getNbActiveFacesInCell();
    const CFuint nbRows = nbFacesInCell + 1;
    const CFuint bStartCellID = blockStart[cellID];
        
    // this block accumulator represents a column block (nbFaces+1 x 1)
    BlockAccumulatorBaseCUDA acc(nbRows, 1, PHYS::NBEQS, &blockJacob[bStartCellID]);
    acc.reset();
    
    for (CFuint f = 0; f < nbFacesInCell; ++f) { 
      const CFint stype = cell.getNeighborType(f);
      if (stype != 0) { // skip all partition faces
	const CFuint stateID =  cell.getNeighborID(f);
	setFluxData(f, stype, stateID, cellID, &kd, currFd, cellFaces);
	
	// compute face quadrature points (centroid)
	CFreal* faceCenters = &midFaceCoord[f*PHYS::DIM];
	computeFaceCentroid<PHYS>(&cell, f, nodes, faceCenters);
	
	// extrapolate solution on quadrature points on both sides of the face
	polyRec.extrapolateOnFace(currFd, faceCenters, uX, uY, uZ, limiter);
		
	/*if (cellID == 0) {    
	  FluxData<PHYS>* data = currFd; 
	  PHYS* model = &pmodel;
	  CudaEnv::CFVec<CFreal,PHYS::NBEQS> m_tmp;
	  CudaEnv::CFVec<CFreal,PHYS::NBEQS> m_tmp2;
	  CudaEnv::CFVec<CFreal,PHYS::NBEQS> m_sumFlux;
	  CudaEnv::CFVec<CFreal,PHYS::DATASIZE> m_pdata;
	  CudaEnv::CFVec<CFreal,PHYS::DIM> m_tempUnitNormal;
	  CudaEnv::CFVecSlice<CFreal,PHYS::DIM> unitNormal(data->getUnitNormal());
	  const CFreal coeff = (data->isOutward()) ? 1. : -1.;
	  m_tempUnitNormal = coeff*unitNormal;
	  CudaEnv::CFVecSlice<CFreal,PHYS::NBEQS> flux(data->getResidual());
	  typename PHYS::UPDATE_VS* updateVS = model->getUpdateVS();
	  
	  updateVS->computePhysicalData(data->getRstate(1), data->getNode(1), &m_pdata[0]);
	  updateVS->getFlux(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
	  flux = 0.5*m_tmp;
	  
	  updateVS->computeEigenValues(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
	  
	  std::cout << std::endl << "stateIDs [" << cellID << "," << stateID << "]" << std::endl;
	  // std::cout.precision(12); std::cout << "centerNodes[0]    = "; printArray<CFreal,PHYS::DIM>(&centerNodes[cellID*PHYS::DIM]);
	  
	  std::cout.precision(12);std::cout << "dataR   = "; printArray<CFreal,PHYS::DATASIZE>(&m_pdata[0]);
	  // std::cout.precision(12);std::cout << "eigenR  = "; printArray<CFreal,PHYS::NBEQS>(&m_tmp[0]);
	  std::cout.precision(12);std::cout << "normalR = "; printArray<CFreal,PHYS::DIM>(&m_tempUnitNormal[0]);
	  
	  CFreal aR = 0.0;
	  for (CFuint i = 0; i < PHYS::NBEQS; ++i) {
	    aR = max(aR, abs(m_tmp[i]));
	  }
	  
	  // left physical data, flux and eigenvalues
	  updateVS->computePhysicalData(data->getRstate(0), data->getNode(0), &m_pdata[0]);
	  updateVS->getFlux(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
	  flux += 0.5*m_tmp;
	  
	  updateVS->computeEigenValues(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
	  std::cout.precision(12);std::cout << "dataL   = "; printArray<CFreal,PHYS::DATASIZE>(&m_pdata[0]);
	  // std::cout.precision(12);std::cout << "eigenL  = "; printArray<CFreal,PHYS::NBEQS>(&m_tmp[0]);
	  std::cout.precision(12);std::cout << "normalL = "; printArray<CFreal,PHYS::DIM>(&m_tempUnitNormal[0]);
	  
	  // compute update coefficient
	  if (!data->isPerturb()) {    
	    const CFreal k = max(m_tmp.max(), 0.)*data->getFaceArea();
	    std::cout << "lambda = " << max(m_tmp.max(), 0.) << ", area =" << data->getFaceArea() << std::endl;
	    std::cout << "updateCoeff [" <<cellID << "] = " << k << std::endl;
	    // data->setUpdateCoeff(k);
	    //abort();
	  }
	  
	  CFreal aL = 0.0;
	  for (CFuint i = 0; i < PHYS::NBEQS; ++i) {
	    aL = max(aL, abs(m_tmp[i]));
	  }
	  
	  const CFreal a = (aR > aL) ? aR : aL; // max(aR,aL);
	  const CFreal aDiff = a*dcof->currentDiffRedCoeff;
	  
	  Framework::VarSetTransformerT<typename PHYS::UPDATE_VS, typename PHYS::SOLUTION_VS, NOTYPE>* up2Sol = 
	    model->getUpdateToSolution();
	  
	  // transform to solution variables
	  up2Sol->transform(data->getRstate(LEFT), &m_tmp[0]);
	  up2Sol->transform(data->getRstate(RIGHT), &m_tmp2[0]);
	  
	  CudaEnv::CFVecSlice<CFreal,PHYS::NBEQS> stateL(&m_tmp[0]);
	  CudaEnv::CFVecSlice<CFreal,PHYS::NBEQS> stateR(&m_tmp2[0]);
	  
	  m_sumFlux = 2.*flux;
	  
	  std::cout.precision(12);;std::cout << "UR       = ";printArray<CFreal,PHYS::NBEQS>(&stateR[0]);
	  std::cout.precision(12);;std::cout << "UL       = ";printArray<CFreal,PHYS::NBEQS>(&stateL[0]);
	  std::cout.precision(12);;std::cout << "_sumFlux = ";printArray<CFreal,PHYS::NBEQS>(&m_sumFlux[0]);
	  std::cout.precision(12);;std::cout << "aDiff    = " << aDiff << std::endl;
	  
	  flux -= (0.5*aDiff)*(stateR - stateL);
	  
	  
	  // // NOTE THE AREA HERE !!!!!!!!!!!!!!!!
	  // flux *= data->getFaceArea();
	  }*/

	
	////////
	

       
	fluxScheme(currFd, &pmodel); // compute the convective flux across the face
        

	for (CFuint iEq = 0; iEq < PHYS::NBEQS; ++iEq) {
	  const CFreal value = currFd->getResidual()[iEq];
	  res[iEq]   -= value;  // update the residual 
	  resBkp[iEq] = value;  // backup the current face-based residual
	}
	
	// update the update coefficient
	updateCoeff[cellID] += currFd->getUpdateCoeff();
	
	// only contribution from internal faces is computed here  
	if (stype > 0) { 
	  currFd->setIsPerturb(true);
	  
	  // flux jacobian computation
	  for (CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar) {
	    // here we perturb the current variable for the left cell state
	    numJacob.perturb(iVar, &currFd->getState(LEFT)[iVar]);
	    
	    // extrapolate solution on quadrature points on both sides of the face
	    const CFreal rstateBkpL = currFd->getRstate(LEFT)[iVar];
	    polyRec.extrapolateOnFace(iVar, currFd, faceCenters, uX, uY, uZ, limiter);
	    fluxScheme(currFd, &pmodel); // compute the convective flux across the face
	    
	    // compute the numerical jacobian of the flux
	    CudaEnv::CFVecSlice<CFreal,PHYS::NBEQS> resPert(currFd->getResidual());
	    numJacob.computeDerivative(&resBkp, &resPert, &fluxDiff);
	    
	    
	    
	    ///
	    /*if (cellID == 0 && stateID == 1 && iVar == 0) {
	      cout << "\n left 0, right 1 \n";
	      cout << iVar << " => resBkp   = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&resBkp[0]);
	      cout << iVar << " => resPert  = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&resPert[0]);
	      cout << iVar << " => fluxDiff = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&fluxDiff[0]);
	    }
	    
	    if (stateID == 0 && cellID == 1 && iVar == 0) {
	      cout << "\n left 1, right 0 \n";
	      cout << iVar << " => resBkp   = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&resBkp[0]);
	      cout << iVar << " => resPert  = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&resPert[0]);
	      cout << iVar << " => fluxDiff = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&fluxDiff[0]);
	      }*/
	    ///
	    
	    
	    
	    // flux is computed with the outward normal, so the sign is correct here
	    // contribution to the row corresponding of the current cell
	    // this subblock gets all contributions from all face cells
	    acc.addValues(0, 0, iVar, &fluxDiff[0]);
	    
	    // contribution to row corresponding to the f+1 cell: 
	    // this is the flux jacobian contribution for the neighbor cells
	    // due to the currently perturbed cell state and is opposite in sign
	    // because the outward normal for neighbors is inward for the current cell
	    fluxDiff *= -1.0;
	    
	    
	    ///
	    /* if (cellID == 0 && stateID == 1 && iVar == 0) {
	      cout << "\n left 0, right 1 \n";
	      cout << iVar << " => resBkp   = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&resBkp[0]);
	      cout << iVar << " => resPert  = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&resPert[0]);
	      cout << iVar << " => fluxDiff = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&fluxDiff[0]);
	    }
	    
	    if (stateID == 0 && cellID == 1 && iVar == 0) {
	      cout << "\n left 1, right 0 \n";
	      cout << iVar << " => resBkp   = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&resBkp[0]);
	      cout << iVar << " => resPert  = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&resPert[0]);
	      cout << iVar << " => fluxDiff = "; std::cout.precision(12);printArray<CFreal,PHYS::NBEQS>(&fluxDiff[0]);
	      }*/
	    ///
	    
	    
	    acc.addValues(f+1, 0, iVar, &fluxDiff[0]);   
	    
	    // restore perturbed states
	    currFd->getRstate(LEFT)[iVar] = rstateBkpL;
	    numJacob.restore(&currFd->getState(LEFT)[iVar]);
	  }
	  	  
	  currFd->setIsPerturb(false);
	}

      }
    }
    //if (abs(res[6]) <= 1e-3){res[6] = 0.0;} 
  } 
  
  CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::computeFluxJacobianCPU() took " << timer.elapsed() << " s\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename SCHEME, typename PHYSICS, typename POLYREC, typename LIMITER, CFuint NB_BLOCK_THREADS>
void FVMCC_ComputeRhsJacobCell<SCHEME,PHYSICS,POLYREC,LIMITER,NB_BLOCK_THREADS>::execute()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  CFTRACEBEGIN;
  
  CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCell::execute() START\n");
  
  initializeComputationRHS();
  
  const CFuint nbCells = this->socket_states.getDataHandle().size();
  cf_assert(nbCells > 0);
  DataHandle<CFreal> updateCoeff = this->socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = this->socket_rhs.getDataHandle();
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();  
  
  SafePtr<SCHEME> lf = this->getMethodData().getFluxSplitter().template d_castTo<SCHEME>();
  SafePtr<POLYREC> pr = this->getMethodData().getPolyReconstructor().template d_castTo<POLYREC>();
  SafePtr<LIMITER> lm = this->getMethodData().getLimiter().template d_castTo<LIMITER>();
  SafePtr<typename PHYSICS::PTERM> phys = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<typename PHYSICS::PTERM>();
  
  typedef typename SCHEME::template  DeviceFunc<GPU, PHYSICS> FluxScheme;  
  typedef typename POLYREC::template DeviceFunc<PHYSICS> PolyRec;  
  typedef typename LIMITER::template DeviceFunc<PHYSICS> Limiter;  
  
  CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
  
  if (this->m_onGPU) {
    
    timer.start();
    // copy of data that change at every iteration
    this->socket_states.getDataHandle().getGlobalArray()->put(); 
    this->m_ghostStates.put();
    
    CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() => CPU-->GPU data transfer took " << timer.elapsed() << " s\n");
    timer.start();
    
    ConfigOptionPtr<POLYREC, NOTYPE, GPU> dcor(pr);
    ConfigOptionPtr<LIMITER, NOTYPE, GPU> dcol(lm);
    ConfigOptionPtr<SCHEME,  NOTYPE, GPU> dcof(lf);
    ConfigOptionPtr<typename PHYSICS::PTERM, NOTYPE, GPU> dcop(phys);
    
    const CFuint blocksPerGrid = CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
    const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
    
    CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() <blocksPerGrid, nThreads> = <" 
	  <<  blocksPerGrid << "," << nThreads << ">\n");
    
    //dim3 blocks(this->m_nbBlocksPerGridX, this->m_nbBlocksPerGridY);
    
    timer.start();
    
    // compute the cell-based gradients
    computeGradientsKernel<PHYSICS, PolyRec> <<<blocksPerGrid,nThreads>>> 
      (dcor.getPtr(),
       nbCells,
       this->socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
       this->socket_nodes.getDataHandle().getGlobalArray()->ptrDev(),
       this->m_centerNodes.ptrDev(), 
       this->m_ghostStates.ptrDev(),
       this->m_ghostNodes.ptrDev(),
       this->socket_uX.getDataHandle().getLocalArray()->ptrDev(),
       this->socket_uY.getDataHandle().getLocalArray()->ptrDev(),
       this->socket_uZ.getDataHandle().getLocalArray()->ptrDev(),
       this->socket_limiter.getDataHandle().getLocalArray()->ptrDev(),
       updateCoeff.getLocalArray()->ptrDev(), 
       rhs.getLocalArray()->ptrDev(),
       normals.getLocalArray()->ptrDev(),
       isOutward.getLocalArray()->ptrDev(),
       this->m_cellInfo.ptrDev(),
       this->m_cellStencil.ptrDev(),
       this->m_cellFaces->getPtr()->ptrDev(),
       this->m_cellNodes->getPtr()->ptrDev(),
       this->m_neighborTypes.ptrDev(),
       this->m_cellConn.ptrDev());
    
    CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() => computeGradientsKernel took " << timer.elapsed() << " s\n");
    
    timer.start();
    
    // compute the limiter in each cell
    computeLimiterKernel<PHYSICS, PolyRec, Limiter> <<<blocksPerGrid,nThreads>>> 
      (dcol.getPtr(), 
       dcor.getPtr(), 
       nbCells,
       this->socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
       this->socket_nodes.getDataHandle().getGlobalArray()->ptrDev(),
       this->m_centerNodes.ptrDev(), 
       this->m_ghostStates.ptrDev(),
       this->m_ghostNodes.ptrDev(),
       this->socket_uX.getDataHandle().getLocalArray()->ptrDev(),
       this->socket_uY.getDataHandle().getLocalArray()->ptrDev(),
       this->socket_uZ.getDataHandle().getLocalArray()->ptrDev(),
       this->socket_limiter.getDataHandle().getLocalArray()->ptrDev(),
       updateCoeff.getLocalArray()->ptrDev(), 
       rhs.getLocalArray()->ptrDev(),
       normals.getLocalArray()->ptrDev(),
       isOutward.getLocalArray()->ptrDev(),
       this->m_cellInfo.ptrDev(),
       this->m_cellStencil.ptrDev(),
       this->m_cellFaces->getPtr()->ptrDev(),
       this->m_cellNodes->getPtr()->ptrDev(),
       this->m_neighborTypes.ptrDev(),
       this->m_cellConn.ptrDev());
    
    CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() => computeLimiterKernel took " << timer.elapsed() << " s\n");
    
    timer.start();
    // compute the flux jacobian in each cell
    CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCell::execute() => Configuring method \n");
    ConfigOptionPtr<NumericalJacobian, PHYSICS, GPU> dcon
      (&this->getMethodData().getNumericalJacobian());
    CFuint startCellID = 0;
    CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCell::execute() => End of Configuring method \n");
    
    CFreal FluxTime = 0.0;
    CFreal UpdateSystemTime = 0.0;
    
    for (CFuint s = 0; s < m_nbCellsInKernel.size(); ++s) {
      CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCell::execute() => loop " << s << " of " << m_nbCellsInKernel.size() << "\n");
      computeFluxJacobianKernel<FluxScheme, PolyRec, Limiter> <<<m_nbKernelBlocks,nThreads>>> 
	(dcof.getPtr(),
	 dcor.getPtr(),
	 dcol.getPtr(),
	 dcon.getPtr(),
	 dcop.getPtr(),
	 nbCells,
	 startCellID,
	 this->socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
	 this->socket_nodes.getDataHandle().getGlobalArray()->ptrDev(),
	 this->m_centerNodes.ptrDev(), 
	 this->m_ghostStates.ptrDev(),
	 this->m_ghostNodes.ptrDev(),
	 m_blockJacobians.ptrDev(), 
	 m_blockStart.ptrDev(),
	 this->socket_uX.getDataHandle().getLocalArray()->ptrDev(),
	 this->socket_uY.getDataHandle().getLocalArray()->ptrDev(),
	 this->socket_uZ.getDataHandle().getLocalArray()->ptrDev(),
	 this->socket_limiter.getDataHandle().getLocalArray()->ptrDev(),
	 updateCoeff.getLocalArray()->ptrDev(), 
	 rhs.getLocalArray()->ptrDev(),
	 normals.getLocalArray()->ptrDev(),
	 isOutward.getLocalArray()->ptrDev(),
	 this->m_cellInfo.ptrDev(),
	 this->m_cellStencil.ptrDev(),
	 this->m_cellFaces->getPtr()->ptrDev(),
	 this->m_cellNodes->getPtr()->ptrDev(),
	 this->m_neighborTypes.ptrDev(),
	 this->m_cellConn.ptrDev());

      FluxTime += timer.elapsed();
      timer.start();
      CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCell::execute() => m_blockJacobians.get() \n");
      m_blockJacobians.get();
      // update the portion of system matrix computed by this kernel
      CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCell::execute() => updateSystemMatrix(" << s <<") \n");
      updateSystemMatrix(s);
      startCellID += m_nbCellsInKernel[s];
      UpdateSystemTime += timer.elapsed();
      timer.start();
    }
    
    CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() => computeFluxJacobianKernel took " << FluxTime << " s\n");
    CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() => updateSystemMatrix took took " << UpdateSystemTime << " s\n");
    
    rhs.getLocalArray()->get();
    updateCoeff.getLocalArray()->get();

    CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() => GPU-->CPU data transfer took " << timer.elapsed() << " s\n");
  }
  else {
    ConfigOptionPtr<SCHEME>  dcof(lf);
    ConfigOptionPtr<POLYREC> dcor(pr);
    ConfigOptionPtr<LIMITER> dcol(lm);
    ConfigOptionPtr<NumericalJacobian, PHYSICS> dcon(&this->getMethodData().getNumericalJacobian());
    ConfigOptionPtr<typename PHYSICS::PTERM> dcop(phys);

    computeFluxJacobianCPU<FluxScheme, PolyRec, Limiter>
      (dcof.getPtr(),
       dcor.getPtr(),
       dcol.getPtr(),
       dcon.getPtr(),
       dcop.getPtr(),
       nbCells,
       this->socket_states.getDataHandle().getGlobalArray()->ptr(), 
       this->socket_nodes.getDataHandle().getGlobalArray()->ptr(),
       this->m_centerNodes.ptr(), 
       this->m_ghostStates.ptr(),
       this->m_ghostNodes.ptr(),
       m_blockJacobians.ptr(), 
       m_blockStart.ptr(),
       this->socket_uX.getDataHandle().getLocalArray()->ptr(),
       this->socket_uY.getDataHandle().getLocalArray()->ptr(),
       this->socket_uZ.getDataHandle().getLocalArray()->ptr(),
       this->socket_limiter.getDataHandle().getLocalArray()->ptr(),
       updateCoeff.getLocalArray()->ptr(), 
       rhs.getLocalArray()->ptr(),
       normals.getLocalArray()->ptr(),
       isOutward.getLocalArray()->ptr(),
       this->m_cellInfo.ptr(),
       this->m_cellStencil.ptr(),
       this->m_cellFaces->getPtr()->ptr(),
       this->m_cellNodes->getPtr()->ptr(),
       this->m_neighborTypes.ptr(),
       this->m_cellConn.ptr());
    
    // update the system matrix
    timer.start();
    updateSystemMatrix(0);
    CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() => updateSystemMatrix took " << timer.elapsed() << " s\n");
  }
  
  timer.start();
  // compute flux jacobians on boundaries
  executeBC();
  CFLog(NOTICE, "FVMCC_ComputeRhsJacobCell::execute() => executeBC() took " << timer.elapsed() << " s\n");

  finalizeComputationRHS();
  
  CFLog(VERBOSE, "FVMCC_ComputeRhsJacobCell::execute() END\n");
  
/*
     const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
     //DataHandle<CFreal> rhs = this->socket_rhs.getDataHandle();
     DataHandle<State*, GLOBAL> states = this->socket_states.getDataHandle();
     //DataHandle<CFreal> updateCoeff = this->socket_updateCoeff.getDataHandle();
     CFreal* state = this->socket_states.getDataHandle().getGlobalArray()->ptr();
   for (CFuint iState = 0; iState < states.size(); ++iState) {
   // if (iState == 2299){
     CFLog(VERBOSE, " \t iState: " << iState << "\t UpdateCoeff: " << updateCoeff[iState] << "\n");
     for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
       //cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << rhs(iState, iEq, nbEqs);
       CFLog(VERBOSE, rhs(iState, iEq, nbEqs) << " \n \t" << state[iState*nbEqs + iEq] << " \n");
       //CFLog(VERBOSE, rhs(iState, iEq, nbEqs) << "\n");
     }
     cout << endl;
   // }
   }
*/



 /* for (int i = 0; i < updateCoeff.size(); ++i) {
    // for (int i = 0; i < 10000; ++i) {
    std::cout << "updateCoeff[" << i << "] = " << updateCoeff[i]  << std::endl;
    std::cout << "rhs[" << i << "] = ";
    for (int j = 0; j < 9; ++j) {
      std::cout << rhs[i*9+j] << " ";
    }
    std::cout << std::endl;
  }*/
  
  // abort();
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

template <CFuint NBEQS>
void printGradients(CFreal* uX, CFreal* uY, CFreal* uZ, CFuint nbCells)
{  
  CFuint idxr = 0;
  for (CFuint cellID = 0; cellID < nbCells; ++cellID) {
    for (CFuint i = 0; i < NBEQS; ++i, ++idxr) {
      std::cout << "cellID["<< cellID << "], "<< i << " => UX (";
      std::cout.precision(12); std::cout << uX[idxr] << ", " << uY[idxr] << ", " << uZ[idxr] << ")\n";
    }
  } 
}

//////////////////////////////////////////////////////////////////////////////

template <CFuint NBEQS>
void printLimiter(CFreal* limiter, CFuint nbCells)
{ 
  CFuint idxl = 0;
  for (CFuint cellID = 0; cellID < nbCells; ++cellID) {
    std::cout << "cellID["<< cellID << "] => LIM (";
    for (CFuint i = 0; i < NBEQS; ++i, ++idxl) {
      std::cout.precision(12); std::cout << limiter[idxl] << " ";
    }
    std::cout << ")\n";
  }
}
 
//////////////////////////////////////////////////////////////////////////////

   } // namespace FiniteVolume
    
  } // namespace Numerics

} // namespace COOLFluiD
