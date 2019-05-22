#include "FluxReconstructionCUDA/ConvRHSFluxReconstructionCUDA.hh"
#include "Framework/MeshData.hh"
#include "Framework/CellConn.hh"
#include "Config/ConfigOptionPtr.hh"
#include "Framework/CudaDeviceManager.hh"
#include "Common/CUDA/CFVec.hh"
#include "Framework/CudaTimer.hh"

//#include "FiniteVolume/FluxData.hh"
#include "FluxReconstructionMethod/KernelData.hh"
#include "FluxReconstructionMethod/CellData.hh"

#include "FluxReconstructionCUDA/FluxReconstructionCUDA.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/VarSetListT.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "NavierStokes/Euler3DCons.hh"

#include "FluxReconstructionMethod/LaxFriedrichsFlux.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

#define FR_EULER_RHS_PROV(__dim__,__svars__,__uvars__,__nbBThreads__,__providerName__) \
MethodCommandProvider<ConvRHSFluxReconstructionCUDA<LaxFriedrichsFlux, \
                      VarSetListT<Euler##__dim__##__svars__, Euler##__dim__##__uvars__>, __nbBThreads__>, \
		      FluxReconstructionSolverData,FluxReconstructionCUDAModule>	\
FR_RhsEuler##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);
// 48 block threads (default)
FR_EULER_RHS_PROV(2D, Cons, Cons, 48, "EulerFRLaxFried2DCons")
FR_EULER_RHS_PROV(3D, Cons, Cons, 48, "EulerFRLaxFried3DCons")
//FR_NS_RHS_PROV(2D, ProjectionCons, ProjectionPrim, 48, "CellLaxFriedMHD2DPrim")
//FR_NS_RHS_PROV(3D, ProjectionCons, ProjectionPrim, 48, "CellLaxFriedMHD3DPrim")
#undef FR_EULER_RHS_PROV

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
      
//template <typename PHYS>
//HOST_DEVICE inline void setFaceNormal(FluxData<PHYS>* fd, CFreal* normal)
//{
//  CudaEnv::CFVecSlice<CFreal,PHYS::DIM> n(normal);
//  const CFreal area = n.norm2();
//  fd->setFaceArea(area);
//  const CFreal ovArea = 1./area;
//  CudaEnv::CFVecSlice<CFreal,PHYS::DIM> un(fd->getUnitNormal());
//  for (CFuint i = 0; i < PHYS::DIM; ++i) {
//    un[i] = n[i]*ovArea;
//  }
//}
      
//////////////////////////////////////////////////////////////////////////////

//template <typename PHYS, typename PTR>
//HOST_DEVICE void setFluxData(const CFuint f, const CFint stype, 
//			     const CFuint stateID, const CFuint cellID, 
//			     KernelData<CFreal>* kd, FluxData<PHYS>* fd,
//			     PTR cellFaces)
//{
//  fd->setStateID(RIGHT, stateID);
//  CFreal* statePtrR = (stype > 0) ? &kd->states[stateID*PHYS::NBEQS] : &kd->ghostStates[stateID*PHYS::NBEQS];  
//  CFreal* nodePtrR = (stype > 0) ? &kd->centerNodes[stateID*PHYS::DIM] : &kd->ghostNodes[stateID*PHYS::DIM];  
//  setState<PHYS>(fd->getState(RIGHT), statePtrR, fd->getNode(RIGHT), nodePtrR);
//  
//  fd->setIsBFace(stype < 0);
//  fd->setStateID(LEFT, cellID);
//  const CFuint faceID = cellFaces[f*kd->nbCells + cellID];
//  fd->setIsOutward(kd->isOutward[faceID] == cellID);
//  
//  CFreal* statePtrL = &kd->states[cellID*PHYS::NBEQS];
//  CFreal* nodePtrL = &kd->centerNodes[cellID*PHYS::DIM];
//  setState<PHYS>(fd->getState(LEFT), statePtrL, fd->getNode(LEFT), nodePtrL);
//  setFaceNormal<PHYS>(fd, &kd->normals[faceID*PHYS::DIM]);
//}

//////////////////////////////////////////////////////////////////////////////

template <typename T, CFuint SIZE>
void print(const std::string& name, T* array) 
{
  CFLog(INFO, name << " = \t");
  for (CFuint i = 0; i < SIZE; ++i) {
    CFLog(INFO, array[i] << " ");
  }
  CFLog(INFO, "\n");
}

//////////////////////////////////////////////////////////////////////////////

//template <typename MODEL>
//HOST_DEVICE void computeFaceCentroid(const CellData::Itr* cell, const CFuint faceIdx, 
//				     const CFreal* nodes, CFreal* midFaceCoord)
//{  
//  CudaEnv::CFVecSlice<CFreal, MODEL::DIM> coord(midFaceCoord);
//  coord = 0.;
//  const CFuint nbFaceNodes = cell->getNbFaceNodes(faceIdx);
//  const CFreal ovNbFaceNodes = 1./(static_cast<CFreal>(nbFaceNodes));
//  for (CFuint n = 0; n < nbFaceNodes; ++n) {
//    const CFuint cellNodeID = cell->getNodeID(faceIdx, n);
//    const CFuint nodeID = cell->getNodeID(faceIdx,n);
//    const CFreal* faceNode = &nodes[nodeID*MODEL::DIM];
//    for (CFuint d = 0; d < MODEL::DIM; ++d) {
//      coord[d] += faceNode[d];
//    }
//  }
//  coord *= ovNbFaceNodes;
//}

//////////////////////////////////////////////////////////////////////////////

//template <typename PHYS, typename POLYREC>
//__global__ void computeGradientsKernel(typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor,
//				       const CFuint nbCells,
//				       CFreal* states, 
//				       CFreal* nodes,
//				       CFreal* centerNodes,
//				       CFreal* ghostStates,
//				       CFreal* ghostNodes,
//				       CFreal* uX,
//				       CFreal* uY,
//				       CFreal* uZ,
//				       CFreal* limiter,
//				       CFreal* updateCoeff, 
//				       CFreal* rhs,
//				       CFreal* normals,
//				       CFint* isOutward,
//				       const CFuint* cellInfo,
//				       const CFuint* cellStencil,
//				       const CFuint* cellFaces,
//				       const CFuint* cellNodes,
//				       const CFint*  neighborTypes,
//				       const Framework::CellConn* cellConn)
//{    
//  // each thread takes care of computing the gradient for one single cell
//  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
//  
//  // __shared__ typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE> s_dcor[32];
//  // typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor = &s_dcor[threadIdx.x];
//  // dcor->init(gdcor);
//  
//  if (cellID < nbCells) {    
//    KernelData<CFreal> kd (nbCells, states, nodes, centerNodes, ghostStates, ghostNodes, updateCoeff, 
//			   rhs, normals, uX, uY, uZ, isOutward);
//    
//    // compute and store cell gradients at once 
//    POLYREC polyRec(dcor);
//    CellData cells(nbCells, cellInfo, cellStencil, cellFaces, cellNodes, neighborTypes, cellConn);
//    CellData::Itr cell = cells.getItr(cellID);
//    polyRec.computeGradients(&states[cellID*PHYS::NBEQS], &centerNodes[cellID*PHYS::DIM], &kd, &cell);
//  }
//}
      
//////////////////////////////////////////////////////////////////////////////

//template <typename PHYS, typename POLYREC, typename LIMITER>
//__global__ void computeLimiterKernel(typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE>* dcol,
//				     typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor,
//				     const CFuint nbCells,
//				     CFreal* states, 
//				     CFreal* nodes,
//				     CFreal* centerNodes,
//				     CFreal* ghostStates,
//				     CFreal* ghostNodes,
//				     CFreal* uX,
//				     CFreal* uY,
//				     CFreal* uZ,
//				     CFreal* limiter,
//				     CFreal* updateCoeff, 
//				     CFreal* rhs,
//				     CFreal* normals,
//				     CFint* isOutward,
//				     const CFuint* cellInfo,
//				     const CFuint* cellStencil,
//				     const CFuint* cellFaces,
//				     const CFuint* cellNodes,
//				     const CFint*  neighborTypes,
//				     const Framework::CellConn* cellConn)
//{    
//  // each thread takes care of computing the gradient for one single cell
//  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
//  
//  // __shared__ typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE> s_dcol[32];
//  // typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE>* dcol = &s_dcol[threadIdx.x];
//  // dcol->init(gdcol);
//  
//  // __shared__ typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE> s_dcor[32];
//  // typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor = &s_dcor[threadIdx.x];
//  // dcor->init(gdcor);
//  
//  if (cellID < nbCells) {    
//    // compute all cell quadrature points at once (size of this array is overestimated)
//    CFreal midFaceCoord[PHYS::DIM*PHYS::DIM*2];
//    
//    CellData cells(nbCells, cellInfo, cellStencil, cellFaces, cellNodes, neighborTypes, cellConn);
//    CellData::Itr cell = cells.getItr(cellID);
//    const CFuint nbFacesInCell = cell.getNbFacesInCell();
//    for (CFuint f = 0; f < nbFacesInCell; ++f) { 
//      computeFaceCentroid<PHYS>(&cell, f, nodes, &midFaceCoord[f*PHYS::DIM]);
//    }
//    
//    // compute cell-based limiter at once
//    KernelData<CFreal> kd (nbCells, states, nodes, centerNodes, ghostStates, ghostNodes, updateCoeff, 
//			   rhs, normals, uX, uY, uZ, isOutward);
//    LIMITER limt(dcol);
//    
//    if (dcor->currRes > dcor->limitRes && (dcor->limitIter > 0 && dcor->currIter < dcor->limitIter)) {	
//      limt.limit(&kd, &cell, &midFaceCoord[0], &limiter[cellID*PHYS::NBEQS]);
//    }
//    else {
//      if (!dcor->freezeLimiter) {
//	// historical modification of the limiter
//	CudaEnv::CFVec<CFreal,PHYS::NBEQS> tmpLimiter;
//	limt.limit(&kd, &cell, &midFaceCoord[0], &tmpLimiter[0]);
//	CFuint currID = cellID*PHYS::NBEQS;
//	for (CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar, ++currID) {
//	  limiter[currID] = min(tmpLimiter[iVar],limiter[currID]);
//	}
//      }
//    }
//  }
//}
  
//////////////////////////////////////////////////////////////////////////////
             
template <typename SCHEME>
__global__ void computeFluxKernel(typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE>* dcop,
				  const CFuint nbCells,
				  CFreal* states, 
				  CFreal* updateCoeff, 
				  CFreal* rhs,
				  const CFuint* cellInfo,
				  const CFuint* cellStencil,
				  const CFuint* cellFaces,
				  const CFuint* cellNodes,
				  const CFint*  neighborTypes,
				  const Framework::CellConn* cellConn)
{    
  // each thread takes care of computing the flux for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  // __shared__ typename SCHEME::BASE::template DeviceConfigOptions<NOTYPE> s_dcof[32];
  // typename SCHEME::BASE::template DeviceConfigOptions<NOTYPE>* dcof = &s_dcof[threadIdx.x];
  // dcof->init(gdcof);
  
  // __shared__ typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE> s_dcor[32];
  // typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor = &s_dcor[threadIdx.x];
  // dcor->init(gdcor);
  
  // __shared__ typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE> s_dcop[32];
  // typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE>* dcop = &s_dcop[threadIdx.x];
  // dcop->init(gdcop);
  
  if (cellID < nbCells) {
    // reset the rhs and update coefficients to 0
    CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> res(&rhs[cellID*SCHEME::MODEL::NBEQS]);
    res = 0.;
    updateCoeff[cellID] = 0.;
    
//    KernelData<CFreal> kd (nbCells, states, nodes, centerNodes, ghostStates, ghostNodes, updateCoeff, 
//			   rhs, normals, uX, uY, uZ, isOutward);
//    
//    // compute and store cell gradients at once 
//    SCHEME fluxScheme(dcof);
//    CFreal midFaceCoord[SCHEME::MODEL::DIM*SCHEME::MODEL::DIM*2];
//    FluxData<typename SCHEME::MODEL> currFd; currFd.initialize();
//    typename SCHEME::MODEL pmodel(dcop);
//    
//    CellData cells(nbCells, cellInfo, cellStencil, cellFaces, cellNodes, neighborTypes, cellConn);
//    CellData::Itr cell = cells.getItr(cellID);
//    
//    // compute the fluxes
//    const CFuint nbFacesInCell = cell.getNbActiveFacesInCell();
//    for (CFuint f = 0; f < nbFacesInCell; ++f) { 
//      const CFint stype = cell.getNeighborType(f);
//      
//      if (stype != 0) { // skip all partition faces
//	// set all flux data for the current face
//	const CFuint stateID = cell.getNeighborID(f);
//	setFluxData(f, stype, stateID, cellID, &kd, &currFd, cellFaces);
//	
//	// compute face quadrature points (centroid)
//	CFreal* faceCenters = &midFaceCoord[f*SCHEME::MODEL::DIM];
//	computeFaceCentroid<typename SCHEME::MODEL>(&cell, f, nodes, faceCenters);
//	
//	// extrapolate solution on quadrature points on both sides of the face
//	polyRec.extrapolateOnFace(&currFd, faceCenters, uX, uY, uZ, limiter);
//	
//	// compute the convective flux across the face
//        fluxScheme.prepareComputation(&currFd, &pmodel);
//	fluxScheme(&currFd, &pmodel);
//	
//	// update the residual
//	CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> ress(currFd.getResidual());
//	res -= ress;
//	
//	// update the update coefficient
//	updateCoeff[cellID] += currFd.getUpdateCoeff();
//      }
//    }
//
//    // get the states in this cell
//    m_cellStates = m_cell->getStates();
//      
//    // set the cell data
//    setCellData();
//      
//    // compute the divergence of the discontinuous flux (-divFD+divhFD)
//    computeDivDiscontFlx(m_divContFlx);
//      
//    // update RHS
//    updateRHS();
     
  }
}
  
//////////////////////////////////////////////////////////////////////////////

//template <typename SCHEME, typename POLYREC, typename LIMITER>
//void computeFluxCPU(CFuint nbThreadsOMP,
//		    typename SCHEME::BASE::template DeviceConfigOptions<NOTYPE>* dcof,
//		    typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor,
//		    typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE>* dcol,
//		    typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE>* dcop,
//		    const CFuint nbCells,
//		    CFreal* states, 
//		    CFreal* nodes,
//		    CFreal* centerNodes,
//		    CFreal* ghostStates,
//		    CFreal* ghostNodes,
//		    CFreal* uX,
//		    CFreal* uY,
//		    CFreal* uZ,
//		    CFreal* limiter,
//		    CFreal* updateCoeff, 
//		    CFreal* rhs,
//		    CFreal* normals,
//		    CFint* isOutward,
//		    const CFuint* cellInfo,
//		    const CFuint* cellStencil,
//		    const CFuint* cellFaces,
//		    const CFuint* cellNodes,
//		    const CFint* neighborTypes,
//		    const Framework::CellConn* cellConn)
//{ 
//  typedef typename SCHEME::MODEL PHYS;
//  
//  FluxData<PHYS> fd;
//#ifndef CF_HAVE_OMP  
//  fd.initialize();
//  FluxData<PHYS>* currFd = &fd;
//  cf_assert(currFd != CFNULL);
//#endif
//  POLYREC polyRec(dcor);
//  SCHEME fluxScheme(dcof);
//  LIMITER limt(dcol);
//  PHYS pmodel(dcop);
//  
//  CellData cells(nbCells, cellInfo, cellStencil, cellFaces, cellNodes, neighborTypes, cellConn);
//  KernelData<CFreal> kd(nbCells, states, nodes, centerNodes, ghostStates, ghostNodes, updateCoeff, 
//			rhs, normals, uX, uY, uZ, isOutward);
//  
//  CFreal midFaceCoord[PHYS::DIM*PHYS::DIM*2];
//  CudaEnv::CFVec<CFreal,PHYS::NBEQS> tmpLimiter;
//
//#ifdef CF_HAVE_OMP
//  //const CFuint nThr = omp_get_num_procs();
//  // omp_set_num_threads(nbThreadsOMP);
//#pragma omp num_thread(nbThreadsOMP) parallel private(polyRec) private(fd)
//{
//  #pragma omp for
//#endif 
//  // compute the cell-based gradients
//  for (CFuint cellID = 0; cellID < nbCells; ++cellID) {
//#ifdef CF_HAVE_OMP
//    fd.initialize();
//    FluxData<PHYS>* currFd = &fd;
//    cf_assert(currFd != CFNULL);
//#endif 
//    CellData::Itr cell = cells.getItr(cellID);
//    polyRec.computeGradients(&states[cellID*PHYS::NBEQS], &centerNodes[cellID*PHYS::DIM], &kd, &cell);
//  }
//#ifdef CF_HAVE_OMP
//}
//#endif
//
//#ifdef CF_HAVE_OMP  
//#pragma omp num_thread(nbThreadsOMP) parallel private(limt) private(kd)
//{
//  #pragma omp for
//#endif 
//  // compute the cell based limiter 
//  for (CFuint cellID = 0; cellID < nbCells; ++cellID) {
//  // for (CellData::Itr cell = cells.begin(); cell <= cells.end(); ++cell) {
//    CellData::Itr cell = cells.getItr(cellID);
//    // compute all cell quadrature points at once (size of this array is overestimated)
//    const CFuint nbFacesInCell = cell.getNbFacesInCell();
//    for (CFuint f = 0; f < nbFacesInCell; ++f) { 
//      computeFaceCentroid<PHYS>(&cell, f, nodes, &midFaceCoord[f*PHYS::DIM]);
//    }
//    
//    //   const CFuint cellID = cell.getCellID();
//    if (dcor->currRes > dcor->limitRes && (dcor->limitIter > 0 && dcor->currIter < dcor->limitIter)) {	
//      // compute cell-based limiter
//      limt.limit(&kd, &cell, &midFaceCoord[0], &limiter[cellID*PHYS::NBEQS]);
//    }
//    else {
//      if (!dcor->freezeLimiter) {
//	// historical modification of the limiter
//	limt.limit(&kd, &cell, &midFaceCoord[0], &tmpLimiter[0]);
//	CFuint currID = cellID*PHYS::NBEQS;
//	for (CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar, ++currID) {
//	  limiter[currID] = min(tmpLimiter[iVar],limiter[currID]);
//	}
//      }
//    }
//  }
//#ifdef CF_HAVE_OMP
//}
//
//#pragma omp num_thread(nbThreadsOMP) parallel private(fd) private(kd) private(fluxScheme) private(pmodel)
//{
//  #pragma omp for
//#endif 
//  // compute the fluxes
//  for (CFuint cellID = 0; cellID < nbCells; ++cellID) {
//  //  for (CellData::Itr cell = cells.begin(); cell <= cells.end(); ++cell) {
//#ifdef CF_HAVE_OMP
//    fd.initialize();
//    FluxData<PHYS>* currFd = &fd;
//    cf_assert(currFd != CFNULL);
//#endif
//    // reset the rhs and update coefficients to 0
//   // const CFuint cellID = cell.getCellID();
//    CudaEnv::CFVecSlice<CFreal,PHYS::NBEQS> res(&rhs[cellID*PHYS::NBEQS]);
//    res = 0.;
//    updateCoeff[cellID] = 0.;
//
//    CellData::Itr cell = cells.getItr(cellID);   
//    const CFuint nbFacesInCell = cell.getNbActiveFacesInCell();
//    for (CFuint f = 0; f < nbFacesInCell; ++f) { 
//      const CFint stype = cell.getNeighborType(f);
//      
//      if (stype != 0) { // skip all partition faces
//	const CFuint stateID =  cell.getNeighborID(f);
//	setFluxData(f, stype, stateID, cellID, &kd, currFd, cellFaces);
//	
//	// compute face quadrature points (centroid)
//	CFreal* faceCenters = &midFaceCoord[f*PHYS::DIM];
//	computeFaceCentroid<PHYS>(&cell, f, nodes, faceCenters);
//	
//	// extrapolate solution on quadrature points on both sides of the face
//	polyRec.extrapolateOnFace(currFd, faceCenters, uX, uY, uZ, limiter);
//        fluxScheme.prepareComputation(currFd, &pmodel);
//	fluxScheme(currFd, &pmodel); // compute the convective flux across the face
//	
//	for (CFuint iEq = 0; iEq < PHYS::NBEQS; ++iEq) {
//	  const CFreal value = currFd->getResidual()[iEq];
//	  res[iEq] -= value;  // update the residual 
//	}
//	
//	// update the update coefficient
//	updateCoeff[cellID] += currFd->getUpdateCoeff();
//      }
//    }
//  }
//#ifdef CF_HAVE_OMP
//} 
//#endif
//}

//////////////////////////////////////////////////////////////////////////////

template <typename SCHEME, typename PHYSICS, CFuint NB_BLOCK_THREADS>
void ConvRHSFluxReconstructionCUDA<SCHEME,PHYSICS,NB_BLOCK_THREADS>::execute()
{
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  CFTRACEBEGIN;
  
  CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() START\n");
  
  initializeComputationRHS();

  const CFuint nbCells = socket_states.getDataHandle().size();
  cf_assert(nbCells > 0);
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle(); 
  
  SafePtr<SCHEME> lf  = getMethodData().getRiemannFlux().d_castTo<SCHEME>();
  SafePtr<typename PHYSICS::PTERM> phys = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<typename PHYSICS::PTERM>();
  
#ifdef CF_HAVE_CUDA
  typedef typename SCHEME::template DeviceFunc<GPU, PHYSICS> FluxScheme;  
#else
  typedef typename SCHEME::template DeviceFunc<CPU, PHYSICS> FluxScheme;
#endif 
  
  if (m_onGPU) {
#ifdef CF_HAVE_CUDA

    CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
    timer.start();
    
    // copy of data that change at every iteration
    socket_states.getDataHandle().getGlobalArray()->put(); 
     
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => CPU-->GPU data transfer took " << timer.elapsed() << " s\n");
    timer.start();
    
    ConfigOptionPtr<SCHEME,  NOTYPE, GPU> dcof(lf);
    ConfigOptionPtr<typename PHYSICS::PTERM, NOTYPE, GPU> dcop(phys);
    
    const CFuint blocksPerGrid = CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
    const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
    
    //dim3 blocks(m_nbBlocksPerGridX, m_nbBlocksPerGridY);
    
    //cudaFuncSetCacheConfig("computeGradientsKernel", cudaFuncCachePreferL1);
    
    // cudaFuncSetCacheConfig("computeFluxKernel", cudaFuncCachePreferL1);
    
    // compute the convective flux in each cell
    computeFluxKernel<FluxScheme> <<<blocksPerGrid,nThreads>>> 
      (dcop.getPtr(),
       nbCells,
       socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
       updateCoeff.getLocalArray()->ptrDev(), 
       rhs.getLocalArray()->ptrDev(),
       m_cellInfo.ptrDev(),
       m_cellStencil.ptrDev(),
       m_cellFaces->getPtr()->ptrDev(),
       m_neighborTypes.ptrDev(),
       m_cellConn.ptrDev());
    
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => computeFluxKernel took " << timer.elapsed() << " s\n");
    
    timer.start();
    rhs.getLocalArray()->get();
    updateCoeff.getLocalArray()->get();
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => GPU-->CPU data transfer took " << timer.elapsed() << " s\n");

#endif
}
  else {
    // AL: useful fo debugging
    // for (CFuint i = 0; i <  m_ghostStates.size()/9; ++i) {
    //   std::cout.precision(12); std::cout << "g" << i << " => ";
    //   for (CFuint j = 0; j < 9; ++j) {
    // 	std::cout << m_ghostStates[i*9+j] << " ";
    //   }
    //   std::cout << "\n";
    // }
    // for (CFuint i = 0; i <  socket_states.getDataHandle().size(); ++i) {
    //   std::cout.precision(12); std::cout << i << " => "<< *socket_states.getDataHandle()[i] <<"\n";
    // }
    
//    ConfigOptionPtr<SCHEME>  dcof(lf);
//    ConfigOptionPtr<POLYREC> dcor(pr);
//    ConfigOptionPtr<LIMITER> dcol(lm);
//    ConfigOptionPtr<typename PHYSICS::PTERM> dcop(phys);
//    
//    computeFluxCPU<FluxScheme, PolyRec, Limiter>
//      (m_nbThreadsOMP,
//       dcof.getPtr(),
//       dcor.getPtr(),
//       dcol.getPtr(),
//       dcop.getPtr(),
//       nbCells,
//       socket_states.getDataHandle().getGlobalArray()->ptr(), 
//       socket_nodes.getDataHandle().getGlobalArray()->ptr(),
//       m_centerNodes.ptr(), 
//       m_ghostStates.ptr(),
//       m_ghostNodes.ptr(),
//       socket_uX.getDataHandle().getLocalArray()->ptr(),
//       socket_uY.getDataHandle().getLocalArray()->ptr(),
//       socket_uZ.getDataHandle().getLocalArray()->ptr(),
//       socket_limiter.getDataHandle().getLocalArray()->ptr(),
//       updateCoeff.getLocalArray()->ptr(), 
//       rhs.getLocalArray()->ptr(),
//       normals.getLocalArray()->ptr(),
//       isOutward.getLocalArray()->ptr(),
//       m_cellInfo.ptr(),
//       m_cellStencil.ptr(),
//       m_cellFaces->getPtr()->ptr(),
//       m_cellNodes->getPtr()->ptr(),
//       m_neighborTypes.ptr(),
//       m_cellConn.ptr());
  }
  
// for (int i = 0; i < updateCoeff.size(); ++i) {
//      std::cout << "updateCoeff[" << i << "] = " << updateCoeff[i]  << std::endl;
//       /* std::cout << "rhs[" << i << "] = ";
//        for (int j = 0; j < 9; ++j) {
//          std::cout << rhs[i*9+j] << " ";
//        }
//        std::cout << std::endl;*/
// } 
//   abort();
  // for (;;) {}
  
  //finalizeComputationRHS();
  
  CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() END\n");
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
