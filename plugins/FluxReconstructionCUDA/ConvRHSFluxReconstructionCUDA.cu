#include "FluxReconstructionCUDA/ConvRHSFluxReconstructionCUDA.hh"
#include "Framework/MeshData.hh"
#include "Framework/CellConn.hh"
#include "Config/ConfigOptionPtr.hh"
#include "Framework/CudaDeviceManager.hh"
#include "Common/CUDA/CFVec.hh"
#include "Framework/CudaTimer.hh"

#include "FluxReconstructionMethod/FluxData.hh"
#include "FluxReconstructionMethod/KernelData.hh"
#include "FluxReconstructionMethod/CellData.hh"

#include "FluxReconstructionCUDA/FluxReconstructionCUDA.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/VarSetListT.hh"
#include "NavierStokes/Euler2DVarSetT.hh"
#include "NavierStokes/Euler2DConsT.hh"

#include "FluxReconstructionMethod/LaxFriedrichsFlux.hh"
#include <stdio.h>

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
                      VarSetListT<Euler##__dim__##__svars__##T, Euler##__dim__##__uvars__##T>, __nbBThreads__>, \
		      FluxReconstructionSolverData,FluxReconstructionCUDAModule>	\
FR_RhsEuler##__dim__##__svars__##__uvars__##__nbBThreads__##Provider(__providerName__);
// 48 block threads (default)
FR_EULER_RHS_PROV(2D, Cons, Cons, 48, "EulerFRLaxFriedrichs2DCons")
//FR_EULER_RHS_PROV(3D, Cons, Cons, 48, "EulerFRLaxFried3DCons")
//FR_NS_RHS_PROV(2D, ProjectionCons, ProjectionPrim, 48, "CellLaxFriedMHD2DPrim")
//FR_NS_RHS_PROV(3D, ProjectionCons, ProjectionPrim, 48, "CellLaxFriedMHD3DPrim")
#undef FR_EULER_RHS_PROV

//////////////////////////////////////////////////////////////////////////////

template <typename PHYS>
HOST_DEVICE inline void setState(CFreal* state, CFreal* statePtr)
{
  // copy the state node data to shared memory
  //for (CFuint i = 0; i < PHYS::DIM; ++i) {node[i] = nodePtr[i];}
  // copy the state data to shared memory
  for (CFuint i = 0; i < PHYS::NBEQS; ++i) {state[i] = statePtr[i];} 
}
      
//////////////////////////////////////////////////////////////////////////////

template <typename PHYS>
HOST_DEVICE void setFluxData(const CFuint stateID, const CFuint cellID, 
			     KernelData<CFreal>* kd, FluxData<PHYS>* fd, const CFuint iSol)
{
  fd->setStateID(LEFT, stateID);
  CFreal* statePtrR = &kd->states[stateID*PHYS::NBEQS];  

  setState<PHYS>(fd->getState(iSol), statePtrR);

  fd->setNbSolPnts(kd->nbSolPnts);
}

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

template <typename SCHEME, typename PHYS>
__global__ void computeStateLocalRHSKernel(typename SCHEME::BASE::template DeviceConfigOptions<NOTYPE>* dcof,
                                  typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE>* dcop,
                                  const CFuint nbCells,
				  CFreal* states, 
                                  CFreal* updateCoeff, 
				  CFreal* rhs,
                                  CFreal* solPntNormals,
                                  const CFuint nbSolPnts,
                                  const CFuint nbrFaces,
				  const CFuint* cellInfo,
                                  const CFuint* stateIDs,
                                  const CFuint* neighbCellIDs,
                                  const CFuint dim,
                                  const CFuint nbrEqs,
                                  const CFuint nbrFlxPnts,
                                  const CFuint nbrSolSolDep,
                                  const CFuint* solSolDep,
                                  const CFuint nbrSolFlxDep,
                                  const CFuint* solFlxDep,
                                  const CFreal* solPolyDerivAtSolPnts,
                                  const CFreal* solPolyValsAtFlxPnts,
                                  const CFuint* flxPntFlxDim,
                                  const CFreal* corrFctDiv)
{    
  // one thread per cell
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
  
  if (cellID < nbCells) 
  { 
    // current kernel data
    KernelData<CFreal> kd (nbCells, states, updateCoeff, rhs, solPntNormals, nbSolPnts);

    // current flux data
    FluxData<typename SCHEME::MODEL> currFd; 
    
    // initialize flux data
    currFd.initialize();
    
    // physical model
    typename SCHEME::MODEL pmodel(dcop);
    SCHEME fluxScheme(dcof);
    
    // current cell data
    CellData cells(nbCells, cellInfo, stateIDs, neighbCellIDs, nbrFaces, nbSolPnts);
    
    // get current cell
    CellData::Itr cell = cells.getItr(cellID);
          
    const CFuint nbFlxPntFlx = SCHEME::MODEL::NBEQS*8;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntFlx;
    
    flxPntFlx = 0.0;

    // loop over sol pnts to compute flux
    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);
      //printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
    
    if (cellID == 0) printf("GPUstate: %f %f %f %f\n", kd.states[0], kd.states[1], kd.states[2], kd.states[3]);
    
      setFluxData(stateID, cellID, &kd, &currFd, iSolPnt);

      const CFuint nbNormals = PHYS::DIM*PHYS::DIM;

      CudaEnv::CFVecSlice<CFreal,nbNormals> n(&(kd.solPntNormals[stateID*nbNormals]));

      CudaEnv::CFVecSlice<CFreal,nbNormals> nFd(currFd.getScaledNormal(iSolPnt));
      
      for (CFuint i = 0; i < nbNormals; ++i) 
      {
        nFd[i] = n[i];
      }
      
      // get the flux
      fluxScheme.prepareComputation(&currFd, &pmodel);
      
      fluxScheme(&currFd, &pmodel, false, iSolPnt);
      
      // loop over sol pnts to compute flux
      for (CFuint iDim = 0; iDim < dim; ++iDim)
      {
        if (cellID == 0) printf("HERE4 iSol: %d, iDim: %d, flux: %f %f %f %f \n", iSolPnt, iDim, currFd.getFlux(iSolPnt, iDim)[0], currFd.getFlux(iSolPnt, iDim)[1], currFd.getFlux(iSolPnt, iDim)[2], currFd.getFlux(iSolPnt, iDim)[3]);
      }
      
      // Loop over solution pnts to count the factor of all sol pnt polys
      for (CFuint jSolPnt = 0; jSolPnt < nbrSolSolDep; ++jSolPnt)
      { 
      
        const CFuint jSolIdx = solSolDep[iSolPnt*nbrSolSolDep+jSolPnt]; //(*m_solSolDep)[iSolPnt][jSolPnt];
        
        // get current vector slice out of rhs
        CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> res(&rhs[stateID*SCHEME::MODEL::NBEQS]);

        // Loop over deriv directions and sum them to compute divergence
        for (CFuint iDir = 0; iDir < dim; ++iDir)
        {
          const CFreal polyCoef = solPolyDerivAtSolPnts[iSolPnt*dim*nbSolPnts+iDir*nbSolPnts+jSolIdx];//(*m_solPolyDerivAtSolPnts)[jSolPnt][iDir][iSolIdx]; 
          
          if (cellID == 0) printf("polyCoef: %f\n", polyCoef);
          
          // Loop over conservative fluxes 
          for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
          {
            // Store divFD in the vector that will be divFC
            res[iEq] -= polyCoef*(currFd.getFlux(iSolPnt, iDir)[iEq]);
            if (cellID == 0) printf("res %f \n", res[iEq]);
	  }
        }
      }
      
      // extrapolate the fluxes to the flux points
      for (CFuint iFlxPnt = 0; iFlxPnt < nbrSolFlxDep; ++iFlxPnt)
      {
        const CFuint flxIdx = solFlxDep[iSolPnt*nbrSolFlxDep+iFlxPnt];
        const CFuint dim = flxPntFlxDim[flxIdx];
        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
        {
          flxPntFlx[flxIdx*nbrEqs+iEq] += solPolyValsAtFlxPnts[flxIdx*nbrFlxPnts+iSolPnt]*currFd.getFlux(iSolPnt, dim)[iEq];
        }
      }
    }
    
    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);
      
      // get current vector slice out of rhs
      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> res(&rhs[stateID*SCHEME::MODEL::NBEQS]);
        
      // add divhFD to the residual updates
      for (CFuint iFlxPnt = 0; iFlxPnt < nbrSolFlxDep; ++iFlxPnt)
      {
        const CFuint flxIdx = solFlxDep[iSolPnt*nbrSolFlxDep+iFlxPnt];

        // get the divergence of the correction function
        const CFreal divh = corrFctDiv[iSolPnt*nbrFlxPnts+iFlxPnt];
 
        // Fill in the corrections
        for (CFuint iVar = 0; iVar < nbrEqs; ++iVar)
        {
          res[iVar] += flxPntFlx[flxIdx*nbrEqs+iVar] * divh; 
        }
      }
    }
    
    
    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      // reset flx pnt fluxes  
      flxPntFlx = 0.0;
        
      // loop over sol pnts to compute flux
      for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
      {
        // get current state ID
        const CFuint neighbStateID = cell.getNeighbStateID(iFace,iSolPnt);
        
        const CFuint neighbCellID = cell.getNeighbCellID(iFace);
    
        setFluxData(neighbStateID, neighbCellID, &kd, &currFd, iSolPnt);

        const CFuint nbNormals = PHYS::DIM*PHYS::DIM;

        CudaEnv::CFVecSlice<CFreal,nbNormals> n(&(kd.solPntNormals[neighbStateID*nbNormals]));

        CudaEnv::CFVecSlice<CFreal,nbNormals> nFd(currFd.getScaledNormal(iSolPnt));
      
        for (CFuint i = 0; i < nbNormals; ++i) 
        {
          nFd[i] = n[i];
        }
      
        // get the flux
        fluxScheme.prepareComputation(&currFd, &pmodel);
      
        fluxScheme(&currFd, &pmodel, true, iSolPnt);
      
        // extrapolate the fluxes to the flux points
        for (CFuint iFlxPnt = 0; iFlxPnt < nbrSolFlxDep; ++iFlxPnt)
        {
          const CFuint flxIdx = solFlxDep[iSolPnt*nbrSolFlxDep+iFlxPnt];
          const CFuint dim = flxPntFlxDim[flxIdx];
          // Loop over conservative fluxes 
          for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
          {
            flxPntFlx[flxIdx*nbrEqs+iEq] += solPolyValsAtFlxPnts[flxIdx*nbrFlxPnts+iSolPnt]*currFd.getFlux(iSolPnt, dim)[iEq];
          }
        }
      }
    }
      
    
    if (cellID == 0) printf("end sol pnt \n",cellID);
    
//
//  // get residual factor
//  const CFreal resFactor = getMethodData().getResFactor();
//
//  // update rhs
//  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//  {
//    CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
//    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//    {
//      rhs[resID+iVar] += resFactor*m_divContFlx[iState][iVar];
//    }
//  }
    
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
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;
  
  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;
  
  // loop over element types, for the moment there should only be one
  const CFuint nbrElemTypes = elemType->size();
  cf_assert(nbrElemTypes == 1);
  
  // get start and end indexes for this type of element
  cf_assert((*elemType)[0].getStartIdx() == 0);
  const CFuint nbCells   = (*elemType)[0].getEndIdx();
  cf_assert(nbCells > 0);
  
  initializeComputationRHS();

  const CFuint nbStates = socket_states.getDataHandle().size();
  cf_assert(nbStates > 0);
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle(); 
  DataHandle<CFreal> solPntNormals = socket_solPntNormals.getDataHandle(); 
  
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
    socket_rhs.getDataHandle().getLocalArray()->put(); 
CFLog(INFO, "nb normals: " << socket_solPntNormals.getDataHandle().size() << ", n0: " << socket_solPntNormals.getDataHandle()[0] << "\n");
socket_solPntNormals.getDataHandle().getLocalArray()->put();
    DataHandle<Framework::State*, Framework::GLOBAL > statesI = socket_states.getDataHandle();
    
    CFLog(INFO, "state: " << *(statesI[0]) << "\n");
     
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => CPU-->GPU data transfer took " << timer.elapsed() << " s\n");
    timer.start();
    
    ConfigOptionPtr<SCHEME,  NOTYPE, GPU> dcof(lf);
    ConfigOptionPtr<typename PHYSICS::PTERM, NOTYPE, GPU> dcop(phys);

    const CFuint blocksPerGrid = CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
    const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
    CFLog(VERBOSE, "blocksPerGrid: " << blocksPerGrid << ", threads: " << nThreads << "\n");

    //dim3 blocks(m_nbBlocksPerGridX, m_nbBlocksPerGridY);
    
    //cudaFuncSetCacheConfig("computeGradientsKernel", cudaFuncCachePreferL1);
    
    // cudaFuncSetCacheConfig("computeFluxKernel", cudaFuncCachePreferL1);

    // compute the convective flux in each cell
    computeStateLocalRHSKernel<FluxScheme,PHYSICS> <<<blocksPerGrid,nThreads>>> 
      (dcof.getPtr(),
       dcop.getPtr(),
       nbCells,
       socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
       updateCoeff.getLocalArray()->ptrDev(), 
       rhs.getLocalArray()->ptrDev(),
       solPntNormals.getLocalArray()->ptrDev(),
       m_nbrSolPnts,
       4,
       m_cellInfo.ptrDev(),
       m_stateIDs.ptrDev(),
       m_neighbCellIDs.ptrDev(),
       m_dim,
       m_nbrEqs,
       m_nbrFlxPnts,
       m_nbrSolSolDep,
       m_solSolDep2.ptrDev(),
       m_nbrFlxDep,
       m_solFlxDep2.ptrDev(),
       m_solPolyDerivAtSolPnts2.ptrDev(),
       m_solPolyValsAtFlxPnts2.ptrDev(),
       m_flxPntFlxDim2.ptrDev(),
       m_corrFctDiv2.ptrDev());
    cudaDeviceSynchronize();
    
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => computeFluxKernel took " << timer.elapsed() << " s\n");
    
    timer.start();
    rhs.getLocalArray()->get();
    updateCoeff.getLocalArray()->get();
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => GPU-->CPU data transfer took " << timer.elapsed() << " s\n");

#endif
}
  else {
      for (CFuint i = 0; i < nbCells; i++)
      {
          CFreal ID = i*5+i*20;
      }
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