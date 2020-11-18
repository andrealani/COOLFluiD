#include "FluxReconstructionCUDA/ConvDiffLLAVRHSFluxReconstructionCUDA.hh"
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
#include "Framework/MathTypes.hh"
#include "NavierStokes/Euler2DVarSetT.hh"
#include "NavierStokes/Euler2DConsT.hh"
#include "NavierStokes/NavierStokes2DVarSetT.hh"
#include "NavierStokes/NavierStokes2DConsT.hh"
#include "NavierStokes/NSVarSetListT.hh"

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

#define FR_NSLLAV_RHS_PROV(__dim__,__svars__,__uvars__,__order__,__nbBThreads__,__providerName__) \
MethodCommandProvider<ConvDiffLLAVRHSFluxReconstructionCUDA<LaxFriedrichsFlux, \
                      VarSetListT<Euler##__dim__##__svars__##T, Euler##__dim__##__uvars__##T>,NSVarSetListT<NavierStokes##__dim__##__svars__##T, NavierStokes##__dim__##__uvars__##T>,__order__,__nbBThreads__>, \
		      FluxReconstructionSolverData,FluxReconstructionCUDAModule>	\
FR_RhsNSLLAV##__dim__##__svars__##__uvars__##__order__##__nbBThreads__##Provider(__providerName__);
// 48 block threads (default)
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 0, 48, "NSLLAVFRLaxFriedrichs2DConsP0")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 1, 48, "NSLLAVFRLaxFriedrichs2DConsP1")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 2, 48, "NSLLAVFRLaxFriedrichs2DConsP2")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 3, 48, "NSLLAVFRLaxFriedrichs2DConsP3")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 4, 48, "NSLLAVFRLaxFriedrichs2DConsP4")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 5, 48, "NSLLAVFRLaxFriedrichs2DConsP5")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 6, 48, "NSLLAVFRLaxFriedrichs2DConsP6")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 7, 48, "NSLLAVFRLaxFriedrichs2DConsP7")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 8, 48, "NSLLAVFRLaxFriedrichs2DConsP8")
FR_NSLLAV_RHS_PROV(2D, Cons, Cons, 9, 48, "NSLLAVFRLaxFriedrichs2DConsP9")
//FR_EULER_RHS_PROV(3D, Cons, Cons, 48, "EulerFRLaxFried3DCons")
//FR_NS_RHS_PROV(2D, ProjectionCons, ProjectionPrim, 48, "CellLaxFriedMHD2DPrim")
//FR_NS_RHS_PROV(3D, ProjectionCons, ProjectionPrim, 48, "CellLaxFriedMHD3DPrim")
#undef FR_NS_RHS_PROV

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

template <typename PHYS, CFuint ORDER>
HOST_DEVICE void setFluxData(const CFuint stateID, const CFuint cellID, 
			     KernelData<CFreal>* kd, FluxData<PHYS,ORDER>* fd, const CFuint iSol)
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

template <typename SCHEME, typename PHYS, typename PHYSNS, CFuint ORDER>
__global__ void computeStateLocalRHSKernel(typename SCHEME::BASE::template DeviceConfigOptions<NOTYPE>* dcof,
                                  typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE>* dcop,
                                  typename PHYSNS::DTERM::template DeviceConfigOptions<NOTYPE>* dcopNS,
                                  typename PHYSNS::PTERM::template DeviceConfigOptions<NOTYPE>* dcopNSConv,
                                  const CFuint nbCells,
                                  const CFreal resFactor,
				  CFreal* states, 
                                  CFreal* gradients,
                                  CFreal* gradientsAV,
                                  CFreal* updateCoeff, 
				  CFreal* rhs,
                                  CFreal* solPntNormals,
                                  CFreal* flxPntNormals,
                                  CFreal* cellVolumes,
                                  CFint* faceDir,
                                  const CFuint nbSolPnts,
                                  const CFuint nbrFaces,
                                  const CFuint* faceFlxPntConn,
                                  const CFuint* stateIDs,
                                  const CFint* neighbCellIDs,
                                  const CFuint* neighbFaceIDs,
                                  const CFuint* innerCellIsLeft,
                                  const CFuint nbrFlxPnts,
                                  const CFuint nbrSolSolDep,
                                  const CFuint* solSolDep,
                                  const CFuint nbrSolFlxDep,
                                  const CFuint* solFlxDep,
                                  const CFuint nbrFlxSolDep,
                                  const CFuint* flxSolDep,
                                  const CFreal* solPolyDerivAtSolPnts,
                                  const CFreal* solPolyValsAtFlxPnts,
                                  const CFuint* flxPntFlxDim,
                                  const CFreal* corrFctDiv,
                                  const CFreal* faceIntCoeff,
                                  const CFreal cflConvDiffRatio,
                                  const CFuint* nbNodeNeighbors,
                                  const CFreal* nodeEpsilons,
                                  const CFuint nbrCornerNodes,
                                  const CFuint* neighbNodeIDs,
                                  const CFuint* faceNeighbNodeIDs,
                                  const CFuint nbFaceNodes,
                                  const CFreal* nodePolyValsAtFlxPnts,
                                  const CFreal* nodePolyValsAtSolPnts,
                                  const bool addUpdCoeff)
{    
  // one thread per cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;

  if (cellID < nbCells) 
  { 
    // current kernel data
    KernelData<CFreal> kd (nbCells, states, updateCoeff, rhs, solPntNormals, flxPntNormals, faceDir, nbSolPnts);

    // current flux data
    FluxData<typename SCHEME::MODEL,ORDER> currFd; 

    // initialize flux data
    currFd.initialize();
    
    // physical model
    typename SCHEME::MODEL pmodel(dcop);
    SCHEME fluxScheme(dcof);
    
    PHYSNS pmodelNS(dcopNS,dcopNSConv);    
    
    // current cell data
    CellData cells(nbCells, stateIDs, neighbCellIDs, neighbFaceIDs, innerCellIsLeft, nbrFaces, nbSolPnts, ORDER);
    
    // get current cell
    CellData::Itr cell = cells.getItr(cellID);
          
    // initialize constants and vectors
    const CFuint nbFlxPntFlx = SCHEME::MODEL::NBEQS*(ORDER+1)*2*PHYS::DIM;//8;
    
    //const CFuint nbFaceFlxPntFlx = SCHEME::MODEL::NBEQS*(ORDER+1);
   
    const CFuint nbrFaceFlxPnts = (ORDER+1);

    const CFuint totNbrFlxPnts = (ORDER+1)*2*PHYS::DIM;

    const CFuint nbNormals = PHYS::DIM*PHYS::DIM;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntFlx;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSol;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx*PHYS::DIM> flxPntGrads;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx*PHYS::DIM> flxPntGradsAV;
    
    CudaEnv::CFVec<CFreal,SCHEME::MODEL::NBEQS*(ORDER+1)*(ORDER+1)*PHYS::DIM> solPntFlx;
    
    CudaEnv::CFVec<CFreal,SCHEME::MODEL::NBEQS> avgSol;
    
    CudaEnv::CFVec<CFreal,SCHEME::MODEL::NBEQS*PHYS::DIM> avgGrad;
    
    CudaEnv::CFVec<CFreal,SCHEME::MODEL::NBEQS*PHYS::DIM> avgGradAV;    
    
    CudaEnv::CFVec<CFreal,SCHEME::MODEL::NBEQS> currFlxPntFlx;
    
    CudaEnv::CFVec<CFreal,(ORDER+1)*(ORDER+1)> solEpsilons;    
    
    flxPntFlx = 0.0;
    
    flxPntSol = 0.0;
    
    flxPntGrads = 0.0;
    
    flxPntGradsAV = 0.0;    
    
    solPntFlx = 0.0;

    avgSol = 0.0;
    
    avgGrad = 0.0;
    
    avgGradAV = 0.0;
    
    currFlxPntFlx = 0.0;
    
    const CFreal currVol = cellVolumes[cellID];
    
    solEpsilons = 0.0;

    // loop over flx pnts to extrapolate the states to the flux points
    for (CFuint iSol = 0; iSol < nbSolPnts; ++iSol)
    {   
      // loop over the sol pnts to compute the states and grads in the flx pnts
      for (CFuint iNode = 0; iNode < nbrCornerNodes; ++iNode)
      {
        // get node local index
        const CFuint nodeIdx = neighbNodeIDs[cellID*nbrCornerNodes+iNode];
      
        solEpsilons[iSol] += nodePolyValsAtSolPnts[iSol*nbrCornerNodes+iNode]*nodeEpsilons[nodeIdx]/nbNodeNeighbors[nodeIdx];
      }
      //if (cellID == 0) printf("eps %e\n",solEpsilons[iSol]);
    }
    
    // loop over sol pnts to compute flux
    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);
    
      setFluxData(stateID, cellID, &kd, &currFd, iSolPnt);

      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS*PHYS::DIM> grad(&gradients[stateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);
      
      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS*PHYS::DIM> gradAV(&gradientsAV[stateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);

      CudaEnv::CFVecSlice<CFreal,nbNormals> n(&(kd.solPntNormals[stateID*nbNormals]));

      CudaEnv::CFVecSlice<CFreal,nbNormals> nFd(currFd.getScaledNormal(iSolPnt));
      
      for (CFuint i = 0; i < nbNormals; ++i) 
      {
        nFd[i] = n[i];
      }
          
      for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
      {
        CFreal nJacob2 = 0.0;

        for (CFuint jDir = 0; jDir < PHYS::DIM; ++jDir)
        {
          nJacob2 += n[iDir*PHYS::DIM+jDir]*n[iDir*PHYS::DIM+jDir];
        }
      }
      
      // get the flux
      fluxScheme.prepareComputation(&currFd, &pmodel);

      fluxScheme(&currFd, &pmodel, iSolPnt);
      
      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[stateID*SCHEME::MODEL::NBEQS]);
            
      for (CFuint iDim = 0; iDim < PHYS::DIM; ++iDim)
      {
        pmodelNS.getUpdateVS()->getFlux(&currState[0],&grad[0],&n[iDim*PHYS::DIM],&solPntFlx[iSolPnt*SCHEME::MODEL::NBEQS*PHYS::DIM+iDim*SCHEME::MODEL::NBEQS]);
        
        for (CFuint iDim2 = 0; iDim2 < PHYS::DIM; ++iDim2)
        {
          for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
          {  
          
            solPntFlx[iSolPnt*SCHEME::MODEL::NBEQS*PHYS::DIM+iDim*SCHEME::MODEL::NBEQS+iEq] += gradAV[iEq*PHYS::DIM+iDim2]*solEpsilons[iSolPnt]*n[iDim*PHYS::DIM+iDim2];
            
            //if (cellID == 0) printf("first iSol: %d, iDir: %d, iEq: %d, gradAV: %e, eps: %e, n: %e\n", iSolPnt, iDim, iEq, gradAV[iEq*PHYS::DIM+iDim2],solEpsilons[iSolPnt],n[iDim*PHYS::DIM+iDim2]);
          }
        }
      }
    }

    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);

      setFluxData(stateID, cellID, &kd, &currFd, iSolPnt);

      // get current vector slice out of rhs
      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> res(&rhs[stateID*SCHEME::MODEL::NBEQS]);

      // Loop over solution pnts to count the factor of all sol pnt polys
      for (CFuint jSolPnt = 0; jSolPnt < nbrSolSolDep; ++jSolPnt)
      { 
        const CFuint jSolIdx = solSolDep[iSolPnt*nbrSolSolDep+jSolPnt]; //(*m_solSolDep)[iSolPnt][jSolPnt];

        // Loop over deriv directions and sum them to compute divergence
        for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
        {
          const CFreal polyCoef = solPolyDerivAtSolPnts[iSolPnt*PHYS::DIM*nbSolPnts+iDir*nbSolPnts+jSolIdx];//(*m_solPolyDerivAtSolPnts)[jSolPnt][iDir][iSolIdx]; 
          
          // Loop over conservative fluxes 
          for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
          {
            // Store divFD in the vector that will be divFC
            res[iEq] -= polyCoef*(currFd.getFlux(jSolIdx, iDir)[iEq] - solPntFlx[jSolIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+iDir*SCHEME::MODEL::NBEQS+iEq])*resFactor;

//if (cellID == 11) printf("State: %d, jSol: %d, iDir: %d, var: %d, flx: %f\n",iSolPnt,jSolIdx,iDir,iEq,polyCoef*(solPntFlx[jSolIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+iDir*SCHEME::MODEL::NBEQS+iEq])*resFactor);  
	  }
        }
      }
    }

    // extrapolate the fluxes to the flux points
    for (CFuint iFlxPnt = 0; iFlxPnt < nbrFlxPnts; ++iFlxPnt)
    {
      const CFuint dim = flxPntFlxDim[iFlxPnt];

      // loop over sol pnts to compute flux
      for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
      {
        const CFuint solIdx = flxSolDep[iFlxPnt*nbrFlxSolDep + iSolPnt];

        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
        {
          flxPntFlx[iFlxPnt*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*(currFd.getFlux(solIdx, dim)[iEq] - solPntFlx[solIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+dim*SCHEME::MODEL::NBEQS+iEq]);

          flxPntSol[iFlxPnt*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*states[cell.getStateID(solIdx)*SCHEME::MODEL::NBEQS+iEq];

          for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
          {
            flxPntGrads[iFlxPnt*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*gradients[cell.getStateID(solIdx)*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir];
            flxPntGradsAV[iFlxPnt*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*gradientsAV[cell.getStateID(solIdx)*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir];
          }          
        }
      }
    }

    // set extrapolated states
    for (CFuint iState = 0; iState < nbrFlxPnts; ++iState)
    {
      for (CFuint iEq = 0; iEq < PHYS::NBEQS; ++iEq) 
      {
        currFd.getLstate(iState)[iEq] = flxPntSol[iState*PHYS::NBEQS+iEq];
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
        const CFreal divh = corrFctDiv[iSolPnt*nbrFlxPnts+flxIdx];

        // Fill in the corrections
        for (CFuint iVar = 0; iVar < SCHEME::MODEL::NBEQS; ++iVar)
        {
          res[iVar] += flxPntFlx[flxIdx*SCHEME::MODEL::NBEQS+iVar] * divh * resFactor;
//if (cellID == 11 && iVar == 2) printf("State: %d, flx: %d, var: %d, update: %e, flux: %e, divh: %e\n",iSolPnt,flxIdx,iVar,flxPntFlx[flxIdx*SCHEME::MODEL::NBEQS+iVar] * divh, flxPntFlx[flxIdx*SCHEME::MODEL::NBEQS+iVar], divh);  
        }
      }
    }

    // reset flx pnt fluxes  
    flxPntFlx = 0.0;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSolNeighb;
        
    CudaEnv::CFVec<CFreal,nbFlxPntFlx*PHYS::DIM> flxPntGradNeighb;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx*PHYS::DIM> flxPntGradNeighbAV;
    
    flxPntSolNeighb = 0.0;
        
    flxPntGradNeighb = 0.0;
    
    flxPntGradNeighbAV = 0.0;
    
    for (CFuint iFlxPnt = 0; iFlxPnt < nbFlxPntFlx; ++iFlxPnt) 
    {
        flxPntSolNeighb[iFlxPnt] = 0.0;
    }

    // current neighb cell data
    CellData cells2(nbCells, stateIDs, neighbCellIDs, neighbFaceIDs, innerCellIsLeft, nbrFaces, nbSolPnts, ORDER);

    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      const CFint neighbCellID = cell.getNeighbCellID(iFace);  

      // get current cell
      CellData::Itr cell2 = cells2.getItr(neighbCellID);
      
      // check if it is a bnd face (if so do nothing): for nbd face the neighbCellID will be -1
      if (neighbCellID != -1)
      {
        CFuint jFaceIdx = 0;
          
        for (CFuint jFace = 0; jFace < nbrFaces; ++jFace)
        {
          if (cell2.getNeighbCellID(jFace) == cellID)
          {
            jFaceIdx = jFace; 
            break;
          }
        }

        CFreal waveSpeedUpd = 0.0;
      
        const CFuint faceID = cell.getNeighbFaceID(iFace);

        const bool isLEFT = (bool) cell.getInnerCellIsLeft(iFace);

      // loop over face flx pnts
      for (CFuint iFlxPnt = 0; iFlxPnt < nbrFaceFlxPnts; ++iFlxPnt)
      { 
        // @TODO check if this also works for non QUADs
        const CFuint flxIdx = faceFlxPntConn[iFace*nbrFaceFlxPnts+iFlxPnt];
        const CFuint jFlxIdx = faceFlxPntConn[jFaceIdx*nbrFaceFlxPnts+nbrFaceFlxPnts-1-iFlxPnt];
        
        
        
        
       
    
    // reset the states in the flx pnts
    CFreal epsL = 0.0;
//    CFreal epsR = 0.0;
    
      //m_cellNodes = m_cells[LEFT]->getNodes();

      for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode)
      {
        const CFuint faceNodeIdx = faceNeighbNodeIDs[faceID*nbFaceNodes+iNode];
          
	for (CFuint iNodeCell = 0; iNodeCell < nbrCornerNodes; ++iNodeCell)
        {
          const CFuint nodeIdx = neighbNodeIDs[cellID*nbrCornerNodes+iNodeCell];
            
          //if(cellID == 0) printf("faceID: %d, faceNodeID: %d, cellNodeIdD: %d\n",faceID,faceNodeIdx,nodeIdx);
          
	  if (faceNodeIdx == nodeIdx)
	  {
	    // get node local index
            //const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[LEFT]->getID(),iNodeCell);
	    
            epsL += nodePolyValsAtFlxPnts[flxIdx*nbrCornerNodes+iNodeCell]*nodeEpsilons[nodeIdx]/nbNodeNeighbors[nodeIdx];
            
            //if(cellID == 0) printf("node eps: %e, polyVal: %e, nbNeighb: %d\n",nodeEpsilons[nodeIdx],nodePolyValsAtFlxPnts[flxIdx*nbrCornerNodes+iNodeCell],nbNodeNeighbors[nodeIdx]);
	  }
	}
      }
      
//      m_cellNodes = m_cells[RIGHT]->getNodes();
//
//      for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
//      {
//	for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
//        {
//	  if ((*m_faceNodes)[iNode]->getLocalID() == (*m_cellNodes)[iNodeCell]->getLocalID())
//	  {
//	    // get node local index
//            const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[RIGHT]->getID(),iNodeCell);
//	    
//            epsR += m_nodePolyValsAtFlxPnts[jFlxIdx][iNodeCell]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
//	  }
//	}
//      }
  
  
  
  
  

        // loop over sol pnts to compute sol at flx pnt
        for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
        {
          const CFuint solIdx = flxSolDep[jFlxIdx*nbrFlxSolDep+iSolPnt]; 

          // Loop over conservative vars 
          for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
          {
            flxPntSolNeighb[flxIdx*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[jFlxIdx*nbSolPnts+solIdx]*states[cell2.getStateID(solIdx)*SCHEME::MODEL::NBEQS+iEq];
                        
            for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
            {
              flxPntGradNeighb[flxIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir] += solPolyValsAtFlxPnts[jFlxIdx*nbSolPnts+solIdx]*gradients[cell2.getStateID(solIdx)*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir];
              flxPntGradNeighbAV[flxIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir] += solPolyValsAtFlxPnts[jFlxIdx*nbSolPnts+solIdx]*gradientsAV[cell2.getStateID(solIdx)*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir];
            }
          }
        }

        for (CFuint iEq = 0; iEq < PHYS::NBEQS; ++iEq) 
        {
          currFd.getRstate(flxIdx)[iEq] = flxPntSolNeighb[flxIdx*PHYS::NBEQS+iEq];
          
          avgSol[iEq] = 0.5*(flxPntSolNeighb[flxIdx*PHYS::NBEQS+iEq] + flxPntSol[flxIdx*SCHEME::MODEL::NBEQS+iEq]);
          
          for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
          {
            avgGrad[iEq*PHYS::DIM+iDir] = 0.5*(flxPntGradNeighb[flxIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir] + flxPntGrads[flxIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir]);
            avgGradAV[iEq*PHYS::DIM+iDir] = 0.5*(flxPntGradNeighbAV[flxIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir] + flxPntGradsAV[flxIdx*SCHEME::MODEL::NBEQS*PHYS::DIM+iEq*PHYS::DIM+iDir]);
          }
        } 

        CudaEnv::CFVecSlice<CFreal,PHYS::DIM> n(&(kd.flxPntNormals[faceID*nbrFaceFlxPnts*PHYS::DIM+iFlxPnt*PHYS::DIM]));

        CudaEnv::CFVecSlice<CFreal,PHYS::DIM> nFd(currFd.getFlxScaledNormal(flxIdx));

        CFreal faceVecAbsSize2 = 0.0;
        
        for (CFuint i = 0; i < PHYS::DIM; ++i) 
        {
          nFd[i] = n[i];
          
          faceVecAbsSize2 += n[i]*n[i];
        }

        // get the flux
        fluxScheme.prepareComputation(&currFd, &pmodel);

        fluxScheme(&currFd, &pmodel, iFlxPnt, flxIdx, faceIntCoeff[iFlxPnt], isLEFT, waveSpeedUpd);

        // add diff contribution to wvspd upd
        const CFreal mu = pmodelNS.getUpdateVS()->getDynViscosity(&avgSol[0]);
        const CFreal rho = pmodelNS.getUpdateVS()->getDensity(&avgSol[0]);
        
        const CFreal factorPr = 0.72;//min(pmodelNS.getUpdateVS()->getModel().getPrandtl(),1.0);
        
        if (addUpdCoeff)
        {
          waveSpeedUpd += (mu/rho/factorPr+epsL)*faceVecAbsSize2*faceIntCoeff[iFlxPnt]/currVol*cflConvDiffRatio;
        }
        else
        {
          waveSpeedUpd += mu/rho/factorPr*faceVecAbsSize2*faceIntCoeff[iFlxPnt]/currVol*cflConvDiffRatio; 
        }

        pmodelNS.getUpdateVS()->getFlux(&avgSol[0],&avgGrad[0],&n[0],&currFlxPntFlx[0]);
        
        // compute artificial part
        // get epsilon
        const CFreal epsilon = epsL;//0.5*(epsL+epsR);
    
        for (CFuint iDim = 0; iDim < PHYS::DIM; ++iDim)
        {
          for (CFuint iVar = 0; iVar < SCHEME::MODEL::NBEQS; ++iVar)
          {
            currFlxPntFlx[iVar] += epsilon*avgGradAV[iVar*PHYS::DIM+iDim]*n[iDim];
            
            //if(cellID == 0) printf("second flx: %d, var: %d, dim: %d, eps: %e, grad: %e, n: %e\n",iFlxPnt,iVar,iDim,epsilon,avgGradAV[iVar*PHYS::DIM+iDim],n[iDim]*faceDir[cellID*totNbrFlxPnts+flxIdx]); 
          }
        }
        
        // extrapolate the fluxes to the flux points
        for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
        {     
          const CFuint solIdx = flxSolDep[flxIdx*nbrFlxSolDep+iSolPnt];

          // get current state ID
          const CFuint stateID = cell.getStateID(solIdx);

          // get current vector slice out of rhs
          CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> res(&rhs[stateID*SCHEME::MODEL::NBEQS]);   

          // divergence of the correction function
          const CFint currFaceDir = faceDir[cellID*totNbrFlxPnts+flxIdx];
          const CFreal divh = corrFctDiv[solIdx*nbrFlxPnts+flxIdx] * currFaceDir;
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < SCHEME::MODEL::NBEQS; ++iVar)
          {
            res[iVar] -= (currFd.getInterfaceFlux(flxIdx)[iVar] - currFlxPntFlx[iVar]) * divh * resFactor;
//if(cellID == 11 && flxIdx == 1 && iVar == 2) printf("State: %d, flx: %d, var: %d, divh: %e. up: %e\n",solIdx,flxIdx,iVar,divh,currFlxPntFlx[iVar] * divh * resFactor); 
          }
        }
      }

      //CFreal* waveSpeedUpd = currFd.getUpdateCoeff();

      for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
      {  
        // get current state ID
        const CFuint stateID = cell.getStateID(iSolPnt);

        updateCoeff[stateID] += waveSpeedUpd*(2.0*ORDER+1);
      }
 
      //currFd.resetUpdateCoeff();
      }
    }
    
//    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
//      {  
//        // get current state ID
//        const CFuint stateID = cell.getStateID(iSolPnt);
//
//        printf("cellID: %d, stateID: %d, resV: %e\n", cellID, stateID, rhs[stateID*SCHEME::MODEL::NBEQS+2]);
//      }
  }
}
  
//////////////////////////////////////////////////////////////////////////////

template <typename SCHEME, typename PHYS, typename PHYSNS, CFuint ORDER>
__global__ void computeGradientsKernel(typename SCHEME::MODEL::PTERM::template DeviceConfigOptions<NOTYPE>* dcop,
                                       typename PHYSNS::DTERM::template DeviceConfigOptions<NOTYPE>* dcopNS,
                                       typename PHYSNS::PTERM::template DeviceConfigOptions<NOTYPE>* dcopNSConv,
                                       const CFuint nbCells,
				       CFreal* states, 
                                       CFreal* gradients,
                                       CFreal* gradientsAV,
                                       CFreal* solPntNormals,
                                       CFreal* flxPntNormals,
                                       CFreal* cellVolumes,
                                       CFreal* volumes,
                                       CFint* faceDir,
                                       CFreal* nodeEpsilons,
                                       CFreal* cellEpsilons,
                                       const CFuint nbSolPnts,
                                       const CFuint nbrFaces,
                                       const CFuint* faceFlxPntConn,
                                       const CFuint* stateIDs,
                                       const CFint* neighbCellIDs,
                                       const CFuint* neighbFaceIDs,
                                       const CFuint* neighbNodeIDs,
                                       const CFuint* innerCellIsLeft,
                                       const CFuint nbrFlxPnts,
                                       const CFuint nbrSolSolDep,
                                       const CFuint* solSolDep,
                                       const CFuint nbrSolFlxDep,
                                       const CFuint* solFlxDep,
                                       const CFuint nbrFlxSolDep,
                                       const CFuint* flxSolDep,
                                       const CFreal* solPolyDerivAtSolPnts,
                                       const CFreal* solPolyValsAtFlxPnts,
                                       const CFuint* flxPntFlxDim,
                                       const CFreal* corrFctDiv,
                                       const CFreal* transformationMatrix,
                                       const CFreal peclet,
                                       const CFreal subcellRes,
                                       const CFreal kappa,
                                       const CFreal s0,
                                       const CFuint monitoredVar,
                                       const CFreal monitoredPhysVar,
                                       const CFuint nbrCornerNodes,
                                       const bool useMax,
                                       const bool flagComputeNbNghb)
{    
  // one thread per cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (cellID < nbCells) 
  { 
    // current kernel data
    //KernelData<CFreal> kd (nbCells, states, updateCoeff, rhs, solPntNormals, flxPntNormals, faceDir, nbSolPnts);

    // physical model
    typename SCHEME::MODEL pmodel(dcop);
    //typename PHYSNS pmodelNS(dcopNS);
    PHYSNS pmodelNS(dcopNS,dcopNSConv);    

    // current cell data
    CellData cells(nbCells, stateIDs, neighbCellIDs, neighbFaceIDs, innerCellIsLeft, nbrFaces, nbSolPnts, ORDER);
    
    // get current cell
    CellData::Itr cell = cells.getItr(cellID);
          
    const CFuint nbFlxPntFlx = SCHEME::MODEL::NBEQS*(ORDER+1)*2*PHYS::DIM;//8;
    
    //const CFuint nbFaceFlxPntFlx = SCHEME::MODEL::NBEQS*(ORDER+1);
   
    const CFuint nbrFaceFlxPnts = (ORDER+1);

    const CFuint totNbrFlxPnts = (ORDER+1)*2*PHYS::DIM;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntFlx;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSol;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSolAV;
    
    CudaEnv::CFVec<CFreal,SCHEME::MODEL::NBEQS*(ORDER+1)*(ORDER+1)> projStates;

    CudaEnv::CFVec<CFreal,SCHEME::MODEL::NBEQS> stateGradVars;
    //typename MathTypes<CFreal, GPU, SCHEME::MODEL::NBEQS>::VEC stateGradVars;
    
    CudaEnv::CFVec<CFreal,SCHEME::MODEL::DATASIZE> pdata;
    
    flxPntFlx = 0.0;
    
    flxPntSol = 0.0;
    
    flxPntSolAV = 0.0;

    stateGradVars = 0.0;
    
    projStates = 0.0;  
    
    pdata = 0.0;
    
    CFreal currVol = cellVolumes[cellID];
   
    ////////COMPUTE PROJECTED STATES//////////////////////////////////////////////////////////////////////
    if (ORDER != 1)
    {
      for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
      {
        for (CFuint iSol = 0; iSol < nbSolPnts; ++iSol)
        {
          // get current state ID
          const CFuint stateID = cell.getStateID(iSol);
        
          CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[stateID*SCHEME::MODEL::NBEQS]);
          
          for (CFuint jSol = 0; jSol < nbSolPnts; ++jSol)
          {
            projStates[jSol*SCHEME::MODEL::NBEQS+iEq] += currState[iEq]*transformationMatrix[jSol*nbSolPnts+iSol];
          }
        }
      }
    }
  else
  {
    for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
    {
      CFreal stateSum = 0.0;
      
      for (CFuint iSol = 0; iSol < nbSolPnts; ++iSol)
      {
        // get current state ID
        const CFuint stateID = cell.getStateID(iSol);
        
        CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[stateID*SCHEME::MODEL::NBEQS]);
        
        stateSum += currState[iEq];
      }

      stateSum /= nbSolPnts;

      for (CFuint iSol = 0; iSol < nbSolPnts; ++iSol)
      {
        projStates[iSol*SCHEME::MODEL::NBEQS+iEq] = stateSum;
        
        //if (cellID == 0) printf("projState %d, %d: %e\n",iSol, iEq, stateSum);
      }
    }
  }
    
    
    
  ////////COMPUTE EPS0//////////////////////////////////////////////////////////////////
    
  // compute a cell average characteristic flow speed. Note that a straight average is used, not a weighted one, maybe change this
  CFreal wavespeed = 0.0;

  for (CFuint iSol = 0; iSol < nbSolPnts; ++iSol)
  {    
    // get current state ID
    const CFuint stateID = cell.getStateID(iSol);
        
    CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[stateID*SCHEME::MODEL::NBEQS]);
    
    pmodel.getUpdateVS()->computePhysicalData(&currState[0], &pdata[0]);

    wavespeed += pdata[EulerTerm::V] + pdata[EulerTerm::A];
  }
  
  wavespeed /= nbSolPnts;
      
  const CFreal oneOverDim = 1./PHYS::DIM;
  
  const CFreal h = pow(currVol,oneOverDim);

  const CFreal eps0 = max(h*wavespeed*(2.0/peclet - subcellRes/peclet),0.0);
  
  //if (cellID == 0) printf("eps0 %e\n",eps0);
  
  
  
  ////////COMPUTE SMOOTHNESS///////////////////////////////////////////////////////////
  
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;

  if (monitoredPhysVar < SCHEME::MODEL::DATASIZE)
  {
    for (CFuint iSol = 0; iSol < nbSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;
      
      // get current state ID
      const CFuint stateID = cell.getStateID(iSol);
        
      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[stateID*SCHEME::MODEL::NBEQS]);
      
      pmodel.getUpdateVS()->computePhysicalData(&currState[0], &pdata[0]);
      
      stateP = pdata[monitoredPhysVar];
      
      pmodel.getUpdateVS()->computePhysicalData(&projStates[iSol*SCHEME::MODEL::NBEQS], &pdata[0]);

      diffStatesPPMinOne = stateP - pdata[monitoredPhysVar];

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < nbSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;
      
      // get current state ID
      const CFuint stateID = cell.getStateID(iSol);
        
      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[stateID*SCHEME::MODEL::NBEQS]);

      stateP = currState[monitoredVar];
      diffStatesPPMinOne = stateP - projStates[iSol*SCHEME::MODEL::NBEQS+monitoredVar];

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  
  CFreal smoothness = 0.0;
      
  if (sNum <= 1.0e-10 || sDenom <= 1.0e-10)
  {
    smoothness = -100.0;
  }
  else
  {
    smoothness = log10(sNum/sDenom);
  }
  
  //if (cellID == 0) printf("s %e\n",smoothness);
  
  ///////COMPUTE EPS///////////////////////////////////////////////////////////////////
  
  CFreal eps = 0.0;
  
  if (smoothness > s0 + kappa)
  {
    eps = eps0;
  }
  else if (smoothness > s0 - kappa)
  {
    eps = eps0*0.5*(1.0 + sin(0.5*3.141592653589793238462643383*(smoothness-s0)/kappa));
  }
  
  //if (cellID == 0) printf("eps %e, s %e, s0 %e, k %e\n",eps,smoothness,s0,kappa);

//  if (m_useWallCutOff)
//  {
//    // Get the wall distance
//    DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
//  
//    CFreal centroidDistance = 0.0;
//      
//    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//    {
//      const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
//      centroidDistance += wallDist[stateID];
//    }
//    
//    centroidDistance /= m_nbrSolPnts;
//    
//    if (centroidDistance < m_wallCutOff) 
//    {
//      if (centroidDistance < 0.5*m_wallCutOff)
//      {
//        m_epsilon = 0.0; 
//      }
//      else
//      {
//        m_epsilon *= 0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(centroidDistance-0.75*m_wallCutOff)/(0.25*m_wallCutOff)));
//      }
//    }
//  }
  
  if (eps < 0.0 || eps != eps) 
  {
    eps = 0.0;
  }
  
  
  ////////STORE EPS///////////////////////////////////////////////////////////////////
  
  for (CFuint iNode = 0; iNode < nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = neighbNodeIDs[cellID*nbrCornerNodes+iNode];

    if (!useMax) 
    {
      //nodeEpsilons[nodeID] += eps;
      atomicAdd(&nodeEpsilons[nodeID],eps);
      cellEpsilons[cellID] = eps;
    }
    else
    {
      const CFreal maxEps = max(eps, cellEpsilons[cellID]);
      //nodeEpsilons[nodeID] += maxEps;
      atomicAdd(&nodeEpsilons[nodeID],maxEps);
      cellEpsilons[cellID] = maxEps;
    }
  }


    // loop over sol pnts to compute flux
    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);

      //typename MathTypes<CFreal, GPU, SCHEME::MODEL::NBEQS>::SLICEVEC currState(&states[stateID*SCHEME::MODEL::NBEQS]);

      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[stateID*SCHEME::MODEL::NBEQS]);

      pmodelNS.getUpdateVS()->setGradientVars(&currState[0],&stateGradVars[0]);

      CudaEnv::CFVecSlice<CFreal,PHYS::DIM*PHYS::DIM> currNormals(&solPntNormals[stateID*PHYS::DIM*PHYS::DIM]);

      // Loop over solution pnts to count the factor of all sol pnt polys
      for (CFuint jSolPnt = 0; jSolPnt < nbrSolSolDep; ++jSolPnt)
      { 
        const CFuint jSolIdx = solSolDep[iSolPnt*nbrSolSolDep+jSolPnt];
        
        // get current j state ID
        const CFuint jStateID = cell.getStateID(jSolIdx);
        
        // get current vector slice out of gradients
        CudaEnv::CFVecSlice<CFreal,PHYS::DIM*SCHEME::MODEL::NBEQS> grad(&gradients[jStateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);
        CudaEnv::CFVecSlice<CFreal,PHYS::DIM*SCHEME::MODEL::NBEQS> gradAV(&gradientsAV[jStateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);

        // Loop over deriv directions and sum them to compute divergence
        for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
        {
          const CFreal polyCoef = solPolyDerivAtSolPnts[jSolIdx*PHYS::DIM*nbSolPnts+iDir*nbSolPnts+iSolPnt];
          
          for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
          {
            for (CFuint jDir = 0; jDir < PHYS::DIM; ++jDir)
            {            
              // Store divFD in the vector that will be divFC
              grad[iEq*PHYS::DIM+jDir] += polyCoef*currNormals[iDir*PHYS::DIM+jDir]*stateGradVars[iEq];//*states[stateID*SCHEME::MODEL::NBEQS+iEq]; 
              gradAV[iEq*PHYS::DIM+jDir] += polyCoef*currNormals[iDir*PHYS::DIM+jDir]*currState[iEq];//*states[stateID*SCHEME::MODEL::NBEQS+iEq]; 

              //if (cellID == 11) printf("after  iSol: %d, iEq: %d, iDir: %d: %f\n", iSolPnt, iEq, jDir, grad[iEq*PHYS::DIM+jDir]); 
	    }
          }
        }
      }
    }

    // extrapolate the fluxes to the flux points
    for (CFuint iFlxPnt = 0; iFlxPnt < nbrFlxPnts; ++iFlxPnt)
    {
      //const CFuint dim = flxPntFlxDim[iFlxPnt];

      // loop over sol pnts to compute flux
      for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
      {
        const CFuint solIdx = flxSolDep[iFlxPnt*nbrFlxSolDep + iSolPnt];

        CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[cell.getStateID(solIdx)*SCHEME::MODEL::NBEQS]);

        pmodelNS.getUpdateVS()->setGradientVars(&currState[0],&stateGradVars[0]);

        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
        {
          flxPntSol[iFlxPnt*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*stateGradVars[iEq];    
          flxPntSolAV[iFlxPnt*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*currState[iEq]; 
        }
      }
    }

    // reset flx pnt fluxes  
    flxPntFlx = 0.0;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSolNeighb;
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSolNeighbAV;
    
    for (CFuint iFlxPnt = 0; iFlxPnt < nbFlxPntFlx; ++iFlxPnt) 
    {
        flxPntSolNeighb[iFlxPnt] = 0.0;
        flxPntSolNeighbAV[iFlxPnt] = 0.0;
    }

    // current neighb cell data
    CellData cells2(nbCells, stateIDs, neighbCellIDs, neighbFaceIDs, innerCellIsLeft, nbrFaces, nbSolPnts, ORDER);

    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      const CFint neighbCellID = cell.getNeighbCellID(iFace);  

      // get current cell
      CellData::Itr cell2 = cells2.getItr(neighbCellID);
      
      if (neighbCellID != -1)
      {
        CFuint jFaceIdx = 0;

        for (CFuint jFace = 0; jFace < nbrFaces; ++jFace)
        {
          if (cell2.getNeighbCellID(jFace) == cellID)
          {
            jFaceIdx = jFace; 
            break;
          }
        }

        const CFuint faceID = cell.getNeighbFaceID(iFace);

        const bool isLEFT = (bool) cell.getInnerCellIsLeft(iFace);

        // loop over face flx pnts
        for (CFuint iFlxPnt = 0; iFlxPnt < nbrFaceFlxPnts; ++iFlxPnt)
        { 
          // @TODO check if this also works for non QUADs
          const CFuint flxIdx = faceFlxPntConn[iFace*nbrFaceFlxPnts+iFlxPnt];
          const CFuint jFlxIdx = faceFlxPntConn[jFaceIdx*nbrFaceFlxPnts+nbrFaceFlxPnts-1-iFlxPnt];

          const CFreal dirFactor = faceDir[cellID*totNbrFlxPnts+flxIdx];

          // loop over sol pnts to compute sol at flx pnt
          for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
          {
            const CFuint solIdx = flxSolDep[jFlxIdx*nbrFlxSolDep+iSolPnt]; 

            CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[cell2.getStateID(solIdx)*SCHEME::MODEL::NBEQS]);

            pmodelNS.getUpdateVS()->setGradientVars(&currState[0],&stateGradVars[0]);

            // Loop over conservative vars 
            for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
            {
              flxPntSolNeighb[flxIdx*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[jFlxIdx*nbSolPnts+solIdx]*stateGradVars[iEq];
              flxPntSolNeighbAV[flxIdx*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[jFlxIdx*nbSolPnts+solIdx]*currState[iEq];
            }
          }

          // extrapolate the fluxes to the flux points
          for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
          {     
            const CFuint solIdx = flxSolDep[flxIdx*nbrFlxSolDep+iSolPnt];

            // get current state ID
            const CFuint stateID = cell.getStateID(solIdx); 

            // get current vector slice out of gradients
            CudaEnv::CFVecSlice<CFreal,PHYS::DIM*SCHEME::MODEL::NBEQS> grad(&gradients[stateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);
            CudaEnv::CFVecSlice<CFreal,PHYS::DIM*SCHEME::MODEL::NBEQS> gradAV(&gradientsAV[stateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);

            // divergence of the correction function
            const CFreal divh = corrFctDiv[solIdx*nbrFlxPnts+flxIdx];
           
            for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
            {
              const CFreal corrFactor = 0.5*(flxPntSolNeighb[flxIdx*SCHEME::MODEL::NBEQS+iEq]-flxPntSol[flxIdx*SCHEME::MODEL::NBEQS+iEq]);
              const CFreal corrFactorAV = 0.5*(flxPntSolNeighbAV[flxIdx*SCHEME::MODEL::NBEQS+iEq]-flxPntSolAV[flxIdx*SCHEME::MODEL::NBEQS+iEq]);

              // Loop over deriv directions and sum them to compute divergence
              for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
              {
                grad[iEq*PHYS::DIM+iDir] += divh*corrFactor*flxPntNormals[faceID*nbrFaceFlxPnts*PHYS::DIM+iFlxPnt*PHYS::DIM+iDir]*dirFactor; 
                gradAV[iEq*PHYS::DIM+iDir] += divh*corrFactorAV*flxPntNormals[faceID*nbrFaceFlxPnts*PHYS::DIM+iFlxPnt*PHYS::DIM+iDir]*dirFactor; 

//              if (cellID == 11) printf("iSol: %d, iEq: %d, iFlx %d, iDir: %d: %e\n", solIdx, iEq, flxIdx, iDir,
//                      divh*corrFactor*flxPntNormals[faceID*nbrFaceFlxPnts*PHYS::DIM+iFlxPnt*PHYS::DIM+iDir]*dirFactor);  
	      }
            }
          }
        }
      }
    }

    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {  
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);

      // get current vector slice out of gradients
      CudaEnv::CFVecSlice<CFreal,PHYS::DIM*SCHEME::MODEL::NBEQS> grad(&gradients[stateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);
      CudaEnv::CFVecSlice<CFreal,PHYS::DIM*SCHEME::MODEL::NBEQS> gradAV(&gradientsAV[stateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);

      const CFreal invacob = 1.0/volumes[stateID];

      for (CFuint jDir = 0; jDir < PHYS::DIM; ++jDir)
      {
        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
        {
          grad[iEq*PHYS::DIM+jDir] *= invacob;
          gradAV[iEq*PHYS::DIM+jDir] *= invacob;
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

template <typename SCHEME, typename PHYSICS, typename PHYSICSNS, CFuint ORDER, CFuint NB_BLOCK_THREADS>
void ConvDiffLLAVRHSFluxReconstructionCUDA<SCHEME,PHYSICS,PHYSICSNS,ORDER,NB_BLOCK_THREADS>::execute()
{
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  CFTRACEBEGIN;
  
  CFLog(VERBOSE, "ConvDiffLLAVRHSFluxReconstructionCUDA::execute() START\n");
  
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
  
  //initializeComputationRHS();

  const CFuint nbStates = socket_states.getDataHandle().size();
  cf_assert(nbStates > 0);

  CFLog(VERBOSE, "nbCells: " << nbCells << ", nbStates: " << nbStates << "\n");

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle(); 
  DataHandle<CFreal> solPntNormals = socket_solPntNormals.getDataHandle(); 
  DataHandle<CFreal> flxPntNormals = socket_flxPntNormals.getDataHandle(); 
  DataHandle<CFint> faceDir = socket_faceDir.getDataHandle(); 
  DataHandle<CFreal> gradients = socket_gradientsCUDA.getDataHandle();
  DataHandle<CFreal> gradientsAV = socket_gradientsAVCUDA.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> cellVolumes = socket_cellVolumes.getDataHandle();
 

  SafePtr<SCHEME> lf  = getMethodData().getRiemannFlux().d_castTo<SCHEME>();
  SafePtr<typename PHYSICS::PTERM> phys = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<typename PHYSICS::PTERM>();

  SafePtr<typename PHYSICSNS::DTERM> physNS = PhysicalModelStack::getActive()->getImplementor()->
    getDiffusiveTerm().d_castTo<typename PHYSICSNS::DTERM>();

  SafePtr<typename PHYSICSNS::PTERM> physNSConv = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<typename PHYSICSNS::PTERM>();
  
#ifdef CF_HAVE_CUDA
  typedef typename SCHEME::template DeviceFunc<GPU, PHYSICS, ORDER> FluxScheme;  
#else
  typedef typename SCHEME::template DeviceFunc<CPU, PHYSICS, ORDER> FluxScheme;
#endif 
  
  // get current iteration
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  // check if LLAV should be frozen
  m_useMax = iter > m_freezeLimiterIter;
  
  if (m_onGPU) 
  {
#ifdef CF_HAVE_CUDA

    //CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
    //timer.start();
      
    // copy of data that change at every iteration
    socket_states.getDataHandle().getGlobalArray()->put(); 
    socket_gradientsCUDA.getDataHandle().getLocalArray()->put();
    socket_gradientsAVCUDA.getDataHandle().getLocalArray()->put();
    socket_rhs.getDataHandle().getLocalArray()->put(); 
    socket_updateCoeff.getDataHandle().getLocalArray()->put();
    
    //CFLog(VERBOSE, "nb normals: " << socket_solPntNormals.getDataHandle().size() << ", n0: " << socket_solPntNormals.getDataHandle()[0] << "\n");

    socket_faceDir.getDataHandle().getLocalArray()->put();
    socket_solPntNormals.getDataHandle().getLocalArray()->put();
    socket_flxPntNormals.getDataHandle().getLocalArray()->put();
    socket_volumes.getDataHandle().getLocalArray()->put();
    socket_cellVolumes.getDataHandle().getLocalArray()->put();

    DataHandle<Framework::State*, Framework::GLOBAL > statesI = socket_states.getDataHandle();
     
    //CFLog(VERBOSE, "ConvDiffLLAVRHSFluxReconstructionCUDA::execute() => CPU-->GPU data transfer took " << timer.elapsed() << " s\n");
    //timer.start();
    
    ConfigOptionPtr<SCHEME,  NOTYPE, GPU> dcof(lf);
    ConfigOptionPtr<typename PHYSICS::PTERM, NOTYPE, GPU> dcop(phys);
    ConfigOptionPtr<typename PHYSICSNS::DTERM, NOTYPE, GPU> dcopNS(physNS);
    ConfigOptionPtr<typename PHYSICSNS::PTERM, NOTYPE, GPU> dcopNSConv(physNSConv);

    const CFuint nThreads = m_nbCellsPerBlock;//512; //CudaEnv::CudaDeviceManager::getInstance().getNThreads();
    const CFuint blocksPerGrid = ceil(nbCells*1.0/nThreads); //CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
    
    CFLog(VERBOSE, "blocksPerGrid: " << blocksPerGrid << ", threads: " << nThreads << "\n");

    // boolean telling whether there is a diffusive term
    const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();

//CFuint megabytesToUse = 24;
//size_t newHeapSize = 1024 * 1000 * megabytesToUse;
//cudaDeviceSetLimit(cudaLimitMallocHeapSize, newHeapSize);
//printf("Adjusted heap size to be %d\n",(int) newHeapSize);

    //dim3 blocks(m_nbBlocksPerGridX, m_nbBlocksPerGridY);
    
    //cudaFuncSetCacheConfig("computeGradientsKernel", cudaFuncCachePreferL1);

    // get residual factor
    const CFreal resFactor = getMethodData().getResFactor();
    
    m_nodeEpsilons = 0.0;
    
    m_nodeEpsilons.put();
    m_cellEpsilons.put();
    
    // cudaFuncSetCacheConfig("computeFluxKernel", cudaFuncCachePreferL1);

    // if there is a diffusive term, compute the gradients
//    if (hasDiffTerm)
//    {
        CFLog(VERBOSE, "grad kernel\n");
      computeGradientsKernel<FluxScheme,PHYSICS,PHYSICSNS,ORDER> <<<blocksPerGrid,nThreads>>>(
                                       dcop.getPtr(),
                                       dcopNS.getPtr(),
                                       dcopNSConv.getPtr(),
                                       nbCells,
				       socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
                                       gradients.getLocalArray()->ptrDev(), 
                                       gradientsAV.getLocalArray()->ptrDev(), 
                                       solPntNormals.getLocalArray()->ptrDev(),
                                       flxPntNormals.getLocalArray()->ptrDev(),
                                       cellVolumes.getLocalArray()->ptrDev(),
                                       volumes.getLocalArray()->ptrDev(),
                                       faceDir.getLocalArray()->ptrDev(),
                                       m_nodeEpsilons.ptrDev(),
                                       m_cellEpsilons.ptrDev(),
                                       m_nbrSolPnts,
                                       4,
                                       m_faceFlxPntConn2.ptrDev(),
                                       m_stateIDs.ptrDev(),
                                       m_neighbCellIDs.ptrDev(),
                                       m_neighbFaceIDs.ptrDev(),
                                       m_neighbNodeIDs.ptrDev(),
                                       m_innerCellIsLeft.ptrDev(),
                                       m_nbrFlxPnts,
                                       m_nbrSolSolDep,
                                       m_solSolDep2.ptrDev(),
                                       m_nbrFlxDep,
                                       m_solFlxDep2.ptrDev(),
                                       m_nbrSolDep,
                                       m_flxSolDep2.ptrDev(),
                                       m_solPolyDerivAtSolPnts2.ptrDev(),
                                       m_solPolyValsAtFlxPnts2.ptrDev(),
                                       m_flxPntFlxDim2.ptrDev(),
                                       m_corrFctDiv2.ptrDev(),
                                       m_transformationMatrix2.ptrDev(),
                                       m_peclet,
                                       m_subcellRes,
                                       m_kappa,
                                       m_s0,
                                       m_monitoredVar,
                                       m_monitoredPhysVar,
                                       m_nbrCornerNodes,
                                       m_useMax,
                                       m_flagComputeNbNghb);
//    }
        CFLog(VERBOSE, "sol kernel\n");

    m_flagComputeNbNghb = false;
        
    // compute the convective flux in each cell
    computeStateLocalRHSKernel<FluxScheme,PHYSICS,PHYSICSNS,ORDER> <<<blocksPerGrid,nThreads>>> 
      (dcof.getPtr(),
       dcop.getPtr(),
       dcopNS.getPtr(),
       dcopNSConv.getPtr(),
       nbCells,
       resFactor,
       socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
       gradients.getLocalArray()->ptrDev(),
       gradientsAV.getLocalArray()->ptrDev(), 
       updateCoeff.getLocalArray()->ptrDev(), 
       rhs.getLocalArray()->ptrDev(),
       solPntNormals.getLocalArray()->ptrDev(),
       flxPntNormals.getLocalArray()->ptrDev(),
       cellVolumes.getLocalArray()->ptrDev(),
       faceDir.getLocalArray()->ptrDev(),
       m_nbrSolPnts,
       4,
       m_faceFlxPntConn2.ptrDev(),
       m_stateIDs.ptrDev(),
       m_neighbCellIDs.ptrDev(),
       m_neighbFaceIDs.ptrDev(),
       m_innerCellIsLeft.ptrDev(),
       m_nbrFlxPnts,
       m_nbrSolSolDep,
       m_solSolDep2.ptrDev(),
       m_nbrFlxDep,
       m_solFlxDep2.ptrDev(),
       m_nbrSolDep,
       m_flxSolDep2.ptrDev(),
       m_solPolyDerivAtSolPnts2.ptrDev(),
       m_solPolyValsAtFlxPnts2.ptrDev(),
       m_flxPntFlxDim2.ptrDev(),
       m_corrFctDiv2.ptrDev(),
       m_faceIntegrationCoefs2.ptrDev(),
       m_cflConvDiffRatio,
       m_nbNodeNeighbors.ptrDev(),
       m_nodeEpsilons.ptrDev(),
       m_nbrCornerNodes,
       m_neighbNodeIDs.ptrDev(),
       m_faceNeighbNodeIDs.ptrDev(),
       m_nbFaceNodes,
       m_nodePolyValsAtFlxPnts2.ptrDev(),
       m_nodePolyValsAtSolPnts2.ptrDev(),
       m_addUpdCoeff);
    
   
    cudaDeviceSynchronize();
    
    
    //for (CFuint i = 0; i < m_solPolyValsAtFlxPnts2.size(); ++i) {CFLog(INFO, "thing: " << m_solPolyValsAtFlxPnts2[i] << "\n");}
    
    //CFLog(INFO, "After Kernel, size: " << socket_states.getDataHandle().size() << "\n");
    
    //CFLog(VERBOSE, "ConvDiffLLAVRHSFluxReconstructionCUDA::execute() => computeFluxKernel took " << timer.elapsed() << " s\n");
    
    //for (CFuint i = 0; i < rhs.size(); ++i) {CFLog(INFO, "res before: " << rhs[i] << "\n");}
    
    //RealVector rhsB;
    //rhsB.resize(rhs.size());
    //for (CFuint i = 0; i < rhs.size(); ++i) {rhsB[i] = rhs[i];}
    
    //timer.start();
    rhs.getLocalArray()->get();
    updateCoeff.getLocalArray()->get();
    gradients.getLocalArray()->get();
    gradientsAV.getLocalArray()->get();
    
    //m_nodeEpsilons.get();
    m_cellEpsilons.get();
    
    //for (CFuint i = 0; i < rhs.size(); ++i) {CFLog(INFO, "res after: " << rhs[i]-rhsB[i] << "\n");}
    //CFLog(VERBOSE, "ConvDiffLLAVRHSFluxReconstructionCUDA::execute() => GPU-->CPU data transfer took " << timer.elapsed() << " s\n");
    //CFLog(INFO, "resSize: " << rhs.size() << "\n");
    //for (CFuint i = 0; i < rhs.size(); ++i)
    //{
      //if (abs(rhs[i]) > 1.0e-10) CFLog(INFO, "res " << i << ": " << rhs[i] << "\n");
    //}

  #endif
  }
  else 
  {
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();

  // loop over element types, for the moment there should only be one
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
      // if the states in the cell are parallel updatable or the gradients need to be computed, set the cell data
      if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
      {
	// set the cell data
	setCellData();
      }
      
      // if the states in the cell are parallel updatable, compute the divergence of the discontinuous flx (-divFD+divhFD)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
	// compute the divergence of the discontinuous flux (-divFD+divhFD)
	computeDivDiscontFlx(m_divContFlx);
      
	// update RHS
        updateRHS();
      } 
      
      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
	computeGradients();
      }
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 35) //true) //
      {
	CFLog(VERBOSE, "ID  = " << (*m_cellStates)[0]->getLocalID() << "\n");
        CFLog(VERBOSE, "coords  = " << (*m_cellStates)[0]->getCoordinates() << "\n");
        CFLog(VERBOSE, "UpdateTotal = \n");
        // get the datahandle of the rhs
        DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            CFLog(VERBOSE, "" << rhs[resID+iVar] << " ");
          }
          CFLog(VERBOSE,"\n");
          DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
          CFLog(VERBOSE, "UpdateCoeff: " << updateCoeff[(*m_cellStates)[iState]->getLocalID()] << "\n");
	  CFLog(VERBOSE, "state " << iState << ": " << *(((*m_cellStates)[iState])->getData()) << "\n");
        }
      }
      
      if(m_cell->getID() == 35 && hasDiffTerm)
      {
	// get the gradients
        DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
	  CFuint solID = ((*m_cellStates)[iState])->getLocalID();
          for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
          {
	    CFLog(VERBOSE, "total gradient " << iGrad << " of  " << iState << ": " << gradients[solID][iGrad] << "\n");
          } 
        }
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
	  CFLog(VERBOSE, "state " << iState << ": " << *(((*m_cellStates)[iState])->getData()) << "\n");
	}
      }
      
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  //// Loop over faces to calculate fluxes and interface fluxes in the flux points
  
  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    CFLog(VERBOSE, "Orient = " << m_orient << "\n");
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      // if one of the neighbouring cells is parallel updatable or if the gradients have to be computed, set the bnd face data
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable() || hasDiffTerm)
      {
	// set the bnd face data
        setFaceData(m_face->getID());//faceID

	// compute the states in the flx pnts
        computeFlxPntStates();

	// compute the interface flux
	computeInterfaceFlxCorrection();
          
	// compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);

        // update the wave speed
        updateWaveSpeed();
      }
	
	// if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
	// compute the correction for the left neighbour
	computeCorrection(LEFT, m_divContFlxL);
	
	// compute the correction for the right neighbour
	computeCorrection(RIGHT, m_divContFlxR);
	
	// update RHS
	updateRHSBothSides();
      }
      
      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
	// compute the face correction term of the corrected gradients
        computeGradientFaceCorrections();
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }

    //DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

    //for (CFuint i = 0; i < rhs.size(); ++i)
    //{
      //if (abs(rhs[i]) > 1.0e-10) CFLog(INFO, "res " << i << ": " << rhs[i] << "\n");
    //}
  }
  
  //finalizeComputationRHS();
  
  CFLog(VERBOSE, "ConvDiffLLAVRHSFluxReconstructionCUDA::execute() END\n");
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
