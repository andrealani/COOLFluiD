#include "FluxReconstructionCUDA/ConvDiffRHSFluxReconstructionCUDA.hh"
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

#define FR_NS_RHS_PROV(__dim__,__svars__,__uvars__,__order__,__nbBThreads__,__providerName__) \
MethodCommandProvider<ConvDiffRHSFluxReconstructionCUDA<LaxFriedrichsFlux, \
                      VarSetListT<Euler##__dim__##__svars__##T, Euler##__dim__##__uvars__##T>,NSVarSetListT<NavierStokes##__dim__##__svars__##T, NavierStokes##__dim__##__uvars__##T>,__order__,__nbBThreads__>, \
		      FluxReconstructionSolverData,FluxReconstructionCUDAModule>	\
FR_RhsNS##__dim__##__svars__##__uvars__##__order__##__nbBThreads__##Provider(__providerName__);
// 48 block threads (default)
FR_NS_RHS_PROV(2D, Cons, Cons, 0, 48, "NSFRLaxFriedrichs2DConsP0")
FR_NS_RHS_PROV(2D, Cons, Cons, 1, 48, "NSFRLaxFriedrichs2DConsP1")
FR_NS_RHS_PROV(2D, Cons, Cons, 2, 48, "NSFRLaxFriedrichs2DConsP2")
FR_NS_RHS_PROV(2D, Cons, Cons, 3, 48, "NSFRLaxFriedrichs2DConsP3")
FR_NS_RHS_PROV(2D, Cons, Cons, 4, 48, "NSFRLaxFriedrichs2DConsP4")
FR_NS_RHS_PROV(2D, Cons, Cons, 5, 48, "NSFRLaxFriedrichs2DConsP5")
FR_NS_RHS_PROV(2D, Cons, Cons, 6, 48, "NSFRLaxFriedrichs2DConsP6")
FR_NS_RHS_PROV(2D, Cons, Cons, 7, 48, "NSFRLaxFriedrichs2DConsP7")
FR_NS_RHS_PROV(2D, Cons, Cons, 8, 48, "NSFRLaxFriedrichs2DConsP8")
FR_NS_RHS_PROV(2D, Cons, Cons, 9, 48, "NSFRLaxFriedrichs2DConsP9")
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
                                  const CFuint nbCells,
                                  const CFreal resFactor,
				  CFreal* states, 
                                  CFreal* gradients,
                                  CFreal* updateCoeff, 
				  CFreal* rhs,
                                  CFreal* solPntNormals,
                                  CFreal* flxPntNormals,
                                  CFint* faceDir,
                                  const CFuint nbSolPnts,
                                  const CFuint nbrFaces,
                                  const CFuint* faceFlxPntConn,
                                  const CFuint* stateIDs,
                                  const CFuint* neighbCellIDs,
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
                                  const CFreal* faceIntCoeff)
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
    
    // current cell data
    CellData cells(nbCells, stateIDs, neighbCellIDs, neighbFaceIDs, innerCellIsLeft, nbrFaces, nbSolPnts, ORDER);
    
    // get current cell
    CellData::Itr cell = cells.getItr(cellID);
          
    const CFuint nbFlxPntFlx = SCHEME::MODEL::NBEQS*(ORDER+1)*2*PHYS::DIM;//8;
    
    const CFuint nbFaceFlxPntFlx = SCHEME::MODEL::NBEQS*(ORDER+1);
   
    const CFuint nbrFaceFlxPnts = (ORDER+1);

    const CFuint totNbrFlxPnts = (ORDER+1)*2*PHYS::DIM;

    const CFuint nbNormals = PHYS::DIM*PHYS::DIM;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntFlx;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSol;
    
    flxPntFlx = 0.0;
    
    flxPntSol = 0.0;

    //CudaEnv::CFVecSlice<CFreal,nbrFaceFlxPnts> intCoeff(currFd.getFaceIntegrationCoef());

    //for (CFuint iFlxPnt = 0; iFlxPnt < nbrFaceFlxPnts; ++iFlxPnt)
    //{
      //intCoeff[iFlxPnt] = faceIntCoeff[iFlxPnt];
    //}
    //currFd.setFaceIntegrationCoef(iFlx,faceIntCoeff[iFlx]);
//if (cellID == 11) printf("iFlx , coeff\n");

 //     if (cellID == 1) printf("iFlx , coeff\n");
//if (cellID == 11) printf("hello %d\n", 0);
    // loop over sol pnts to compute flux
    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);
      //printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
    //if (cellID == 0) printf("GPUstate: %f %f %f %f\n", kd.states[0], kd.states[1], kd.states[2], kd.states[3]);
    
      setFluxData(stateID, cellID, &kd, &currFd, iSolPnt);

      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS*PHYS::DIM> grad(&gradients[stateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);

      CudaEnv::CFVecSlice<CFreal,nbNormals> n(&(kd.solPntNormals[stateID*nbNormals]));

      CudaEnv::CFVecSlice<CFreal,nbNormals> nFd(currFd.getScaledNormal(iSolPnt));
      
      CFuint k = 0;

      for (CFuint i = 0; i < nbNormals; ++i) 
      {
        nFd[i] = n[i];
      }

      // get the flux
      fluxScheme.prepareComputation(&currFd, &pmodel);

      fluxScheme(&currFd, &pmodel, iSolPnt);

/// add diff fluxes here to currFD!
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
            res[iEq] -= polyCoef*(currFd.getFlux(jSolIdx, iDir)[iEq])*resFactor;

//if (cellID == 11 && abs(polyCoef*(currFd.getFlux(jSolIdx, iDir)[iEq])) > 1e-8) printf("State: %d, jSol: %d, iDir: %d, var: %d, up: %f, poly: %f, flx: %f\n",iSolPnt,jSolIdx,iDir,iEq,polyCoef*(currFd.getFlux(iSolPnt, iDir)[iEq]),polyCoef,currFd.getFlux(jSolIdx, iDir)[iEq]);  
	  }
        }
      }
    }

    CudaEnv::CFVec<CFreal,nbFlxPntFlx*nbNormals> flxPntGrads;

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
          flxPntFlx[iFlxPnt*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*currFd.getFlux(solIdx, dim)[iEq];

          flxPntSol[iFlxPnt*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*states[cell.getStateID(solIdx)*SCHEME::MODEL::NBEQS+iEq];

          for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
          {

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
//if (cellID==11 && abs(flxPntFlx[flxIdx*SCHEME::MODEL::NBEQS+iVar] * divh) > 1e-8) printf("State: %d, flx: %d, var: %d, update: %f\n",iSolPnt,flxIdx,iVar,flxPntFlx[flxIdx*SCHEME::MODEL::NBEQS+iVar] * divh);  
        }
      }
    }

    // reset flx pnt fluxes  
    flxPntFlx = 0.0;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSolNeighb;
    
    for (CFuint iFlxPnt = 0; iFlxPnt < nbFlxPntFlx; ++iFlxPnt) {flxPntSolNeighb[iFlxPnt] = 0.0;}

    // current neighb cell data
    CellData cells2(nbCells, stateIDs, neighbCellIDs, neighbFaceIDs, innerCellIsLeft, nbrFaces, nbSolPnts, ORDER);

    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      const CFuint neighbCellID = cell.getNeighbCellID(iFace);  

      // get current cell
      CellData::Itr cell2 = cells2.getItr(neighbCellID);

      CFuint jFaceIdx = 100;

      for (CFuint jFace = 0; jFace < nbrFaces; ++jFace)
      {
        if (cell2.getNeighbCellID(jFace) == cellID)
        {
          jFaceIdx = jFace; 
          break;
        }
      }

      CFreal waveSpeedUpd = 0.0;

      if (jFaceIdx != 100)
      {
        const CFuint faceID = cell.getNeighbFaceID(iFace);

        const bool isLEFT = (bool) cell.getInnerCellIsLeft(iFace);

      // loop over face flx pnts
      for (CFuint iFlxPnt = 0; iFlxPnt < nbrFaceFlxPnts; ++iFlxPnt)
      { 
        // @TODO check if this also works for non QUADs
        const CFuint flxIdx = faceFlxPntConn[iFace*nbrFaceFlxPnts+iFlxPnt];
        const CFuint jFlxIdx = faceFlxPntConn[jFaceIdx*nbrFaceFlxPnts+nbrFaceFlxPnts-1-iFlxPnt];

        // loop over sol pnts to compute sol at flx pnt
        for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
        {
          const CFuint solIdx = flxSolDep[jFlxIdx*nbrFlxSolDep+iSolPnt]; 

          // Loop over conservative vars 
          for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
          {
            flxPntSolNeighb[flxIdx*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[jFlxIdx*nbSolPnts+solIdx]*states[cell2.getStateID(solIdx)*SCHEME::MODEL::NBEQS+iEq];
          }
        }

        for (CFuint iEq = 0; iEq < PHYS::NBEQS; ++iEq) 
        {
          currFd.getRstate(flxIdx)[iEq] = flxPntSolNeighb[flxIdx*PHYS::NBEQS+iEq];
        } 

        CudaEnv::CFVecSlice<CFreal,PHYS::DIM> n(&(kd.flxPntNormals[faceID*nbrFaceFlxPnts*PHYS::DIM+iFlxPnt*PHYS::DIM]));

        CudaEnv::CFVecSlice<CFreal,PHYS::DIM> nFd(currFd.getFlxScaledNormal(flxIdx));

        for (CFuint i = 0; i < PHYS::DIM; ++i) 
        {
          nFd[i] = n[i];
        }

        // get the flux
        fluxScheme.prepareComputation(&currFd, &pmodel);

        fluxScheme(&currFd, &pmodel, iFlxPnt, flxIdx, faceIntCoeff[iFlxPnt], isLEFT, waveSpeedUpd);

        // extrapolate the fluxes to the flux points
        for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
        {     
          const CFuint solIdx = flxSolDep[flxIdx*nbrFlxSolDep+iSolPnt];

          // get current state ID
          const CFuint stateID = cell.getStateID(solIdx);

          // get current vector slice out of rhs
          CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> res(&rhs[stateID*SCHEME::MODEL::NBEQS]);   
          
          // divergence of the correction function
          const CFreal divh = corrFctDiv[solIdx*nbrFlxPnts+flxIdx] * faceDir[cellID*totNbrFlxPnts+flxIdx];

          // Fill in the corrections
          for (CFuint iVar = 0; iVar < SCHEME::MODEL::NBEQS; ++iVar)
          {
            res[iVar] -= currFd.getInterfaceFlux(flxIdx)[iVar] * divh * resFactor;
//if (cellID==768) printf("resID: %d, State: %d, flx: %d, var: %d, updateFace: %f, flx %f, divh %f\n",stateID*SCHEME::MODEL::NBEQS+iVar,solIdx,flxIdx,iVar,-currFd.getInterfaceFlux(flxIdx)[iVar] * divh, currFd.getInterfaceFlux(flxIdx)[iVar],divh); 
          }
        }
      }

      //CFreal* waveSpeedUpd = currFd.getUpdateCoeff();

      for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
      {  
        // get current state ID
        const CFuint stateID = cell.getStateID(iSolPnt);

        updateCoeff[stateID] += waveSpeedUpd;
      }
 
      //currFd.resetUpdateCoeff();
      }
    }
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
                                       CFreal* solPntNormals,
                                       CFreal* flxPntNormals,
                                       CFint* faceDir,
                                       const CFuint nbSolPnts,
                                       const CFuint nbrFaces,
                                       const CFuint* faceFlxPntConn,
                                       const CFuint* stateIDs,
                                       const CFuint* neighbCellIDs,
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
                                       const CFreal* corrFctDiv)
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
    
    const CFuint nbFaceFlxPntFlx = SCHEME::MODEL::NBEQS*(ORDER+1);
   
    const CFuint nbrFaceFlxPnts = (ORDER+1);

    const CFuint totNbrFlxPnts = (ORDER+1)*2*PHYS::DIM;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntFlx;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSol;

    CudaEnv::CFVec<CFreal,SCHEME::MODEL::NBEQS> stateGradVars;
    //typename MathTypes<CFreal, GPU, SCHEME::MODEL::NBEQS>::VEC stateGradVars;
    
    flxPntFlx = 0.0;
    
    flxPntSol = 0.0;

    stateGradVars = 0.0;

    // loop over sol pnts to compute flux
    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);

      //typename MathTypes<CFreal, GPU, SCHEME::MODEL::NBEQS>::SLICEVEC currState(&states[stateID*SCHEME::MODEL::NBEQS]);

      CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[stateID*SCHEME::MODEL::NBEQS]);

      pmodelNS.getUpdateVS()->setGradientVars(&currState[0],&stateGradVars[0]);

      CudaEnv::CFVecSlice<CFreal,PHYS::DIM*PHYS::DIM> currNormals(&solPntNormals[stateID*PHYS::DIM*PHYS::DIM]);

      // get current vector slice out of gradients
      CudaEnv::CFVecSlice<CFreal,PHYS::DIM*SCHEME::MODEL::NBEQS> grad(&gradients[stateID*SCHEME::MODEL::NBEQS*PHYS::DIM]);

      // Loop over solution pnts to count the factor of all sol pnt polys
      for (CFuint jSolPnt = 0; jSolPnt < nbrSolSolDep; ++jSolPnt)
      { 
        const CFuint jSolIdx = solSolDep[iSolPnt*nbrSolSolDep+jSolPnt]; //(*m_solSolDep)[iSolPnt][jSolPnt];

        // Loop over deriv directions and sum them to compute divergence
        for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
        {
          const CFreal polyCoef = solPolyDerivAtSolPnts[iSolPnt*PHYS::DIM*nbSolPnts+iDir*nbSolPnts+jSolIdx];

          for (CFuint jDir = 0; jDir < PHYS::DIM; ++jDir)
          {
            // Loop over conservative fluxes 
            for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
            {
              //if (cellID == 11) printf("iSol: %d, iEq: %d, iDir: %d: %f\n", iSolPnt, iEq, jDir, polyCoef*currNormals[iDir*PHYS::DIM+jDir]*stateGradVars[iEq]); 

              // Store divFD in the vector that will be divFC
              grad[iEq*PHYS::DIM+jDir] += polyCoef*currNormals[iDir*PHYS::DIM+jDir]*stateGradVars[iEq];//*states[stateID*SCHEME::MODEL::NBEQS+iEq]; 

              //if (cellID == 11) printf("after  iSol: %d, iEq: %d, iDir: %d: %f\n", iSolPnt, iEq, jDir, grad[iEq*PHYS::DIM+jDir]); 
	    }
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

        CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> currState(&states[cell.getStateID(solIdx)*SCHEME::MODEL::NBEQS]);

        pmodelNS.getUpdateVS()->setGradientVars(&currState[0],&stateGradVars[0]);

        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
        {
          flxPntSol[iFlxPnt*SCHEME::MODEL::NBEQS+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*stateGradVars[iEq];          
        }
      }
    }

    // reset flx pnt fluxes  
    flxPntFlx = 0.0;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSolNeighb;
    
    for (CFuint iFlxPnt = 0; iFlxPnt < nbFlxPntFlx; ++iFlxPnt) {flxPntSolNeighb[iFlxPnt] = 0.0;}

    // current neighb cell data
    CellData cells2(nbCells, stateIDs, neighbCellIDs, neighbFaceIDs, innerCellIsLeft, nbrFaces, nbSolPnts, ORDER);

    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      const CFuint neighbCellID = cell.getNeighbCellID(iFace);  

      // get current cell
      CellData::Itr cell2 = cells2.getItr(neighbCellID);

      CFuint jFaceIdx = 100;

      for (CFuint jFace = 0; jFace < nbrFaces; ++jFace)
      {
        if (cell2.getNeighbCellID(jFace) == cellID)
        {
          jFaceIdx = jFace; 
          break;
        }
      }

      if (jFaceIdx != 100)
      {
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

            // divergence of the correction function
            const CFreal divh = corrFctDiv[solIdx*nbrFlxPnts+flxIdx];
           
            for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
            {
              const CFreal corrFactor = 0.5*(flxPntSolNeighb[flxIdx*SCHEME::MODEL::NBEQS+iEq]-flxPntSol[flxIdx*SCHEME::MODEL::NBEQS+iEq]);

              // Loop over deriv directions and sum them to compute divergence
              for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
              {
                for (CFuint jDir = 0; jDir < PHYS::DIM; ++jDir)
                {
                  // Store divFD in the vector that will be divFC
                  grad[iEq*PHYS::DIM+jDir] += divh*corrFactor*flxPntNormals[faceID*nbrFaceFlxPnts*PHYS::DIM+iFlxPnt*PHYS::DIM+jDir]*dirFactor; 

              //if (cellID == 11) printf("iSol: %d, iEq: %d, iFlx %d, iDir: %d: %e\n", iSolPnt, iEq, flxIdx, jDir, divh*corrFactor*flxPntNormals[faceID*nbrFaceFlxPnts*PHYS::DIM+iFlxPnt*PHYS::DIM+jDir]*dirFactor);  
	        }
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

      // Loop over deriv directions and sum them to compute divergence
      for (CFuint iDir = 0; iDir < PHYS::DIM; ++iDir)
      {
        CFreal nJacob2 = 0.0;

        for (CFuint jDir = 0; jDir < PHYS::DIM; ++jDir)
        {
          nJacob2 += solPntNormals[stateID*PHYS::DIM*PHYS::DIM+iDir*PHYS::DIM+jDir]*solPntNormals[stateID*PHYS::DIM*PHYS::DIM+iDir*PHYS::DIM+jDir];
        }

        const CFreal invJacob = 1/pow(nJacob2,0.5);

        for (CFuint jDir = 0; jDir < PHYS::DIM; ++jDir)
        {
          // Loop over conservative fluxes 
          for (CFuint iEq = 0; iEq < SCHEME::MODEL::NBEQS; ++iEq)
          {
            // Store divFD in the vector that will be divFC
            grad[iEq*PHYS::DIM+jDir] *= invJacob;  
	  }
        }
      }
        //if (cellID == 11) printf("iSol: %d, invJacob: %e\n", iSolPnt, temp);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

template <typename SCHEME, typename PHYSICS, typename PHYSICSNS, CFuint ORDER, CFuint NB_BLOCK_THREADS>
void ConvDiffRHSFluxReconstructionCUDA<SCHEME,PHYSICS,PHYSICSNS,ORDER,NB_BLOCK_THREADS>::execute()
{
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  
  CFTRACEBEGIN;
  
  CFLog(VERBOSE, "ConvDiffRHSFluxReconstructionCUDA::execute() START\n");
  
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
  
  if (m_onGPU) 
  {
#ifdef CF_HAVE_CUDA

    //CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
    //timer.start();
    
    // copy of data that change at every iteration
    socket_states.getDataHandle().getGlobalArray()->put(); 
    socket_gradientsCUDA.getDataHandle().getLocalArray()->put();
    socket_rhs.getDataHandle().getLocalArray()->put(); 
    socket_updateCoeff.getDataHandle().getLocalArray()->put();

    
    //CFLog(VERBOSE, "nb normals: " << socket_solPntNormals.getDataHandle().size() << ", n0: " << socket_solPntNormals.getDataHandle()[0] << "\n");

    socket_faceDir.getDataHandle().getLocalArray()->put();
    socket_solPntNormals.getDataHandle().getLocalArray()->put();
    socket_flxPntNormals.getDataHandle().getLocalArray()->put();


    DataHandle<Framework::State*, Framework::GLOBAL > statesI = socket_states.getDataHandle();
     
    //CFLog(VERBOSE, "ConvDiffRHSFluxReconstructionCUDA::execute() => CPU-->GPU data transfer took " << timer.elapsed() << " s\n");
    //timer.start();
    
    ConfigOptionPtr<SCHEME,  NOTYPE, GPU> dcof(lf);
    ConfigOptionPtr<typename PHYSICS::PTERM, NOTYPE, GPU> dcop(phys);
    ConfigOptionPtr<typename PHYSICSNS::DTERM, NOTYPE, GPU> dcopNS(physNS);
    ConfigOptionPtr<typename PHYSICSNS::PTERM, NOTYPE, GPU> dcopNSConv(physNSConv);

    const CFuint blocksPerGrid = CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
    const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
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
    
    // cudaFuncSetCacheConfig("computeFluxKernel", cudaFuncCachePreferL1);

    // if there is a diffusive term, compute the gradients
    if (hasDiffTerm)
    {
      computeGradientsKernel<FluxScheme,PHYSICS,PHYSICSNS,ORDER> <<<blocksPerGrid,nThreads>>>(
                                       dcop.getPtr(),
                                       dcopNS.getPtr(),
                                       dcopNSConv.getPtr(),
                                       nbCells,
				       socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
                                       gradients.getLocalArray()->ptrDev(), 
                                       solPntNormals.getLocalArray()->ptrDev(),
                                       flxPntNormals.getLocalArray()->ptrDev(),
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
                                       m_corrFctDiv2.ptrDev());
    }

    // compute the convective flux in each cell
    computeStateLocalRHSKernel<FluxScheme,PHYSICS,PHYSICSNS,ORDER> <<<blocksPerGrid,nThreads>>> 
      (dcof.getPtr(),
       dcop.getPtr(),
       dcopNS.getPtr(),
       nbCells,
       resFactor,
       socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
       gradients.getLocalArray()->ptrDev(), 
       updateCoeff.getLocalArray()->ptrDev(), 
       rhs.getLocalArray()->ptrDev(),
       solPntNormals.getLocalArray()->ptrDev(),
       flxPntNormals.getLocalArray()->ptrDev(),
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
       m_faceIntegrationCoefs2.ptrDev());
   
    cudaDeviceSynchronize();
    
    //for (CFuint i = 0; i < m_solPolyValsAtFlxPnts2.size(); ++i) {CFLog(INFO, "thing: " << m_solPolyValsAtFlxPnts2[i] << "\n");}
    
    //CFLog(INFO, "After Kernel, size: " << socket_states.getDataHandle().size() << "\n");
    
    //CFLog(VERBOSE, "ConvDiffRHSFluxReconstructionCUDA::execute() => computeFluxKernel took " << timer.elapsed() << " s\n");
    
    //for (CFuint i = 0; i < rhs.size(); ++i) {CFLog(INFO, "res before: " << rhs[i] << "\n");}
    
    //RealVector rhsB;
    //rhsB.resize(rhs.size());
    //for (CFuint i = 0; i < rhs.size(); ++i) {rhsB[i] = rhs[i];}
    
    //timer.start();
    rhs.getLocalArray()->get();
    updateCoeff.getLocalArray()->get();
    gradients.getLocalArray()->get();
    
    //for (CFuint i = 0; i < rhs.size(); ++i) {CFLog(INFO, "res after: " << rhs[i]-rhsB[i] << "\n");}
    //CFLog(VERBOSE, "ConvDiffRHSFluxReconstructionCUDA::execute() => GPU-->CPU data transfer took " << timer.elapsed() << " s\n");
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
  
  CFLog(VERBOSE, "ConvDiffRHSFluxReconstructionCUDA::execute() END\n");
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
