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
                                  const CFreal resFactor,
				  CFreal* states, 
                                  CFreal* updateCoeff, 
				  CFreal* rhs,
                                  CFreal* solPntNormals,
                                  CFreal* flxPntNormals,
                                  CFint* faceDir,
                                  const CFuint nbSolPnts,
                                  const CFuint nbrFaces,
                                  const CFuint* faceFlxPntConn,
				  const CFuint* cellInfo,
                                  const CFuint* stateIDs,
                                  const CFuint* neighbCellIDs,
                                  const CFuint* neighbFaceIDs,
                                  const CFuint dim,
                                  const CFuint nbrEqs,
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
    FluxData<typename SCHEME::MODEL> currFd; 

    // initialize flux data
    currFd.initialize();
    
    // physical model
    typename SCHEME::MODEL pmodel(dcop);
    SCHEME fluxScheme(dcof);
    
    // current cell data
    CellData cells(nbCells, cellInfo, stateIDs, neighbCellIDs, neighbFaceIDs, nbrFaces, nbSolPnts);
    
    // get current cell
    CellData::Itr cell = cells.getItr(cellID);
          
    const CFuint nbFlxPntFlx = SCHEME::MODEL::NBEQS*8;
    
    const CFuint nbFaceFlxPntFlx = SCHEME::MODEL::NBEQS*2;
   
    const CFuint nbrFaceFlxPnts = 2;

    const CFuint totNbrFlxPnts = 8;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntFlx;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSol;
    
    flxPntFlx = 0.0;
    
    flxPntSol = 0.0;

    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
//
      currFd.setFaceIntegrationCoef(iFlx,faceIntCoeff[iFlx]);
//if (cellID == 11) printf("iFlx , coeff\n");
    }
 //     if (cellID == 1) printf("iFlx , coeff\n");
//if (cellID == 11) printf("hello %d\n", 0);
    // loop over sol pnts to compute flux
    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
//ok
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);
      //printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
    //if (cellID == 0) printf("GPUstate: %f %f %f %f\n", kd.states[0], kd.states[1], kd.states[2], kd.states[3]);
    
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
//      if (cellID == 11) printf("before LF cell\n");
      fluxScheme(&currFd, &pmodel, false, 1, iSolPnt, cellID);
//      if (cellID == 11) printf("after LF cell\n");
//      // loop over sol pnts to compute flux
//      for (CFuint iDim = 0; iDim < dim; ++iDim)
//      {
//        if (cellID == 0) printf("HERE4 iSol: %d, iDim: %d, flux: %f %f %f %f \n", iSolPnt, iDim, currFd.getFlux(iSolPnt, iDim)[0], currFd.getFlux(iSolPnt, iDim)[1], currFd.getFlux(iSolPnt, iDim)[2], currFd.getFlux(iSolPnt, iDim)[3]);
//      }
    }

    for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
    {
      // get current state ID
      const CFuint stateID = cell.getStateID(iSolPnt);

      setFluxData(stateID, cellID, &kd, &currFd, iSolPnt);

      // Loop over solution pnts to count the factor of all sol pnt polys
      for (CFuint jSolPnt = 0; jSolPnt < nbrSolSolDep; ++jSolPnt)
      { 
        const CFuint jSolIdx = solSolDep[iSolPnt*nbrSolSolDep+jSolPnt]; //(*m_solSolDep)[iSolPnt][jSolPnt];

        // get current vector slice out of rhs
        CudaEnv::CFVecSlice<CFreal,SCHEME::MODEL::NBEQS> res(&rhs[stateID*SCHEME::MODEL::NBEQS]);
//if (cellID == 11) printf("resID: %d\n", stateID);
        // Loop over deriv directions and sum them to compute divergence
        for (CFuint iDir = 0; iDir < dim; ++iDir)
        {
          const CFreal polyCoef = solPolyDerivAtSolPnts[iSolPnt*dim*nbSolPnts+iDir*nbSolPnts+jSolIdx];//(*m_solPolyDerivAtSolPnts)[jSolPnt][iDir][iSolIdx]; 

          //if (cellID == 0) printf("polyCoef: %f\n", polyCoef);
          
          // Loop over conservative fluxes 
          for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
          {
//if (cellID == 11 && abs(polyCoef*(currFd.getFlux(jSolIdx, iDir)[iEq])) > 1e-8) printf("State: %d, jSol: %d, iDir: %d, var: %d, res before: %f\n",iSolPnt,jSolIdx,iDir,iEq,res[iEq]);

            // Store divFD in the vector that will be divFC
            res[iEq] -= polyCoef*(currFd.getFlux(jSolIdx, iDir)[iEq])*resFactor;

//if (cellID == 11 && abs(polyCoef*(currFd.getFlux(jSolIdx, iDir)[iEq])) > 1e-8) printf("State: %d, jSol: %d, iDir: %d, var: %d, up: %f, poly: %f, flx: %f\n",iSolPnt,jSolIdx,iDir,iEq,polyCoef*(currFd.getFlux(iSolPnt, iDir)[iEq]),polyCoef,currFd.getFlux(jSolIdx, iDir)[iEq]);  
            //if (cellID == 0) printf("res %f \n", res[iEq]);
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
        for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
        {
          flxPntFlx[iFlxPnt*nbrEqs+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*currFd.getFlux(solIdx, dim)[iEq];

          flxPntSol[iFlxPnt*nbrEqs+iEq] += solPolyValsAtFlxPnts[iFlxPnt*nbSolPnts+solIdx]*kd.states[iEq];          
        }
      }
    }

    // set extrapolated states
    for (CFuint iState = 0; iState < nbrFlxPnts; ++iState)
    {
      for (CFuint iEq = 0; iEq < PHYS::NBEQS; ++iEq) 
      {
        currFd.getLstate(iState)[iEq] = flxPntSol[iState*PHYS::NBEQS+iEq];
        //printf("stateL %d: %f\n", iEq, currFd.getLstate(iState)[iEq]);
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
        for (CFuint iVar = 0; iVar < nbrEqs; ++iVar)
        {
          res[iVar] += flxPntFlx[flxIdx*nbrEqs+iVar] * divh * resFactor;
//if (cellID==11 && abs(flxPntFlx[flxIdx*nbrEqs+iVar] * divh) > 1e-8) printf("State: %d, flx: %d, var: %d, update: %f\n",iSolPnt,flxIdx,iVar,flxPntFlx[flxIdx*nbrEqs+iVar] * divh);  
//if (cellID==11 && abs(flxPntFlx[flxIdx*nbrEqs+iVar] * divh) > 1e-8) printf("res after: %f\n",res[iVar]);  
        }
      }
    }

    // reset flx pnt fluxes  
    flxPntFlx = 0.0;
    
    //const CFuint nbFaceFlxPnts = nbFlxPntFlx/nbrFaces;
    
    CudaEnv::CFVec<CFreal,nbFlxPntFlx> flxPntSolNeighb;
    
    for (CFuint iFlxPnt = 0; iFlxPnt < nbFlxPntFlx; ++iFlxPnt) {flxPntSolNeighb[iFlxPnt] = 0.0;}

    // current neighb cell data
    CellData cells2(nbCells, cellInfo, stateIDs, neighbCellIDs, neighbFaceIDs, nbrFaces, nbSolPnts);

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
//if (cellID == 11) printf("bbbb\n");
      if (jFaceIdx != 100)
      {
        const CFuint faceID = cell.getNeighbFaceID(iFace);

      // loop over face flx pnts
      for (CFuint iFlxPnt = 0; iFlxPnt < nbrFaceFlxPnts; ++iFlxPnt)
      { 
        const CFuint flxIdx = faceFlxPntConn[iFace*nbrFaceFlxPnts+iFlxPnt];
        const CFuint jFlxIdx = faceFlxPntConn[jFaceIdx*nbrFaceFlxPnts+iFlxPnt];

        // loop over sol pnts to compute sol at flx pnt
        for (CFuint iSolPnt = 0; iSolPnt < nbrFlxSolDep; ++iSolPnt)
        {

          const CFuint solIdx = flxSolDep[jFlxIdx*nbrFlxSolDep+iSolPnt]; 
           //printf("flxIdx: %d\n", flxIdx);

          // Loop over conservative vars 
          for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
          {
            //printf("state %d before: %d, solID: %d, state: %f\n", iEq, cell2.getStateID(solIdx), flxIdx*nbrEqs+iEq, 0.0);
            flxPntSolNeighb[flxIdx*nbrEqs+iEq] += solPolyValsAtFlxPnts[jFlxIdx*nbSolPnts+solIdx]*states[cell2.getStateID(solIdx)*SCHEME::MODEL::NBEQS+iEq];
          }
        }

        // get current state ID
        //const CFuint neighbStateID = cell.getNeighbStateID(iFace,iSolPnt);
        
        //const CFuint neighbCellID = cell.getNeighbCellID(iFace);

        for (CFuint iEq = 0; iEq < PHYS::NBEQS; ++iEq) 
        {
          currFd.getRstate(flxIdx)[iEq] = flxPntSolNeighb[flxIdx*PHYS::NBEQS+iEq];
          //printf("stateR %d: %f\n", iEq, currFd.getRstate(flxIdx)[iEq]);
        } 

        CudaEnv::CFVecSlice<CFreal,PHYS::DIM> n(&(kd.flxPntNormals[faceID*nbrFaceFlxPnts*PHYS::DIM+iFlxPnt*PHYS::DIM]));

        CudaEnv::CFVecSlice<CFreal,PHYS::DIM> nFd(currFd.getFlxScaledNormal(flxIdx));
//if (cellID == 11) printf("faceID: %d, flx idx: %d, normal: %f, %f\n",faceID,flxIdx,n[0],n[1]);
        for (CFuint i = 0; i < PHYS::DIM; ++i) 
        {
          nFd[i] = n[i];
        }

//      if (cellID == 11) printf("aaa\n");
        // get the flux
        fluxScheme.prepareComputation(&currFd, &pmodel);
        //printf("flxIdx: %d\n", flxIdx);
//if (cellID == 11) printf("before LF face\n");

        fluxScheme(&currFd, &pmodel, true, iFlxPnt, flxIdx, cellID);
//if (cellID == 11) printf("after LF face\n");

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
//if (cellID == 11) printf("State %d, flx %d, divh %f\n",solIdx,flxIdx,divh);
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < nbrEqs; ++iVar)
          {
//if (cellID==1752) printf("resBefore: %f\n",res[iVar]);
            res[iVar] -= currFd.getInterfaceFlux(flxIdx)[iVar] * divh * resFactor;
//if (cellID == 11) printf("var %d, fI %f\n",iVar,currFd.getInterfaceFlux(flxIdx)[iVar]);
//if (cellID==1752) printf("resID: %d, State: %d, flx: %d, var: %d, updateFace: %f, res: %f\n",stateID*SCHEME::MODEL::NBEQS+iVar,solIdx,flxIdx,iVar,-currFd.getInterfaceFlux(flxIdx)[iVar] * divh, res[iVar]); 
          }
        }
      }

      CFreal* waveSpeedUpd = currFd.getUpdateCoeff();

      for (CFuint iSolPnt = 0; iSolPnt < nbSolPnts; ++iSolPnt)
      {  
        // get current state ID
        const CFuint stateID = cell.getStateID(iSolPnt);

        updateCoeff[stateID] += *waveSpeedUpd;

//if (cellID == 11) printf("iSol: %d, upd %f\n",iSolPnt,updateCoeff[stateID]);
      }
 
      currFd.resetUpdateCoeff();

//if (cellID == 11) printf("end face \n",cellID);
      }
    }
    //if (cellID == 0) printf("end sol pnt \n",cellID);
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
  
  //initializeComputationRHS();

  const CFuint nbStates = socket_states.getDataHandle().size();
  cf_assert(nbStates > 0);

  CFLog(INFO, "nbCells: " << nbCells << ", nbStates: " << nbStates << "\n");

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle(); 
  DataHandle<CFreal> solPntNormals = socket_solPntNormals.getDataHandle(); 
  DataHandle<CFreal> flxPntNormals = socket_flxPntNormals.getDataHandle(); 
  DataHandle<CFint> faceDir = socket_faceDir.getDataHandle();  

  SafePtr<SCHEME> lf  = getMethodData().getRiemannFlux().d_castTo<SCHEME>();
  SafePtr<typename PHYSICS::PTERM> phys = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<typename PHYSICS::PTERM>();
  
#ifdef CF_HAVE_CUDA
  typedef typename SCHEME::template DeviceFunc<GPU, PHYSICS> FluxScheme;  
#else
  typedef typename SCHEME::template DeviceFunc<CPU, PHYSICS> FluxScheme;
#endif 
  
  if (m_onGPU) 
  {
#ifdef CF_HAVE_CUDA

    CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
    timer.start();
    
    // copy of data that change at every iteration
    socket_states.getDataHandle().getGlobalArray()->put(); 
    socket_rhs.getDataHandle().getLocalArray()->put(); 
    socket_updateCoeff.getDataHandle().getLocalArray()->put();
    socket_faceDir.getDataHandle().getLocalArray()->put();
    
    CFLog(VERBOSE, "nb normals: " << socket_solPntNormals.getDataHandle().size() << ", n0: " << socket_solPntNormals.getDataHandle()[0] << "\n");
 
    socket_solPntNormals.getDataHandle().getLocalArray()->put();
    socket_flxPntNormals.getDataHandle().getLocalArray()->put();
    DataHandle<Framework::State*, Framework::GLOBAL > statesI = socket_states.getDataHandle();
     
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => CPU-->GPU data transfer took " << timer.elapsed() << " s\n");
    timer.start();
    
    ConfigOptionPtr<SCHEME,  NOTYPE, GPU> dcof(lf);
    ConfigOptionPtr<typename PHYSICS::PTERM, NOTYPE, GPU> dcop(phys);

    const CFuint blocksPerGrid = CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
    const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
    CFLog(VERBOSE, "blocksPerGrid: " << blocksPerGrid << ", threads: " << nThreads << "\n");

//CFuint megabytesToUse = 24;
//size_t newHeapSize = 1024 * 1000 * megabytesToUse;
//cudaDeviceSetLimit(cudaLimitMallocHeapSize, newHeapSize);
//printf("Adjusted heap size to be %d\n",(int) newHeapSize);

    //dim3 blocks(m_nbBlocksPerGridX, m_nbBlocksPerGridY);
    
    //cudaFuncSetCacheConfig("computeGradientsKernel", cudaFuncCachePreferL1);

    // get residual factor
    const CFreal resFactor = getMethodData().getResFactor();
    
    // cudaFuncSetCacheConfig("computeFluxKernel", cudaFuncCachePreferL1);
    CFLog(INFO, "Before Kernel: " << m_faceIntegrationCoefs2[0] << ", " << m_faceIntegrationCoefs2[1] << "\n");
    // compute the convective flux in each cell
    computeStateLocalRHSKernel<FluxScheme,PHYSICS> <<<blocksPerGrid,nThreads>>> 
      (dcof.getPtr(),
       dcop.getPtr(),
       nbCells,
       resFactor,
       socket_states.getDataHandle().getGlobalArray()->ptrDev(), 
       updateCoeff.getLocalArray()->ptrDev(), 
       rhs.getLocalArray()->ptrDev(),
       solPntNormals.getLocalArray()->ptrDev(),
       flxPntNormals.getLocalArray()->ptrDev(),
       faceDir.getLocalArray()->ptrDev(),
       m_nbrSolPnts,
       4,
       m_faceFlxPntConn2.ptrDev(),
       m_cellInfo.ptrDev(),
       m_stateIDs.ptrDev(),
       m_neighbCellIDs.ptrDev(),
       m_neighbFaceIDs.ptrDev(),
       m_dim,
       m_nbrEqs,
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
    
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => computeFluxKernel took " << timer.elapsed() << " s\n");
    
    //for (CFuint i = 0; i < rhs.size(); ++i) {CFLog(INFO, "res before: " << rhs[i] << "\n");}
    
    //RealVector rhsB;
    //rhsB.resize(rhs.size());
    //for (CFuint i = 0; i < rhs.size(); ++i) {rhsB[i] = rhs[i];}
    
    timer.start();
    rhs.getLocalArray()->get();
    updateCoeff.getLocalArray()->get();
    
    //for (CFuint i = 0; i < rhs.size(); ++i) {CFLog(INFO, "res after: " << rhs[i]-rhsB[i] << "\n");}
    CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() => GPU-->CPU data transfer took " << timer.elapsed() << " s\n");
    //CFLog(INFO, "resSize: " << rhs.size() << "\n");
    for (CFuint i = 0; i < rhs.size(); ++i)
    {
      //if (abs(rhs[i]) > 1.0e-10) CFLog(INFO, "res " << i << ": " << rhs[i] << "\n");
    }

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
if (m_cells[LEFT ]->getID()==11 || m_cells[RIGHT ]->getID()==11) CFLog(INFO, "faceID: " << m_face->getID() << "\n");
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

DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
//CFLog(INFO, "resSize: " << rhs.size() << "\n");
for (CFuint i = 0; i < rhs.size(); ++i)
 {
    //if (abs(rhs[i]) > 1.0e-10) CFLog(INFO, "res " << i << ": " << rhs[i] << "\n");
  }
  }
  
  //finalizeComputationRHS();
  
  CFLog(VERBOSE, "ConvRHSFluxReconstructionCUDA::execute() END\n");
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
