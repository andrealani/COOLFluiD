#include <fstream>
#include <iostream>

#include "Common/PE.hh"
#include "Common/BadValueException.hh"
#include "Common/CFPrintContainer.hh"
#include "Framework/CudaDeviceManager.hh"
#include "Common/CUDA/CFVec.hh"
#include "Framework/CudaTimer.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/SocketBundleSetter.hh"

#include "FiniteVolume/CellCenterFVM.hh"

#include "RadiativeTransfer/RadiativeTransfer.hh"
#include "RadiativeTransfer/Solvers/FiniteVolumeDOM/RadiativeTransferFVDOMCUDA.hh"
#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"

// the following depends on the GPU model, should be user-defined
#define SHARED_ARRAY_SIZE 1536 

/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RadiativeTransferFVDOMCUDA, 
		      DataProcessingData, 
		      RadiativeTransferModule>
radiativeTransferFVDOMCUDAProvider("RadiativeTransferFVDOMCUDA");

//////////////////////////////////////////////////////////////////////////////

//this useful code sample is taken from 
 //http://stackoverflow.com/questions/12626096/why-has-atomicadd-not-been-implemented-for-doubles
inline __device__ double atomicAdd(double* address, double val)
{
  unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
		    __double_as_longlong(val +
					 __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

//////////////////////////////////////////////////////////////////////////////

__global__ void getFieldOpacitiesKernel(const bool useExponentialMethod,
					const CFuint TID, 
					const CFuint PID,
					const CFuint nbTemp,
					const CFuint nbPress,
					const CFuint nbBins,
					const CFuint ib,
					const CFuint nbEqs,
					const CFuint nbCells,
					const CFreal* Ttable,
					const CFreal* Ptable,
					const CFreal* states,
					const CFreal* volumes,
					const CFreal* opacities,
					const CFreal* radSource,
					CFreal* fieldSource,
					CFreal* fieldAbsor,
					CFreal* fieldAbSrcV,
					CFreal* fieldAbV)
{    
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (cellID < nbCells) {
    fieldSource[cellID] = 0.;
    if(useExponentialMethod) {
      fieldAbsor[cellID]  = 0.;
    }
    else {
      fieldAbSrcV[cellID] = 0.;
      fieldAbV[cellID]    = 0.;
    }
    
    //Get the field pressure and T commented because now we impose a temperature profile
    const CFuint sIdx = cellID*nbEqs; 
    const CFreal p = states[sIdx + PID];
    const CFreal T = states[sIdx + TID];
    const CFreal patm = p/101325.; //converting from Pa to atm
    
    CFreal val1 = 0;
    CFreal val2 = 0;
    
    RadiativeTransferFVDOM::DeviceFunc interp;
    interp.tableInterpolate(nbBins, nbTemp, nbPress, Ttable, Ptable,
			    opacities, radSource, T, patm, ib, val1, val2); 
    
    if(useExponentialMethod){
      if (val1 <= 1e-30 || val2 <= 1e-30 ){
	fieldSource[cellID] = 1e-30;
	fieldAbsor[cellID]  = 1e-30;
      }
      else {
	fieldSource[cellID] = val2/val1;
	fieldAbsor[cellID]  = val1;
      }
    }
    else{
      if (val1 <= 1e-30 || val2 <= 1e-30 ){
	fieldSource[cellID] = 1e-30;
	fieldAbV[cellID]    = 1e-30*volumes[cellID]; // Volume converted from m^3 into cm^3
      }
      else {
	fieldSource[cellID] = val2/val1;
	fieldAbV[cellID]    = val1*volumes[cellID];
      }      
      fieldAbSrcV[cellID]   = fieldSource[cellID]*fieldAbV[cellID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

__global__ void computeQKernelExponentialDirBigMem(const CFuint dStart,
						   const CFuint dEnd,
						   const CFuint nbCells,
						   const CFuint* cellFaces,
						   const CFint* faceCell,
						   const CFuint* nbFacesInCell,
						   const CFint* isOutward,
						   const CFint* advanceOrder,
						   const CFreal* weightIn,
						   const CFreal* volumes,
						   const CFreal* fieldSource,
						   const CFreal* fieldAbsor,
						   const CFreal* normals,
						   const CFreal* mdirs,
						   CFreal* In,
						   CFreal* divq,
						   CFreal* qx, CFreal* qy, CFreal* qz)
{ 
  // mdirs is used often, by all threads and should be shared
  //  __shared__ CFreal mdirs[256]; // overallocated memory
  const CFuint nbDirs = dEnd - dStart; 
  //CFint tID = threadIdx.x; 
  //while (tID < nbDirs) {
  //  mdirs[dStart+tID] = mdirsIn[dStart+tID];
  //  tID += blockDim.x;
  // }
  // __syncthreads();
  
  // each thread takes care of computing the gradient for one single cell
  const CFuint dirID = threadIdx.x + blockIdx.x*blockDim.x;
  if (dirID < nbDirs) {
    RadiativeTransferFVDOM::DeviceFunc fun;
    const CFuint d = dStart + dirID;
    const CFreal weight = weightIn[d];
    const CFuint d3 = d*3;
    const CFreal mdirs30 = mdirs[d3];
    const CFreal mdirs31 = mdirs[d3+1];
    const CFreal mdirs32 = mdirs[d3+2]; 
    const CFuint startCell = dirID*nbCells;
    for (CFuint m = 0; m < nbCells; ++m) {
      CFreal inDirDotnANeg = 0.;
      CFreal Ic            = 0.;
      CFreal dirDotnANeg   = 0.;
      CFreal Lc            = 0.;
      CFreal halfExp       = 0.;
      
      // allocate the cell entity
      const CFuint iCell   = abs(advanceOrder[startCell+m]);
      const CFuint nbFaces = nbFacesInCell[iCell];
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) { 
	const CFuint faceID = cellFaces[iFace*nbCells + iCell];
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFuint startID = faceID*3;
	const CFreal dotProdInFace = fun.getDirDotNA(d,&mdirs[0],&normals[startID]);
	const CFreal dirDotNA = dotProdInFace*factor;
	
	if(dirDotNA < 0.) {
	  dirDotnANeg += dirDotNA;
	  const CFint fcellID = faceCell[faceID*2]; 
	  const CFint neighborCellID = (fcellID == iCell) ? faceCell[faceID*2+1] : fcellID;
	  const CFreal source = (neighborCellID >=0) ? In[neighborCellID*nbDirs+dirID] : fieldSource[iCell];
	  inDirDotnANeg += source*dirDotNA;
	   
	  /*if (iCell==100 && dirID == 0) {
	    printf ("source   : %6.6f \n", source);
	    printf ("dirDotNA : %6.6f  \n", dirDotNA);
	    printf ("inDirDotnANeg : %6.6f \n",inDirDotnANeg);
	    printf ("factor   : %6.6f  \n", factor);
	    }*/
	}
      } 
      
      Lc        = volumes[iCell]/(- dirDotnANeg); 
      halfExp   = exp(-0.5*Lc*fieldAbsor[iCell]);
      const CFreal InCell = (inDirDotnANeg/dirDotnANeg)*halfExp*halfExp + 
	(1. - halfExp*halfExp)*fieldSource[iCell];
      Ic = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*fieldSource[iCell];

      CFreal inDirDotnA = inDirDotnANeg;
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = cellFaces[iFace*nbCells + iCell];
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFuint startID = faceID*3;
	const CFreal dotProdInFace = fun.getDirDotNA(d,&mdirs[0],&normals[startID]);
	const CFreal dirDotNA = dotProdInFace*factor;
	if (dirDotNA > 0.) {
	  inDirDotnA += InCell*dirDotNA;
	}
      }
      
      const CFuint cellIDin = iCell*nbDirs + dirID; // dirID*nbCells+iCell;
      In[cellIDin] = InCell; 
      const CFreal IcWeight = Ic*weight;
      
      // here no atomics are needed since qi/divq is a (logically) 2D array, with one qi/diq array for each threads
      qx[cellIDin]    += mdirs30*IcWeight;
      qy[cellIDin]    += mdirs31*IcWeight;
      qz[cellIDin]    += mdirs32*IcWeight;
      divq[cellIDin]  += inDirDotnA*weight;
     
      /*if (iCell==100 && dirID == 0) {
        printf ("IcWeight    : %6.6f \n", IcWeight);
        printf ("inDirDotnA  : %6.6f \n",inDirDotnA);
        printf ("InCell      : %6.6f \n", InCell);
        printf ("cellIDin    : %d  \n", cellIDin);
        const CFreal qxIcell = qx[cellIDin];
        printf ("qx[iCell]   : %6.6f  \n", qxIcell);
        const CFreal divqIcell = divq[cellIDin];
        printf ("divq[iCell] : %6.6f  \n", divqIcell);
        const CFreal In0 = In[cellIDin];
        printf ("In[iCell]   : %6.6f  \n", In0);
        printf ("d3          : %d  \n", d3);
        printf ("mdirs[d3]   : %6.6f  \n", mdirs30);
        printf ("mdirs[d3+1] : %6.6f  \n", mdirs31);
        printf ("mdirs[d3+2] : %6.6f  \n", mdirs32);
      }*/
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

__global__ void computeQKernelExponentialDirAtomic(const CFuint dStart,
						   const CFuint dEnd,
						   const CFuint nbCells,
						   const CFuint* cellFaces,
						   const CFint* faceCell,
						   const CFuint* nbFacesInCell,
						   const CFint* isOutward,
						   const CFint* advanceOrder,
						   const CFreal* weightIn,
						   const CFreal* volumes,
						   const CFreal* fieldSource,
						   const CFreal* fieldAbsor,
						   const CFreal* normals,
						   const CFreal* mdirs,
						   CFreal* In,
						   CFreal* divq,
						   CFreal* qx, CFreal* qy, CFreal* qz)
{ 
  // mdirs is used often, by all threads and should be shared 
  //__shared__ CFreal mdirs[256]; // overallocated memory
  const CFuint nbDirs = dEnd - dStart; 
 // CFint tID = threadIdx.x; 
 // while (tID < nbDirs) {
 //   mdirs[dStart+tID] = mdirsIn[dStart+tID];
 //   tID += blockDim.x;
 // }
 // __syncthreads();
  
  // each thread takes care of computing the gradient for one single cell
  const CFuint dirID = threadIdx.x + blockIdx.x*blockDim.x;
  if (dirID < nbDirs) {
    RadiativeTransferFVDOM::DeviceFunc fun;
    const CFuint d = dStart + dirID;
    const CFuint d3 = d*3;
    const CFreal mdirs30 = mdirs[d3];
    const CFreal mdirs31 = mdirs[d3+1];
    const CFreal mdirs32 = mdirs[d3+2];
    const CFreal weight = weightIn[d];
    const CFuint startCell = dirID*nbCells;
     
    for (CFuint m = 0; m < nbCells; ++m) {
      CFreal inDirDotnANeg = 0.;
      CFreal Ic            = 0.;
      CFreal dirDotnANeg   = 0.;
      CFreal Lc            = 0.;
      CFreal halfExp       = 0.;
     
      // allocate the cell entity
      const CFuint iCell   = abs(advanceOrder[startCell+m]);
      const CFuint nbFaces = nbFacesInCell[iCell];
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) { 
	const CFuint faceID = cellFaces[iFace*nbCells + iCell];
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFuint startID = faceID*3;
	const CFreal dotProdInFace = fun.getDirDotNA(d,&mdirs[0],&normals[startID]);
	const CFreal dirDotNA = dotProdInFace*factor;
	
	if(dirDotNA < 0.) {
	  dirDotnANeg += dirDotNA;
	  const CFint fcellID = faceCell[faceID*2]; 
	  const CFint neighborCellID = (fcellID == iCell) ? faceCell[faceID*2+1] : fcellID;
	  const CFreal source = (neighborCellID >=0) ? In[neighborCellID*nbDirs+dirID] : fieldSource[iCell];
	  inDirDotnANeg += source*dirDotNA;
	  
	  /*if (iCell==100 && dirID == 0) {
	    printf ("source   : %6.6f \n", source);
	    printf ("dirDotNA : %6.6f  \n", dirDotNA);
	    printf ("inDirDotnANeg : %6.6f \n",inDirDotnANeg);
	    printf ("factor   : %6.6f  \n", factor);
	    }*/
	  
	}
      } 
      
      Lc        = volumes[iCell]/(- dirDotnANeg); 
      halfExp   = exp(-0.5*Lc*fieldAbsor[iCell]);
      const CFreal InCell = (inDirDotnANeg/dirDotnANeg)*halfExp*halfExp + 
	(1. - halfExp*halfExp)*fieldSource[iCell];
      Ic = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*fieldSource[iCell];
      
      CFreal inDirDotnA = inDirDotnANeg;
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = cellFaces[iFace*nbCells + iCell];
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFuint startID = faceID*3;
	const CFreal dotProdInFace = fun.getDirDotNA(d,&mdirs[0],&normals[startID]);
	const CFreal dirDotNA = dotProdInFace*factor;
	if (dirDotNA > 0.) {
	  inDirDotnA += InCell*dirDotNA;
	}
      }
      
      const CFuint cellIDin = iCell*nbDirs + dirID; // dirID*nbCells+iCell;
      In[cellIDin] = InCell;
      const CFreal IcWeight = Ic*weight;
      const CFreal m0w = mdirs30*IcWeight;
      const CFreal m1w = mdirs31*IcWeight;
      const CFreal m2w = mdirs32*IcWeight;
      const CFreal inw = inDirDotnA*weight;    
  
      // here atomics are needed since qi/divq are single arrays shared by all threads
      atomicAdd(&qx[iCell], m0w);
      atomicAdd(&qy[iCell], m1w);          
      atomicAdd(&qz[iCell], m2w);
      atomicAdd(&divq[iCell], inw);
      
      /*if (iCell==100 && dirID == 0) {
	printf ("IcWeight    : %6.6f \n", IcWeight);
	printf ("inDirDotnA  : %6.6f \n",inDirDotnA);
	printf ("InCell      : %6.6f \n", InCell);
      	printf ("cellIDin    : %d  \n", cellIDin);
	const CFreal qxIcell = qx[iCell];
	printf ("qx[iCell]   : %6.6f  \n", qxIcell);
	const CFreal divqIcell = divq[iCell];
	printf ("divq[iCell] : %6.6f  \n", divqIcell);
	const CFreal In0 = In[cellIDin];
	printf ("In[iCell]   : %6.6f  \n", In0);
	printf ("d3          : %d  \n", d3);
	printf ("mdirs[d3]   : %6.6f  \n", mdirs30);
	printf ("mdirs[d3+1] : %6.6f  \n", mdirs31);
	printf ("mdirs[d3+2] : %6.6f  \n", mdirs32);
      }*/
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

__global__ void computeQKernelExponentialDir(const CFuint dStart,
					     const CFuint dEnd,
					     const CFuint startCellID,
					     const CFuint endCellID,
					     const CFuint nbCells,
					     const CFuint* cellFaces,
					     const CFint* faceCell,
					     const CFuint* nbFacesInCell,
					     const CFint* isOutward,
					     const CFint* advanceOrder,
					     const CFreal* weightIn,
					     const CFreal* volumes,
					     const CFreal* fieldSource,
					     const CFreal* fieldAbsor,
					     const CFreal* normals,
					     const CFreal* mdirs,
					     CFreal* In,
					     CFreal* II,
					     CFreal* divq,
					     CFreal* qx, CFreal* qy, CFreal* qz)
{ 
  // mdirs is used often and should be shared 
  // __shared__ CFreal mdirs[2048]; // overallocated memory
  
  const CFuint nbDirs = dEnd - dStart; 
  // while (tID < nbDirs) {
  //   mdirs[dStart+tID] = mdirsIn[dStart+tID];
  //   tID += blockIdx.x*blockDim.x; // this could be buggy
  // }
  // __syncthreads();
  
  // this assumes 49,152 bytes of shared memory
  __shared__ CFreal qxSh[SHARED_ARRAY_SIZE];
  __shared__ CFreal qySh[SHARED_ARRAY_SIZE]; 
  __shared__ CFreal qzSh[SHARED_ARRAY_SIZE];
  __shared__ CFreal divqSh[SHARED_ARRAY_SIZE];
  
  CFint tID = threadIdx.x; 
  while (tID < SHARED_ARRAY_SIZE) {
    qxSh[tID] = qySh[tID] = qzSh[tID] = divqSh[tID] = 0.;
    tID += blockDim.x;
  }
  __syncthreads();
  
 // printf ("BBBBBBBB   : %d %d \n", blockDim.x, blockIdx.x);
  
  CFuint iCellList[SHARED_ARRAY_SIZE];
  
  // each thread takes care of computing the gradient for one single cell
  const CFuint dirID = threadIdx.x + blockIdx.x*blockDim.x;
  if (dirID < nbDirs) {
    RadiativeTransferFVDOM::DeviceFunc fun;
    const CFuint d = dStart + dirID;
    const CFreal weight = weightIn[d];
    const CFuint startCell = dirID*nbCells;
    CFuint iCellSh = 0;
    for (CFuint m = startCellID; m < endCellID; ++m, ++iCellSh) {
      CFreal inDirDotnANeg = 0.;
      CFreal Ic            = 0.;
      CFreal dirDotnANeg   = 0.;
      CFreal Lc            = 0.;
      CFreal halfExp       = 0.;
      
      // allocate the cell entity
      const CFuint iCell   = abs(advanceOrder[startCell+m]);
      iCellList[iCellSh] = iCell;
      const CFuint nbFaces = nbFacesInCell[iCell];
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) { 
	const CFuint faceID = cellFaces[iFace*nbCells + iCell];
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFuint startID = faceID*3;
	const CFreal dotProdInFace = fun.getDirDotNA(d,&mdirs[0],&normals[startID]);
	const CFreal dirDotNA = dotProdInFace*factor;
	
	if(dirDotNA < 0.) {
	  dirDotnANeg += dirDotNA;
	  const CFint fcellID = faceCell[faceID*2]; 
	  const CFint neighborCellID = (fcellID == iCell) ? faceCell[faceID*2+1] : fcellID;
	  const CFreal source = (neighborCellID >=0) ? In[nbCells*dirID+neighborCellID] : fieldSource[iCell];
	  inDirDotnANeg += source*dirDotNA;
	  
	  /*if (iCell==100) {
	    printf ("source   : %6.6f \n", source);
	    printf ("dirDotNA : %6.6f  \n", dirDotNA);
	    printf ("inDirDotnANeg : %6.6f \n",inDirDotnANeg);
	    printf ("factor   : %6.6f  \n", factor);
	    }*/
	}
      } 
      const CFuint cellIDin = nbCells*dirID+iCell;
      Lc        = volumes[iCell]/(- dirDotnANeg); 
      halfExp   = exp(-0.5*Lc*fieldAbsor[iCell]);
      In[cellIDin] = (inDirDotnANeg/dirDotnANeg)*halfExp*halfExp + (1. - halfExp*halfExp)*fieldSource[iCell];
      Ic        = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*fieldSource[iCell];
      
      /*if (iCell==100) {
	printf ("Lc   : %6.6f  \n", Lc);
	printf ("halfExp : %6.6f  \n", halfExp);
	printf ("In : %6.6f \n",In[cellIDin]);
	printf ("Ic   : %6.6f  \n", Ic);
	printf ("weight   : %6.6f \n", weight);
	}*/

      CFreal inDirDotnA = inDirDotnANeg;
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = cellFaces[iFace*nbCells + iCell];
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFuint startID = faceID*3;
	const CFreal dotProdInFace = fun.getDirDotNA(d,&mdirs[0],&normals[startID]);
	const CFreal dirDotNA = dotProdInFace*factor;
	if (dirDotNA > 0.) {
	  inDirDotnA += In[cellIDin]*dirDotNA;
	}
      }
      
      const CFreal IcWeight = Ic*weight;
      II[cellIDin] += IcWeight;
      const CFuint d3 = d*3;
      atomicAdd(&qxSh[iCellSh], mdirs[d3]*IcWeight);
      atomicAdd(&qySh[iCellSh], mdirs[d3+1]*IcWeight);          
      atomicAdd(&qzSh[iCellSh], mdirs[d3+2]*IcWeight);
      atomicAdd(&divqSh[iCellSh], inDirDotnA*weight);
    }
  }
  
  __syncthreads();
  
  tID = threadIdx.x; 
  const CFuint nbCellsSh = endCellID - startCellID;
  // printf ("start, end, nbCellsSh   : %d %d %d\n", startCellID, endCellID, nbCellsSh);
 
  // at this point should be iCellSh == nbCellsSh
  while (tID < nbCellsSh) {
    // the following could be done outside the loop
    const CFuint iCell = iCellList[tID];
    atomicAdd(&qx[iCell], qxSh[tID]);
    atomicAdd(&qy[iCell], qySh[tID]);
    atomicAdd(&qz[iCell], qzSh[tID]);
    atomicAdd(&divq[iCell], divqSh[tID]);
    tID += blockDim.x; 
  }
}
  
//////////////////////////////////////////////////////////////////////////////
 
__global__ void computeQKernelNoExponentialDir(const CFuint dStart,
					       const CFuint dEnd,
					       const CFuint nbCells,
					       const CFuint* cellFaces,
					       const CFint* faceCell,
					       const CFuint* nbFacesInCell,
					       const CFint* isOutward,
					       const CFint* advanceOrder,
					       const CFreal* weightIn,
					       const CFreal* volumes,
                                               const CFreal* fieldSource,
                                               const CFreal* fieldAbSrcV,
                                               const CFreal* fieldAbV,
					       const CFreal* normals,
					       const CFreal* mdirsIn,
					       CFreal* In,
					       //					       CFreal* II,
					       CFreal* divQ,
					       CFreal* qx, CFreal* qy, CFreal* qz)
{ 
  // mdirs is used often and should be shared 
  __shared__ CFreal mdirs[2048]; // overallocated memory
  CFint tID = threadIdx.x; 
  const CFuint nbDirs = dEnd - dStart; 
  while (tID < nbDirs) {
    mdirs[dStart+tID] = mdirsIn[dStart+tID];
    tID += blockDim.x; // this could be buggy
  }
  __syncthreads();
  
  // each thread takes care of computing the gradient for one single cell
  const CFuint dirID = threadIdx.x + blockIdx.x*blockDim.x;
  if (dirID < nbDirs) {
    RadiativeTransferFVDOM::DeviceFunc fun;
    const CFuint d = dStart + dirID;
    const CFreal weight = weightIn[d];
    const CFuint startCell = dirID*nbCells;
    for (CFuint m = 0; m < nbCells; ++m) {    
      // allocate the cell entity
      const CFuint iCell   = abs(advanceOrder[startCell+m]);
      CFreal inDirDotnANeg = 0.;
      CFreal Ic            = 0.;
      CFreal dirDotnAPos   = 0.;
      
      const CFuint nbFaces = nbFacesInCell[iCell];
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = cellFaces[iFace*nbCells + iCell];
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFuint startID = faceID*3;
	const CFreal dotProdInFace = fun.getDirDotNA(d,&mdirs[0],&normals[startID]);
	const CFreal dirDotNA = dotProdInFace*factor;
	
	if (dirDotNA >= 0.){
	  dirDotnAPos += dirDotNA;
	}
	else {
	  const CFint fcellID = faceCell[faceID*2]; 
	  const CFint neighborCellID = (fcellID == iCell) ? faceCell[faceID*2+1] : fcellID;
	  const CFreal source = (neighborCellID >=0) ? In[nbCells*dirID+neighborCellID] : fieldSource[iCell];
	  inDirDotnANeg += source*dirDotNA;
	}
      } 
      const CFuint cellIDin = nbCells*dirID+iCell;
      In[cellIDin] = (fieldAbSrcV[iCell] - inDirDotnANeg)/(fieldAbV[iCell] + dirDotnAPos);
      Ic = In[cellIDin];
      
      qx[iCell] += Ic*mdirs[d*3]*weight;
      qy[iCell] += Ic*mdirs[d*3+1]*weight;
      qz[iCell] += Ic*mdirs[d*3+2]*weight;
      
      CFreal inDirDotnA = inDirDotnANeg;
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = cellFaces[iFace*nbCells + iCell];
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFuint startID = faceID*3;
	const CFreal dotProdInFace = fun.getDirDotNA(d,&mdirs[0],&normals[startID]);
	const CFreal dirDotNA = dotProdInFace*factor;
	if (dirDotNA > 0.) {
	  inDirDotnA += In[cellIDin]*dirDotNA;
	}
      }
      
      divQ[iCell] += inDirDotnA*weight;
      //  II[nbCells*dirID+iCell] += Ic*weight;
    }  
  }
}
      
//////////////////////////////////////////////////////////////////////////////
 
__global__ void computeQKernelExponential(const CFuint d,
					  const CFuint nbCells,
					  const CFreal weightIn,
					  const CFuint* cellFaces,
					  const CFint* faceCell,
					  const CFuint* nbFacesInCell,
					  const CFint* isOutward,
					  const CFint* advanceOrder,
					  const CFreal* volumes,
					  const CFreal* fieldSource,
					  const CFreal* fieldAbsor,
					  const CFreal* dotProdInFace,
					  const CFreal* mdirs,
					  CFreal* In,
					  CFreal* II,
					  CFreal* divQ,
					  CFreal* qx, CFreal* qy, CFreal* qz)
{    
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  __shared__ CFreal weight;
  weight = weightIn;
  __syncthreads();
  
  if (cellID < nbCells) {
    CFreal inDirDotnANeg = 0.;
    CFreal Ic            = 0.;
    CFreal dirDotnANeg   = 0.;
    CFreal Lc            = 0.;
    CFreal halfExp       = 0.;
        
    // allocate the cell entity
    const CFuint iCell   = abs(advanceOrder[d*nbCells+cellID]);
    
    const CFuint nbFaces = nbFacesInCell[iCell];
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) { 
      const CFuint faceID = cellFaces[iFace*nbCells + iCell];
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = dotProdInFace[faceID]*factor;
      
      if(dirDotNA < 0.) {
	dirDotnANeg += dirDotNA;
	
        const CFint fcellID = faceCell[faceID*2]; 
        const CFint neighborCellID = (fcellID == iCell) ? faceCell[faceID*2+1] : fcellID;
	const CFreal source = (neighborCellID >=0) ? In[neighborCellID] : fieldSource[iCell];
        inDirDotnANeg += source*dirDotNA;
      

	if (iCell==100) {
	  printf ("source   : %6.6f \n", source);
	  printf ("dirDotNA : %6.6f  \n", dirDotNA);
	  printf ("inDirDotnANeg : %6.6f \n",inDirDotnANeg);
	  printf ("factor   : %6.6f  \n", factor);
	}
	
	
      }
    } 
    Lc        = volumes[iCell]/(- dirDotnANeg); 
    halfExp   = exp(-0.5*Lc*fieldAbsor[iCell]);
    In[iCell] = (inDirDotnANeg/dirDotnANeg)*halfExp*halfExp + (1. - halfExp*halfExp)*fieldSource[iCell];
    Ic        = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*fieldSource[iCell];
    
    if (iCell==100) {
      printf ("Lc   : %6.6f  \n", Lc);
      printf ("halfExp : %6.6f  \n", halfExp);
      printf ("In : %6.6f \n",In[iCell]);
      printf ("Ic   : %6.6f  \n", Ic);
      printf ("weight   : %6.6f \n", weight);
    }
    
    qx[iCell] += Ic*mdirs[d*3]*weight;
    qy[iCell] += Ic*mdirs[d*3+1]*weight;
    qz[iCell] += Ic*mdirs[d*3+2]*weight;
    
    CFreal inDirDotnA = inDirDotnANeg;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID = cellFaces[iFace*nbCells + iCell];
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = dotProdInFace[faceID]*factor;
      if (dirDotNA > 0.) {
	inDirDotnA += In[iCell]*dirDotNA;
      }
    }
    
    divQ[iCell] += inDirDotnA*weight;
    II[iCell]   += Ic*weight;
  }  
}
      
//////////////////////////////////////////////////////////////////////////////

__global__ void computeQKernelNoExponential(const CFuint d, 
					    const CFuint nbCells,
					    const CFreal weightIn,
					    const CFuint* cellFaces,
					    const CFint* faceCell,
					    const CFuint* nbFacesInCell,
					    const CFint* isOutward,
					    const CFint* advanceOrder,
					    const CFreal* volumes,
					    const CFreal* fieldSource,
					    const CFreal* fieldAbSrcV,
					    const CFreal* fieldAbV,
					    const CFreal* dotProdInFace,
					    const CFreal* mdirs,
					    CFreal* In,
					    CFreal* II,
					    CFreal* divQ,
					    CFreal* qx, CFreal* qy, CFreal* qz)
{    
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  __shared__ CFreal weight;
  weight = weightIn;
  __syncthreads();

  if (cellID < nbCells) {
    
    // allocate the cell entity
    const CFuint iCell = abs(advanceOrder[d*nbCells+cellID]);
    CFreal inDirDotnANeg = 0.;
    CFreal Ic            = 0.;
    CFreal dirDotnAPos   = 0.;

    const CFuint nbFaces = nbFacesInCell[iCell];
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID = cellFaces[iFace*nbCells + iCell];
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = dotProdInFace[faceID]*factor;
      
      if (dirDotNA >= 0.){
	dirDotnAPos += dirDotNA;
      }
      else {
	const CFint fcellID = faceCell[faceID*2]; 
        const CFint neighborCellID = (fcellID == iCell) ? faceCell[faceID*2+1] : fcellID;
	const CFreal source = (neighborCellID >=0) ? In[neighborCellID] : fieldSource[iCell];
        inDirDotnANeg += source*dirDotNA;
      }
    } 
    In[iCell] = (fieldAbSrcV[iCell] - inDirDotnANeg)/(fieldAbV[iCell] + dirDotnAPos);
    Ic = In[iCell];
    
    qx[iCell] += Ic*mdirs[d*3]*weight;
    qy[iCell] += Ic*mdirs[d*3+1]*weight;
    qz[iCell] += Ic*mdirs[d*3+2]*weight;
    
    CFreal inDirDotnA = inDirDotnANeg;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID = cellFaces[iFace*nbCells + iCell];
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = dotProdInFace[faceID]*factor;
      if (dirDotNA > 0.) {
	inDirDotnA += In[iCell]*dirDotNA;
      }
    }
    
    divQ[iCell] += inDirDotnA*weight;
    II[iCell] += Ic*weight;
  }  
}
      
//////////////////////////////////////////////////////////////////////////////

RadiativeTransferFVDOMCUDA::RadiativeTransferFVDOMCUDA(const std::string& name) :
  RadiativeTransferFVDOM(name),
  m_faceCell(),
  m_nbFacesInCell(),
  m_InDir(),
  m_qxDir(),
  m_qyDir(),
  m_qzDir(),
  m_divqDir()  
{
  addConfigOptionsTo(this);
  
  m_qAlgoName = "Atomic";
  setParameter("QAlgoName", &m_qAlgoName);
}
      
//////////////////////////////////////////////////////////////////////////////

RadiativeTransferFVDOMCUDA::~RadiativeTransferFVDOMCUDA()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOMCUDA::defineConfigOptions(Config::OptionList& options)
{  
  options.addConfigOption< string >
    ("QAlgoName",
     "Name of the algorithm (Atomic, BigMem) to use for computing Q on the GPU.");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOMCUDA::setup()
{
  CFAUTOTRACE;
    
  RadiativeTransferFVDOM::setup();
  
  CFLog(VERBOSE, "RadiativeTransferFVDOMCUDA::setup() => START\n");
  
  // store invariant data on GPU
  SafePtr<ConnectivityTable<CFuint> > cellFaces = 
    MeshDataStack::getActive()->getConnectivity("cellFaces");
  cellFaces->getPtr()->put();
  m_dirs.put();
  m_weight.put();
  
  /*CFLog(INFO, "dirs: ");
  for (CFuint i = 0; i < m_dirs.size(); ++i) {
   CFLog(INFO, m_dirs[i] << " "); 
  }

  CFLog(INFO, "weights: ");
  for (CFuint i = 0; i < m_weight.size(); ++i) {
   CFLog(INFO, m_weight[i] << " ");  
  }*/

  m_fieldSource.put();  // to be removed 
  m_fieldAbsor.put();  // to be removed
  m_fieldAbSrcV.put(); // to be removed
  m_fieldAbV.put();  // to be removed
  if (!m_loopOverBins) {
    m_In.put(); // to be removed
    m_II.put(); // to be removed
  }
  m_opacities.put();
  m_radSource.put();
  m_Ttable.put();
  m_Ptable.put();
  m_advanceOrder.put(); // this can be a very big storage
  
  // AL: redundant
  const CFuint nbCells = m_fieldSource.size();
  if(m_useExponentialMethod){
    for (CFuint i=0;i<nbCells;++i) {
      m_fieldSource[i] = m_fieldAbsor[i] = m_In[i] = m_II[i] = 0;
    }
  }
  else {
    for (CFuint i=0;i<nbCells;++i) {
      m_fieldSource[i] = m_fieldAbSrcV[i] = m_fieldAbV[i] = m_In[i] = m_II[i] = 0;
    }
  }
  
  const CFuint totalNbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  m_faceCell.resize(totalNbFaces*2);
  m_faceCell = -1;
    
  m_nbFacesInCell.resize(nbCells);
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    const CFuint nbFaces = cellFaces->nbCols(iCell);
    m_nbFacesInCell[iCell] = nbFaces;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID2 = (*cellFaces)(iCell, iFace)*2;
      if (m_faceCell[faceID2] == -1) {
	m_faceCell[faceID2] = iCell;
      }
      else {
	m_faceCell[faceID2+1] = iCell;
      }
    }
  }
  
  m_faceCell.put();
  m_nbFacesInCell.put();
  
  if (m_loopOverBins) {  
    // AL: check if this creates a memory leak at the exit when !m_loopOverBins
    const CFuint nbDirs = m_advanceOrder.size()/nbCells;
    m_InDir.resize(nbCells*nbDirs);
    
    if (m_qAlgoName == "BigMem") {
      m_qxDir.resize(nbCells*nbDirs);
      m_qyDir.resize(nbCells*nbDirs);
      m_qzDir.resize(nbCells*nbDirs);
      m_divqDir.resize(nbCells*nbDirs);
      
      for (CFuint i = 0; i < m_InDir.size(); ++i) {
	m_qxDir[i] = m_qyDir[i] = m_qzDir[i] = m_divqDir[i] = 0.;
      }
      
      m_qxDir.put();
      m_qyDir.put();
      m_qzDir.put();
      m_divqDir.put();
    }
    
    for (CFuint i = 0; i < m_InDir.size(); ++i) {
      m_InDir[i] = 0.;
    }
    m_InDir.put();
  }
  
  CFLog(VERBOSE, "RadiativeTransferFVDOMCUDA::setup() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOMCUDA::unsetup()
{
  CFAUTOTRACE;
  
  RadiativeTransferFVDOM::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////
 
void RadiativeTransferFVDOMCUDA::loopOverDirs(const CFuint startBin, 
					      const CFuint endBin, 
					      const CFuint startDir,
					      const CFuint endDir)
{
  CFLog(VERBOSE, "RadiativeTransferFVDOMCUDA::loopOverDirs() => START\n");
  
  SafePtr<ConnectivityTable<CFuint> > cellFaces = 
    MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();
  
  const CFuint nbCells = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  states.getGlobalArray()->put();
  isOutward.getLocalArray()->put(); 
  volumes.getLocalArray()->put(); 
  
  const CFuint blocksPerGrid = 
    CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
  const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
    
  for (CFuint d = startDir; d < endDir; ++d) {
    CFLog(INFO, "( dir: " << d << " ), ( bin: ");
    const CFuint bStart = (d != startDir) ? 0 : startBin;
    const CFuint bEnd   = (d != m_startEndDir.second) ? m_nbBins : endBin;
    
    // precompute dot products for all faces and directions (a part from the sign)
    computeDotProdInFace(d, m_dotProdInFace);
    m_dotProdInFace.put();
    
    for (CFuint ib = startBin; ib < endBin; ++ib) {
      CFLog(INFO, ib << " ");
      
      // precompute the radiation properties for all cells
      getFieldOpacitiesKernel<<<blocksPerGrid,nThreads>>>
	(m_useExponentialMethod, 
	 m_TID, m_PID, m_nbTemp, m_nbPress, m_nbBins,
	 ib, nbEqs, nbCells, m_Ttable.ptrDev(), m_Ptable.ptrDev(), 
	 states.getGlobalArray()->ptrDev(),
	 volumes.getLocalArray()->ptrDev(),
	 m_opacities.ptrDev(),
	 m_radSource.ptrDev(),
	 m_fieldSource.ptrDev(),
	 m_fieldAbsor.ptrDev(),
	 m_fieldAbSrcV.ptrDev(),
	 m_fieldAbV.ptrDev());  
      
      m_fieldSource.get();
      m_fieldAbsor.get();
      m_fieldAbSrcV.get();
      m_fieldAbV.get();
      
      // compute the radiative heat flux: unfortunately, this cannot be done
      // on the GPU since the computeQ algorithm is serial (the result on one 
      // cell depends on previously computed cells) and must follow the predefined 
      // advance ordering
      (m_useExponentialMethod) ? 
	computeQExponential(ib,startDir,d) : computeQNoExponential(ib,startDir,d);
    }
    CFLog(INFO, ")\n");
  }
  
  CFLog(VERBOSE, "RadiativeTransferFVDOMCUDA::loopOverDirs() =>END\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void RadiativeTransferFVDOMCUDA::loopOverBins(const CFuint startBin, 
					      const CFuint endBin, 
					      const CFuint startDir,
					      const CFuint endDir)
{
  CFLog(VERBOSE, "RadiativeTransferFVDOMCUDA::loopOverBins() => START\n");
  
  SafePtr<ConnectivityTable<CFuint> > cellFaces = 
    MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  const CFuint nbCells = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // store useful data on the GPU
  states.getGlobalArray()->put();
  isOutward.getLocalArray()->put(); 
  volumes.getLocalArray()->put(); 
  socket_normals.getDataHandle().getLocalArray()->put();
  
  const CFuint blocksPerGrid = 
    CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
  const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
  
  // if more than one iteration is needed, the initialization has to be done here
  // for the moment it is done the RadiativeTransferFVDOM::setup();
  // divq = 0.; qx = 0.; qy = 0.; qz = 0.; 
  socket_divq.getDataHandle().getLocalArray()->put();
  socket_qx.getDataHandle().getLocalArray()->put();
  socket_qy.getDataHandle().getLocalArray()->put();
  socket_qz.getDataHandle().getLocalArray()->put();

  Stopwatch<WallTime> stp;
  stp.start();
  
  for(CFuint ib = startBin; ib < endBin; ++ib) {
    CFLog(INFO, "( bin: " << ib << " ), ( dir: ");
    
    if (m_oldAlgo) { 
      // precompute the radiation properties for all cells
      getFieldOpacitiesKernel<<<blocksPerGrid,nThreads>>>
	(m_useExponentialMethod, 
	 m_TID, m_PID, m_nbTemp, m_nbPress, m_nbBins,
	 ib, nbEqs, nbCells, m_Ttable.ptrDev(), m_Ptable.ptrDev(), 
	 states.getGlobalArray()->ptrDev(),
	 volumes.getLocalArray()->ptrDev(),
	 m_opacities.ptrDev(),
	 m_radSource.ptrDev(),
	 m_fieldSource.ptrDev(),
	 m_fieldAbsor.ptrDev(),
	 m_fieldAbSrcV.ptrDev(),
	 m_fieldAbV.ptrDev());  
    }
    
    /* m_fieldSource.get();
       m_fieldAbsor.get();
       m_fieldAbSrcV.get();
       m_fieldAbV.get();
       
       for (CFuint k = 0; k < m_fieldSource.size(); ++k)
       CFLog(INFO, "m_fieldSource[" <<k << "] => (" << m_fieldSource[k] << "\n");
       for (CFuint k = 0; k < m_fieldAbsor.size(); ++k)
       CFLog(INFO, "m_fieldAbsor[" <<k << "] => (" << m_fieldAbsor[k] << "\n");
       for (CFuint k = 0; k < m_fieldAbSrcV.size(); ++k)
       CFLog(INFO, "m_fieldAbSrcV[" <<k << "] => (" << m_fieldAbSrcV[k] << "\n");
       for (CFuint k = 0; k < m_fieldAbV.size(); ++k)
       CFLog(INFO, "m_fieldAbV[" <<k << "] => (" << m_fieldAbV[k] << "\n");
       
       AL: this works perfectly
       exit(1);*/
    
    const CFuint dStart = (ib != startBin) ? 0 : startDir;
    const CFuint dEnd = (ib != m_startEndBin.second)? m_nbDirs : endDir;
    const CFuint nbDirs = dEnd - dStart;
    const CFuint blocksPerDir = 
      CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbDirs);
    const CFuint nThreadsDir = 
      std::min((CFuint)CudaEnv::CudaDeviceManager::getInstance().getNThreads(), nbDirs);
    
    if (m_useExponentialMethod) {
      if (m_qAlgoName == "Atomic") {
	// AL: this algorithm is my preferred choice, since needs much less memory
	//     but uses atomics
	computeQKernelExponentialDirAtomic<<<blocksPerDir,nThreadsDir>>> 
	  (dStart, dEnd, nbCells,
	   cellFaces->getPtr()->ptrDev(),
	   m_faceCell.ptrDev(),
	   m_nbFacesInCell.ptrDev(),
	   isOutward.getLocalArray()->ptrDev(),
	   m_advanceOrder.ptrDev(),
	   m_weight.ptrDev(),
	   volumes.getLocalArray()->ptrDev(),
	   m_fieldSource.ptrDev(),
	   m_fieldAbsor.ptrDev(),
	   normals.getLocalArray()->ptrDev(),
	   m_dirs.ptrDev(),
	   m_InDir.ptrDev(),
	   divQ.getLocalArray()->ptrDev(),
	   qx.getLocalArray()->ptrDev(),
	   qy.getLocalArray()->ptrDev(),
	   qz.getLocalArray()->ptrDev());
	
	// this will allow for printing while exiting right after 
	// cudaDeviceReset(); exit(1);
      }
      else if (m_qAlgoName == "BigMem") {
	// AL: this algorithm is similar to the other one but can avoid using 
	//     atomics by storing much more memory
	
	// both algorithms have similar performance problems...
	computeQKernelExponentialDirBigMem<<<blocksPerDir,nThreadsDir>>> 
	  (dStart, dEnd, nbCells,
	   cellFaces->getPtr()->ptrDev(),
	   m_faceCell.ptrDev(),
	   m_nbFacesInCell.ptrDev(),
	   isOutward.getLocalArray()->ptrDev(),
	   m_advanceOrder.ptrDev(),
	   m_weight.ptrDev(),
	   volumes.getLocalArray()->ptrDev(),
	   m_fieldSource.ptrDev(),
	   m_fieldAbsor.ptrDev(),
	   normals.getLocalArray()->ptrDev(),
	   m_dirs.ptrDev(),
	   m_InDir.ptrDev(),
	   m_divqDir.ptrDev(),
	   m_qxDir.ptrDev(),
	   m_qyDir.ptrDev(),
	   m_qzDir.ptrDev());
	
	// this will allow for printing while exiting right after 
	// cudaDeviceReset(); exit(1);
      }
      
      /*const CFuint nbKernelCalls = std::ceil((CFreal)nbCells/(CFreal)SHARED_ARRAY_SIZE); 
	
	CFLog(VERBOSE, "RadiativeTransferFVDOMCUDA::loopOverBins() => nbBlocks["
	<< blocksPerDir <<"], nbThreads[" << nThreadsDir 
	<< "], nbKernelCalls = " << nbKernelCalls <<"\n");
	
	CFuint startCellID = 0;
	CFuint endCellID = std::min(nbCells, (CFuint) SHARED_ARRAY_SIZE);
	for (CFuint k = 0; k < nbKernelCalls; ++k) {
	computeQKernelExponentialDir<<<blocksPerDir,nThreadsDir>>> 
	(dStart, dEnd, startCellID, endCellID, nbCells,
	cellFaces->getPtr()->ptrDev(),
	m_faceCell.ptrDev(),
	m_nbFacesInCell.ptrDev(),
	isOutward.getLocalArray()->ptrDev(),
	m_advanceOrder.ptrDev(),
	m_weight.ptrDev(),
	volumes.getLocalArray()->ptrDev(),
	m_fieldSource.ptrDev(),
	m_fieldAbsor.ptrDev(),
	normals.getLocalArray()->ptrDev(),
	m_dirs.ptrDev(),
	divQ.getLocalArray()->ptrDev(),
	qx.getLocalArray()->ptrDev(),
	qy.getLocalArray()->ptrDev(),
	qz.getLocalArray()->ptrDev());
	
	startCellID = endCellID;
	endCellID = std::min(nbCells, endCellID + SHARED_ARRAY_SIZE);
	}*/
    }
    else {
      computeQKernelNoExponentialDir<<<blocksPerDir,nThreadsDir>>> 
	(dStart, dEnd, nbCells,
	 cellFaces->getPtr()->ptrDev(),
	 m_faceCell.ptrDev(),
	 m_nbFacesInCell.ptrDev(),
	 isOutward.getLocalArray()->ptrDev(),
	 m_advanceOrder.ptrDev(),
	 m_weight.ptrDev(),
	 volumes.getLocalArray()->ptrDev(),
	 m_fieldSource.ptrDev(),
	 m_fieldAbSrcV.ptrDev(),
	 m_fieldAbV.ptrDev(),
	 normals.getLocalArray()->ptrDev(),
	 m_dirs.ptrDev(),
	 m_InDir.ptrDev(), 
	 divQ.getLocalArray()->ptrDev(),
	 qx.getLocalArray()->ptrDev(),
	 qy.getLocalArray()->ptrDev(),
	 qz.getLocalArray()->ptrDev());
    }
    
    for (CFuint d = dStart; d < dEnd; ++d) {
      CFLog(INFO, d << " ");
      //(m_useExponentialMethod) ? 
      //computeQExponential(ib,d) : computeQNoExponential(ib,d);
    }
    CFLog(INFO, ")\n");
  }
  
  CFLog(INFO, "RadiativeTransferFVDOMCUDA::loopOverBins() => loop took " << stp.read() << "s \n");
  
  stp.start();
  
  if (m_qAlgoName == "Atomic") {
    socket_divq.getDataHandle().getLocalArray()->get();
    socket_qx.getDataHandle().getLocalArray()->get();
    socket_qy.getDataHandle().getLocalArray()->get();
    socket_qz.getDataHandle().getLocalArray()->get();
  }
  else if (m_qAlgoName == "BigMem") {
    m_qxDir.get(); 
    m_qyDir.get();
    m_qzDir.get();
    m_divqDir.get();
   
    // this is only working on single GPU
    // neds to be fixed for multi-GPU  
    const CFuint dStart = 0; 
    const CFuint dEnd = m_nbDirs; 
    for (CFuint c = 0; c < nbCells; ++c) {
      qx[c] = qy[c] = qz[c] = divQ[c] = 0.;
      const CFuint startc = c*m_nbDirs;
      for (CFuint d = dStart; d < dEnd; ++d) {
	const CFuint startd = startc+d;
	qx[c]   += m_qxDir[startd];
        qy[c]   += m_qyDir[startd];
        qz[c]   += m_qzDir[startd];
	divQ[c] += m_divqDir[startd];
      }
    }
  }
  
  CFLog(INFO, "RadiativeTransferFVDOMCUDA::loopOverBins() => GPU-CPU transfer took " << stp.read() << "s \n");
  
  /* for (CFuint k = 0; k < socket_divq.getDataHandle().size(); ++k) {
     CFLog(INFO, "divQ[" <<k << "] => (" << socket_divq.getDataHandle()[k] << "\n");
     }
     for (CFuint k = 0; k < socket_qx.getDataHandle().size(); ++k) {
     CFLog(INFO, "qx[" <<k << "] => (" << socket_qx.getDataHandle()[k] << "\n");
     }
     for (CFuint k = 0; k < socket_qy.getDataHandle().size(); ++k) {
     CFLog(INFO, "qy[" <<k << "] => (" << socket_qy.getDataHandle()[k] << "\n");
     }
     for (CFuint k = 0; k < socket_qz.getDataHandle().size(); ++k) {
     CFLog(INFO, "qz[" <<k << "] => (" << socket_qz.getDataHandle()[k] << "\n");
     }
     exit(1);*/// 
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

