#include <paralution.hpp>
#include "ParalutionMatrix.hh"

#include "Framework/BlockAccumulatorBaseCUDA.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/CudaDeviceManager.hh"
#include "Common/CUDA/CFVec.hh"
#include "FiniteVolume/CellData.hh"
#include "Common/CUDA/CudaEnv.hh"
#include "Framework/CudaDeviceManager.hh"


using namespace COOLFluiD::Framework;
//using namespace COOLFluiD::Common;
//using namespace COOLFluiD::Config;
using namespace std;

namespace COOLFluiD {

    namespace Paralution {

/*
__global__ void addDiagBlockGPU(CFreal *val, CFint *rowoff, CFint *col, CFreal *accDevPtr, CFuint *nbs, CFint *IDs)
{
  const int nbj = threadIdx.x; //+ blockIdx.x*blockDim.x;
  const int nbi = threadIdx.y; // + blockIdx.x*blockDim.x;
  //printf("nbi %d \t nbj %d \n", nbi, nbj);
  //if(cellID!=0){printf("inside addDiagBlockGPU %d \n", cellID);}
  //if (cellID < nbCells) { 
    const CFint nb = nbs[0];

    CFuint RowPositionDiag = rowoff[IDs[0]*nb];           //In this case we are looking for the diagonal block
    CFuint RowPositionPlusOneDiag = rowoff[IDs[0]*nb + 1];
    CFuint mmDiag = (RowPositionPlusOneDiag-RowPositionDiag)/nb;
    CFuint IndexCSRDiag = -1;

          IndexCSRDiag = RowPositionDiag; //+mii*nb;

    //for (CFint nbj=0; nbj<nb; nbj++){
      //for (CFint nbi=0; nbi<nb; nbi++){
        val[IndexCSRDiag+nbi*nb*mmDiag+nbj] += accDevPtr[nbi*nb + nbj];
     // }
    //}
  //}
} 
*/

__global__ void addDiagBlockGPU(CFreal *val, CFint *rowoff, CFreal *accDevPtr, CFuint nb, CFuint startCellID, CFint nbCells)
{
  const int localID = threadIdx.x + blockIdx.x*blockDim.x;
  const int cellID = localID + startCellID;
  if (cellID < nbCells) { 
    CFuint RowPositionDiag = rowoff[cellID*nb];
    CFuint RowPositionPlusOneDiag = rowoff[cellID*nb + 1];
    CFuint mmDiag = (RowPositionPlusOneDiag-RowPositionDiag)/nb;
    for(CFint nbi=0; nbi<nb; nbi++){
      for(CFint nbj=0; nbj<nb; nbj++){
          val[RowPositionDiag+nbi*nb*mmDiag+nbj] += accDevPtr[localID*nb*nb + nbi*nb + nbj];
      }
    }
  }
}

////////////////////////////////////////////////////////////////////

void ParalutionMatrix::resetToZeroEntriesGPU(){
   //std::cout << "resetToZeroEntriesGPU \n";
   //CFuint nbThreadPerBlock = 64;
   //CFuint nbBlocks = _size/64;
cudaMemset(_valDev, 0.0, _size*sizeof(CFreal));
std::fill_n(_diagAcc, diagAccSize, 0);
   //ZeroEntriesGPU<<< nbBlocks,nbThreadPerBlock >>> (_valDev, _size);
}

////////////////////////////////////////////////////////////////////

void ParalutionMatrix::addValuesGPU(const Framework::BlockAccumulator& acc)
{
   CFuint nb = acc.getNB();
   CFreal* accPtr = const_cast<Framework::BlockAccumulator&>(acc).getPtr();
   CFint IDs = const_cast<std::vector<CFint>&>(acc.getIN())[0]; //Array storing cellID

   CFuint cellIndex = IDs*nb*nb;
   for(CFint nbi=0; nbi<nb; nbi++){
     for(CFint nbj=0; nbj<nb; nbj++){
        _diagAcc[cellIndex + nbi*nb + nbj] += accPtr[nbi*nb + nbj];
     }
   } 
}  


void ParalutionMatrix::updateDiagBlocks(CFuint nbCells, CFuint nbEqs)
{
   //Copy the array to GPU

   CFuint startCellID = 0;
   for (CFuint s = 0; s < _sizeb; ++s) {
      CudaEnv::copyHost2Dev(_diagAccDev, &_diagAcc[startCellID*nbEqs*nbEqs], _nThreads*_nbKernelBlocks*nbEqs*nbEqs);
      addDiagBlockGPU <<<_nbKernelBlocks,_nThreads>>> (_valDev, _rowoffDev, _diagAccDev, nbEqs, startCellID, nbCells);
 
      startCellID += _nThreads*_nbKernelBlocks; //m_nbCellsInKernel[s];
   }

}

/*
void ParalutionMatrix::addValuesGPU(const Framework::BlockAccumulator& acc)
{
   //const CFuint blocksPerGrid = CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
   //const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();

   //Always m=n=1!! -> More efficient algorithm

   //Read the data from the BlockAccumulator
    CFuint size = acc.size();
   // CFuint m = acc.getM();
   // CFuint n = acc.getN();
    CFuint nb = acc.getNB();
   CFint nbID = const_cast<std::vector<CFint>&>(acc.getIM())[0]; //Array storing indexes of neigbours
   CFint IDs = const_cast<std::vector<CFint>&>(acc.getIN())[0]; //Array storing cellID
   CFreal* accPtr = const_cast<Framework::BlockAccumulator&>(acc).getPtr();

   //Create the device pointers
   CFreal* accDevPtr;
   //CFint* nbIDDev;
   CFint* IDsDev;
   CFuint* nbDev;

   //Allocate
   CudaEnv::allocDev(accDevPtr, size);
   CudaEnv::allocDev(nbDev, 1);
   CudaEnv::allocDev(IDsDev, 1);

   //Copy the data
   CudaEnv::copyHost2Dev(accDevPtr, accPtr, size);
   CudaEnv::copyHost2Dev(nbDev, &nb, 1);
   CudaEnv::copyHost2Dev(IDsDev, &IDs, 1);
   //CudaEnv::copyHost2Dev(mDev, &m, 1);
   //CudaEnv::copyHost2Dev(nDev, &n, 1);

   //Call kernel to add values
   //if(IDs != nbID){printf("IDs %d \t nbIDs %d \n", IDs, nbID);}
   dim3 threads(nb, nb);
   addDiagBlockGPU<<<1,threads>>> (_valDev, _rowoffDev, _colDev, accDevPtr, nbDev, IDsDev);

   // Deallocate memory GPU
   CudaEnv::free(accDevPtr);
   CudaEnv::free(nbDev);
   CudaEnv::free(IDsDev);
   //CudaEnv::free(mDev);
   //CudaEnv::free(nDev);
}
*/


    }   // namespace COOLFluiD

}  // namespace Paralution
