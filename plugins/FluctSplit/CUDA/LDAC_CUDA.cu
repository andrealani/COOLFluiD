#include "FluctSplit/CUDA/LDAC_CUDA.hh"
#include "Framework/CudaDeviceManager.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////
      
__global__ void addToKPlus(int msize, double* a,  double* b, double* c, double* d) 
{
  int tid = threadIdx.x; // + blockIdx.x * blockDim.x;
  if (tid < msize) {
    d[tid] = a[tid] + b[tid] + c[tid];
    // tid += blockDim.x*gridDim.x;`
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void addToKplusCUDA(int msize, 
		    double* dev_a,  double* dev_b, double* dev_c, double* dev_d, 
		    double* k0,  double* k1,  double* k2, double* sumKplus)
{
  using namespace std;
  
  CUDA_CHECK(cudaMemcpy(dev_a, k0, msize*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(dev_b, k1, msize*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(dev_c, k2, msize*sizeof(double), cudaMemcpyHostToDevice));
  
  addToKPlus<<<1,msize>>>(msize, dev_a, dev_b, dev_c, dev_d);
  
  CUDA_CHECK(cudaMemcpy(sumKplus, dev_d, msize*sizeof(double), cudaMemcpyDeviceToHost));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
