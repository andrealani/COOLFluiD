// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_hh
#define COOLFluiD_CudaEnv_hh

//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cuda_runtime.h>

//////////////////////////////////////////////////////////////////////////////

// Error handling macro for CUDA functions
#ifndef DNDEBUG
#define CUDA_CHECK(call)						\
    if((call) != cudaSuccess) {						\
      cudaError_t err = cudaGetLastError();				\
      std::cerr << "CUDA error calling \""#call"\", code is " << err << std::endl; \
    abort(); }
#define cfassert(expression)      \
  {				    \
    bool test = (expression);	    \
    				    \
    if (!test) {							\
      std::cout << "CFASSERT failed at line = " << __LINE__ << ", file = " << __FILE__ << "\n"; \
      exit(-1); /* So that destructors are run. */			\
    }									\
  }
#else 
#define CUDA_CHECK(call) call
#define cfassert(expression)
#endif

//////////////////////////////////////////////////////////////////////////////
 
namespace COOLFluiD {
  
  namespace CudaEnv {
      
//////////////////////////////////////////////////////////////////////////////
     
 // allocation of host buffer in pinned memory
  template <typename T>
  static inline void allocHost(T*& vec, size_t n) 
  {
    CUDA_CHECK(cudaHostAlloc((void**)&vec,  n*sizeof(T), cudaHostAllocDefault));
  }
  
  // allocation of device buffer
  template <typename T>
  static inline void allocDev(T*& vec_d, size_t n) 
  {
    CUDA_CHECK(cudaMalloc((void**)&vec_d,  n*sizeof(T)));
  }
  
  // deallocation of host buffer in pinned memory
  template <typename T>
  static inline void freeHost(T*& vec) 
  {
    CUDA_CHECK(cudaFreeHost(vec));
  }
  
  // deallocation of buffers
  template <typename T>
  static inline void free(T*& vec) 
  {
    CUDA_CHECK(cudaFree(vec));
  }
  
  // synchronous copy from CPU to GPU
  template <typename T>
  static inline void copyHost2Dev(T* vec_d, T* vec, size_t n)
  {
    CUDA_CHECK(cudaMemcpy(vec_d, vec, n*sizeof(T),cudaMemcpyHostToDevice));
  }
  
  // synchronous copy from GPU to CPU
  template <typename T>
  static inline void copyDev2Host(T* vec, T* vec_d, size_t n)
  {
     CUDA_CHECK(cudaMemcpy(vec, vec_d, n*sizeof(T),cudaMemcpyDeviceToHost));
  }
  
  // asynchronous copy from CPU to GPU
  template <typename T>
  static inline void copyAsyncHost2Dev(T* vec_d, T* vec, size_t n, cudaStream_t& stream)
  {
    CUDA_CHECK(cudaMemcpyAsync(vec_d, vec, n*sizeof(T),cudaMemcpyHostToDevice, stream));
  }
  
  // asynchronous copy from GPU to CPU
  template <typename T>
  static inline void copyAsyncDev2Host(T* vec, T* vec_d, size_t n, cudaStream_t& stream)
  {
    CUDA_CHECK(cudaMemcpyAsync(vec, vec_d, n*sizeof(T),cudaMemcpyDeviceToHost, stream));
  }
  
  // create stream
  static void createStream(cudaStream_t* stream) { CUDA_CHECK( cudaStreamCreate(stream)); }
  
  // destroy stream
  static void destroyStream(cudaStream_t& stream) { CUDA_CHECK( cudaStreamDestroy(stream)); }
  
  // synchronize stream
  static void synchronizeStream(cudaStream_t& stream) { CUDA_CHECK( cudaStreamSynchronize(stream)); }
    
//////////////////////////////////////////////////////////////////////////////

} // end namespace CudaEnv 
  
} // end  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
