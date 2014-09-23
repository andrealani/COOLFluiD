// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_CudaLock_hh
#define COOLFluiD_CudaEnv_CudaLock_hh 

//////////////////////////////////////////////////////////////////////////////

#include "CudaDeviceManager.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
          
//////////////////////////////////////////////////////////////////////////////

/// This class implements a lock as explained in "CUDA by example" 
/// by J. Sanders and E. Kandrot
///
/// @author Andrea Lani
struct CudaLock {

  /// Constructor
  CudaLock (void) 
  {
    int state = 0;
    CudaEnv::allocDev(mutex, 1);
    CudaEnv::copyHost2Dev(mutex, &state, 1); 
  }
  
  /// Destructor
  ~CudaLock (void) 
  {
    CudaEnv::free(mutex);
  }
  
  /// Lock the given datastructure
  __device__ void lock (void) {while (atomicCAS(mutex, 0, 1) != 0);}
  
  /// Unlock and restore the value of the mutex to 0
  __device__ void unlock (void) {atomicExch(mutex, 0);}
  
  int* mutex;
};

//////////////////////////////////////////////////////////////////////////////
    
} // end namespace CudaEnv 
  
} // end  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
