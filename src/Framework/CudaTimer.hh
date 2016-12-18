// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_CudaTimer_hh
#define COOLFluiD_CudaEnv_CudaTimer_hh 

//////////////////////////////////////////////////////////////////////////////

#include "CudaDeviceManager.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
          
//////////////////////////////////////////////////////////////////////////////

/// This class implements a stopwatch for CUDA
///
/// @author Andrea Lani
class CudaTimer {
public:
  /// member returning a reference to a statically allocated object:
  /// makes this class a singleton  
  static CudaTimer& getInstance() {static CudaTimer cTimer; return cTimer;}  
  
  /// start timing  
  void start() 
  { 
   CUDA_CHECK(cudaEventCreate( &m_start ));
   CUDA_CHECK(cudaEventCreate( &m_stop ));
   CUDA_CHECK(cudaEventRecord( m_start, 0 ));
  }
  
  /// tell elapsed time
  double elapsed()
  { 
    CUDA_CHECK(cudaEventRecord( m_stop, 0 ));
    CUDA_CHECK(cudaEventSynchronize( m_stop ));
    float elapsed = 0.; 
    CUDA_CHECK(cudaEventElapsedTime( &elapsed, m_start, m_stop ));
    return (double)elapsed/1000.;
  }
  
private: // members
  
  /// constructor
  CudaTimer() {}
  
  /// destructor
  ~CudaTimer() {}
  
private: // data
  
  /// start event
  cudaEvent_t  m_start;

  /// stop event
  cudaEvent_t  m_stop;
};

//////////////////////////////////////////////////////////////////////////////

} // end namespace CudaEnv 
  
} // end  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
