// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CudaEnv_CudaDeviceManager_hh
#define COOLFluiD_CudaEnv_CudaDeviceManager_hh 

//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>

#include "Config/ConfigObject.hh"
#include "Common/NonCopyable.hh"
#include "Common/CUDA/CudaEnv.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
          
//////////////////////////////////////////////////////////////////////////////

/// This class provides and stores info about the device
///
/// @author Andrea Lani
class CudaDeviceManager : public Config::ConfigObject,
			  public Common::NonCopyable<CudaDeviceManager> {
  
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /// member returning a reference to a statically allocated object:
  /// makes this class a singleton  
  static CudaDeviceManager& getInstance() 
  {
    static CudaDeviceManager cDeviceManager; return cDeviceManager;
  }  
  
  /// print the device properties  
  void printProperties(int dev);
  
  /// get the number of threads
  int getNThreads() const {return NTHREADS_PER_BLOCK;}
  
  /// get the number of blocks
  int getNBlocks() const {return NBLOCKS;}
  
  /// get the number of blocks per grid
  int getBlocksPerGrid(int N) 
  {
    return std::min(NBLOCKS, (N + NTHREADS_PER_BLOCK-1)/NTHREADS_PER_BLOCK);
  }
  
  /// Configure alluser-defined  parameters in this object
  void configure ( Config::ConfigArgs& args );
  
private: // members

  /// constructor
  CudaDeviceManager(); 
  
  /// destructor
  ~CudaDeviceManager() {}
  
private: // data
  
  /// maximum number of threads per block
  CFuint NTHREADS_PER_BLOCK;
  
  /// maximum number of blocks
  CFuint NBLOCKS;
  
  /// device properties
  cudaDeviceProp m_prop;

};

//////////////////////////////////////////////////////////////////////////////

} // end namespace CudaEnv 
  
} // end  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
