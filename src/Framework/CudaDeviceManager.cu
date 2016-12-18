// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/CudaDeviceManager.hh"
#include "Common/CFLog.hh"
#include "Common/PE.hh"
	
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CudaEnv {
          
//////////////////////////////////////////////////////////////////////////////

void CudaDeviceManager::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("NbThreadsPerBlock", "Number if threads per block");
  options.addConfigOption< CFuint >("NbBlocks", "Number of blocks");
}
    
//////////////////////////////////////////////////////////////////////////////

CudaDeviceManager::CudaDeviceManager() : Config::ConfigObject("CudaDeviceManager"), 
					 Common::NonCopyable<CudaDeviceManager>()
{
  addConfigOptionsTo(this);

  NTHREADS_PER_BLOCK = 0;
  setParameter("NbThreadsPerBlock", &NTHREADS_PER_BLOCK);
  
  NBLOCKS = 0;
  setParameter("NbBlocks", &NBLOCKS);
}
    
//////////////////////////////////////////////////////////////////////////////
 
void CudaDeviceManager::configure ( Config::ConfigArgs& args )
{
  CFLog(VERBOSE, "CudaDeviceManager::configure() START\n");
  
  ConfigObject::configure(args);
  
  CFLog(VERBOSE, "CudaDeviceManager::configure() => cudaGetDeviceCount()\n");
  
  // the first time getInstance() is called, the device is initialized
  // for the moment we assume one device
  int count = 0; 
  int dev   = -1;
  int rank = Common::PE::GetPE().GetRank("Default");
  cudaGetDeviceCount(&count);
 
  const CFuint nbProc = Common::PE::GetPE().GetProcessorCount("Default"); 
  if (nbProc > 1) {
    dev = rank % count;
    cfassert(dev >= 0); 
    cudaSetDevice(dev);
    //    std::cout << "P" << rank << " => device/count => " << dev << "/" << count << std::endl;
    // std::cout << "CudaDeviceManager::configure() => infinite loop ...\n"; 
    // for (;;) {}
  } else {
     /* the code below works for serial GPU simulations */
    cudaGetDevice(&dev);
    cfassert(dev >= 0);  
    cudaSetDevice(dev);
  }
  
  CFLog(VERBOSE, "CudaDeviceManager::configure() => cudaGetDeviceProperties()\n");
  cudaGetDeviceProperties(&m_prop, dev); 
  
  CFLog(VERBOSE, "CudaDeviceManager::configure() => printProperties()\n");
  printProperties(dev);
  
  NTHREADS_PER_BLOCK = (NTHREADS_PER_BLOCK > 0) ? NTHREADS_PER_BLOCK :  m_prop.maxThreadsPerBlock;
  cfassert(NTHREADS_PER_BLOCK > 0);
  
  NBLOCKS = (NBLOCKS > 0) ? NBLOCKS : m_prop.maxGridSize[0];
  cfassert(NBLOCKS > 0);
  
  CFLog(INFO, "CudaDeviceManager::configure() => NTHREADS_PER_BLOCK = " << NTHREADS_PER_BLOCK 
	<< ", NBLOCKS = " << NBLOCKS << " END\n");
  
  CFLog(VERBOSE, "CudaDeviceManager::configure() END\n");
}
    
//////////////////////////////////////////////////////////////////////////////
 
void CudaDeviceManager::printProperties(int dev)
{
  using namespace std;
  
  CFLog(INFO, "##### CudaDeviceManager::printProperties() for device [" << dev << "] #####\n");
  CFLog(INFO, "name = " << m_prop.name << "\n");
  CFLog(INFO, "capability = " << m_prop.major << "." << m_prop.minor << "\n"); 
  CFLog(INFO, "clock rate = " << m_prop.clockRate <<"\n"); 
  CFLog(INFO, "total global mem = " << m_prop.totalGlobalMem <<"\n"); 
  CFLog(INFO, "total constant mem = " << m_prop.totalConstMem <<"\n");
  CFLog(INFO, "overlap execution and transfer = "); 
  if (m_prop.deviceOverlap) {
    CFLog(INFO, "ENABLED\n");
  }
  else {
    CFLog(INFO, "DISABLED\n");
  }
  CFLog(INFO, "can map host memory = "); 
  if (m_prop.canMapHostMemory) {
    CFLog(INFO, "ENABLED\n");
  }
  else {
    CFLog(INFO, "DISABLED\n");
  }
  CFLog(INFO, "texture alignment = " << m_prop.textureAlignment << "\n");
  CFLog(INFO, "multiprocessor count = " << m_prop.multiProcessorCount << "\n");
  CFLog(INFO, "shared mem per block = " << m_prop.sharedMemPerBlock << "\n");
  CFLog(INFO, "registers per block = " << m_prop.regsPerBlock << "\n");
  CFLog(INFO, "threads in warp = " << m_prop.warpSize << "\n");
  CFLog(INFO, "max threads per block = " << m_prop.maxThreadsPerBlock << "\n");
  CFLog(INFO, "max threads dimensions = " << m_prop.maxThreadsDim[0] << " " 
       <<  m_prop.maxThreadsDim[1] << " " <<  m_prop.maxThreadsDim[2] << "\n");
  CFLog(INFO, "max grid dimensions = " << m_prop.maxGridSize[0] << " " 
       <<  m_prop.maxGridSize[1] << " " <<  m_prop.maxGridSize[2] << "\n");
  CFLog(INFO,  "############################################################\n\n"); 
}

//////////////////////////////////////////////////////////////////////////////

} // end namespace CudaEnv 
  
} // end  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
