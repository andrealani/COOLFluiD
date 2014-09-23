#ifndef COOLFluiD_Config_ConfigOptionPtr_hh
#define COOLFluiD_Config_ConfigOptionPtr_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Config {
          
//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores configurable options options data to be exchanged with 
 * a device
 *
 * @author Andrea Lani
 *
 */
template <typename T, typename P = NOTYPE, DeviceType DT = CPU> 
class ConfigOptionPtr {
public:
  
  /// Constructor
  ConfigOptionPtr(Common::SafePtr<T> obj) : m_dco()
  {
    obj->copyConfigOptions(&m_dco);
  }
  
  /// Destructor
  ~ConfigOptionPtr() {}
  
  ///get the raw pointer
  typename T::template DeviceConfigOptions<P>* getPtr() {return &m_dco;} 
  
private:
  
  /// pointer to the raw data
  typename T::template DeviceConfigOptions<P> m_dco;
};
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores configurable options options data to be exchanged with 
 * a device. Partial specialization for GPU.
 *
 * @author Andrea Lani
 *
 */
#ifdef CF_HAVE_CUDA
template <typename T, typename P> 
class ConfigOptionPtr<T, P, GPU> {
public:
  
  /// Constructor
  ConfigOptionPtr(Common::SafePtr<T> obj)
  {
    // allocate data pointer on device
    CudaEnv::allocDev(m_dco, 1); 
    // copy the configurable data to the device
    obj->copyConfigOptionsToDevice(m_dco);
  }
  
  /// Destructor
  ~ConfigOptionPtr() 
  {
    // deallocate data pointer on device
    CudaEnv::free(m_dco);
  }
  
  ///get the raw pointer
  typename T::template DeviceConfigOptions<P>* getPtr() const {return m_dco;} 
  
private:
  
  /// pointer to the raw data
  typename T::template DeviceConfigOptions<P>* m_dco;
};
#endif
    
    //////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_ConfigOptionPtr_hh
