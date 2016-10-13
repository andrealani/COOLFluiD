#ifndef COOLFluiD_Physics_MultiFluidMHD_ReacMFMHDTerm_hh
#define COOLFluiD_Physics_MultiFluidMHD_ReacMFMHDTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a source term for MultiFluidMHD
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez
 *
 */
class ReacMFMHDTerm : public Framework::BaseTerm {
public:


#ifdef CF_HAVE_CUDA
    
   //Nested class defining local options
   template <typename P = NOTYPE >  
   class DeviceConfigOptions{
   public:
       //Constructor
       HOST_DEVICE DeviceConfigOptions() {}
       //Destructor
       HOST_DEVICE virtual ~DeviceConfigOptions() {}

       //Initialize with another DeviceConfigOption object
       HOST_DEVICE void init(DeviceConfigOptions<P> *const in){}
       
    };
      
    //Por ahora no he visto nada que haga falta copiar el GPU (04/08)    


    //Copy the configuration to the protected variables (CPU)
    void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco){}

       //Copy the local configuration to the DEVICE
    void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco){}

#endif





  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   *
   * Tau is the characteristic chemistry time
   */
  enum {TAU=0};

  /**
   * Constructor without arguments
   */
  ReacMFMHDTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ReacMFMHDTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Resize the physical data
   */
  virtual void resizePhysicalData(RealVector& physicalData);

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 1;
  }

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "ReacMFMHDTerm";
  }
  
}; // end of class ReacMFMHDTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_ReacMFMHDTerm_hh
