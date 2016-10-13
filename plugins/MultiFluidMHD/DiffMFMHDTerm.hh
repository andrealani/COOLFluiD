#ifndef COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDTerm_hh
#define COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDTerm_hh

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
 * This class represents the interface for a Diffusive model for the MultifluidMHD.
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez
 *
 */
class DiffMFMHDTerm : public Framework::BaseTerm {
  
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
      
    //For now, there is no need to copy the variables to de GPU (04/08)    


    //Copy the configuration to the protected variables (CPU)
    void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco){}

       //Copy the local configuration to the DEVICE
    void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco){}

#endif


  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor without arguments
   */
  DiffMFMHDTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~DiffMFMHDTerm();

  /**
   * Physical data size to be adapted 
   */
  CFuint getDataSize() const
  {
    //to test the Braginskii transport
    if (m_braginskiiTransport) {
      return 10; 				//2 Viscosities + 7 ThermConductiv (ion) + 1 ThermConductiv (neutral)
    }
    else{
      return 2*m_nbSpecies;
    }
  }
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();
   
  /**
   * Get the name
   */
  static std::string getName()
  {
    return "DiffMFMHDTerm";
  }
  
  RealVector& getDynViscosityDim()
  {
    return m_dynViscosityVec;
  }  
  
  RealVector& getThermConductivityDim()
  {
    return m_thermConductivityVec;
  }   
  
  /**
   * @return the number of species
   */
  virtual CFuint getNbSpecies() const {return m_nbSpecies;}   
  
  /**
   * @return is using the braginskii model
  */
  bool isBraginskii() const {return m_braginskiiTransport;}

  /**
   * @return is using the braginskii model
  */
  bool isExtendedDomain() const {return m_isExtended;}

  /**
   * @return TopHeight
  */
  CFreal getTopHeight() const {return m_topHeight;}

  /**
   * @return m_y0
  */
  CFreal getDampingHeight() const {return m_y0;}
  
  void computeNonInducedEMField(CFreal xCoord, CFreal yCoord);
  
   /**
   * Get the magnetic dipole field and dipole moment values
   */
  RealVector& getNonInducedEMField(CFreal x, CFreal y)
  {
    computeNonInducedEMField(x,y);
    return _NonInducedEMField;
  } 

  RealVector& getIncreasedDynViscosityDim()
  {
    return m_IncreasedDynViscosityVec;
  }
private:
  
  /// dimensional coefficient
  RealVector m_dynViscosityVec;
  
  /// dimensional coefficient
  RealVector m_thermConductivityVec;  
  
  /// number of species
  CFuint m_nbSpecies;
  
  /// Flag to use the Braginskii properties
  bool m_braginskiiTransport;
  
  /// dimensional coefficient to store the options input
  std::vector<CFreal> m_dynViscosity;
  
  /// dimensional coefficient to store the options input
  std::vector<CFreal> m_thermConductivity;  
  
  /// Non Induced part of Electromagnetic field
  RealVector _NonInducedEMField;

  ///introduced non induced electromagnetic field
  std::vector<CFreal> _nonInducedEMField;  

  /// Flag to use extended domain
  bool m_isExtended;

  /// Height of the upper boundary of the PHYSICAL domain
  CFreal m_topHeight;

  /// Increased viscosity in the damped domain
  std::vector<CFreal> m_IncreasedDynViscosity;

  /// dimensional coefficient of extended domain
  RealVector m_IncreasedDynViscosityVec;

  /// Height of the center of the tanh used to damp the viscosity smoothly
  CFreal m_y0;

}; // end of class DiffMFMHDTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_DiffMFMHDTerm_hh
