#ifndef COOLFluiD_Physics_MHD_MHDTerm_hh
#define COOLFluiD_Physics_MHD_MHDTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics { 

    namespace MHD { 
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MHDTerm.
 * 
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */
class MHDTerm : public Framework::BaseTerm {
public:
  
  /**
   * Enumerator defining the potential B type
   */
  enum PotentialBType {NONE=0, DIPOLE=1, PFSS=2};
  
#ifdef CF_HAVE_CUDA
  /// nested class defining local options
  template <typename P = NOTYPE>
  class DeviceConfigOptions {
  public:
    /// constructor
    HOST_DEVICE DeviceConfigOptions() {}
    
    /// destructor
    HOST_DEVICE virtual ~DeviceConfigOptions() {}
    
    /// initialize with another object of the same kind
    HOST_DEVICE void init(DeviceConfigOptions<P> *const in) 
    {
      potentialBType = in->potentialBType;
      gamma          = in->gamma;
      refSpeed       = in->refSpeed;
      mX             = in->mX;
      mY             = in->mY;
      mZ             = in->mZ;
    }
    
    /// potential B type
    PotentialBType potentialBType; 
    
    /// gamma
    CFreal gamma; 
    
    /// reference speed for the projection scheme
    CFreal refSpeed;
    
    /// x-component of the magnetic dipole moment (mX)
    CFreal mX;
    
    /// y-component of the magnetic dipole moment (mY)
    CFreal mY;
    
    /// z-component of the magnetic dipole moment (mZ)
    CFreal mZ;
  };
  
  /// copy the local configuration options to the Framework::DEVICE
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {
    CudaEnv::copyHost2Dev(&dco->potentialBType, &_pBtype, 1);
    CudaEnv::copyHost2Dev(&dco->gamma, &_gamma, 1);
    CudaEnv::copyHost2Dev(&dco->refSpeed, &_refSpeed, 1);
    CudaEnv::copyHost2Dev(&dco->mX, &_mX, 1);
    CudaEnv::copyHost2Dev(&dco->mY, &_mY, 1);
    CudaEnv::copyHost2Dev(&dco->mZ, &_mZ, 1);
  }  
  
  /// copy the local configuration options to the Framework::DEVICE
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    dco->potentialBType = _pBtype;
    dco->gamma = _gamma;
    dco->refSpeed = _refSpeed;
    dco->mX = _mX;
    dco->mY = _mY;
    dco->mZ = _mZ;
  }     
#endif
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the 
   * physical data
   */
  enum {RHO=0, P=1, A=2, B=3, BX=4, BY=5, 
	BZ=6, V=7, VX=8, VY=9, VZ=10, GAMMA=11, XP=12, YP=13, ZP=14};
  
  /**
   * Constructor without arguments
   */
  MHDTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~MHDTerm();
  
  /**
   * Set physical data
   */
  virtual void setupPhysicalData();
  
  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 15;
  }
  
  /**
   * Get the reference speed for projection scheme
   */
  CFreal getRefSpeed() const 
  {
    return _refSpeed;
  }

  /**
   * Get the mass of the external object 
   */
  CFreal getMass() const
  {
    return _mass;
  }

  /**
   * Get the number of l modes to be used in PFSS reconstruction
   */
  CFuint getNbLModes() const
  {
    return _nbLModes;
  }

  /**
   * Get the reference proton density to non-dimensionalize certain source terms
   */
  CFreal getNRef() const
  {
    return _nRef;
  }

  /**
   * Get the reference magnetic field to non-dimensionalize certain source terms
   */
  CFreal getBRef() const
  {
    return _BRef;
  }

  /**
   * Get the reference length (i.e. radius of the external object) to non-dimensionalize certain source terms
   */
  CFreal getLRef() const
  {
    return _lRef;
  }

  /**
   * Get the reference temperature to non-dimensionalize the equations for the solar wind problem
   */
  CFreal getTRef() const
  {
    return _TRef;
  }
    
  /**
   * Get the dissipation coefficient for projection scheme
   */
  CFreal getDissipationCoefficient() const 
  {
    return _dissipCoeff;
  }

  /**
   * Get the x-component of the magnetic dipole moment (mX)
   */
  CFreal getMX() const 
  {
    return _mX;
  }

  /**
   * Get the y-component of the magnetic dipole moment (mY)
   */
  CFreal getMY() const 
  {
    return _mY;
  }

  /**
   * Get the z-component of the magnetic dipole moment (mZ)
   */
  CFreal getMZ() const 
  {
    return _mZ;
  }

  /**
   * Get the potential magnetic field type to model the coronal
   * magnetic field initially: Dipole or PFSS 
   */
  std::string getPotentialBType() const
  {
    return _potentialBType;
  }

  /**
   * Get the name of correction for projection scheme
   */
  std::string getCorrectionType() const 
  {
    return _correctionType;
  }

  /**
   * Get the name of the output file for divB errors
   */
  std::string getNameOutputFile() const 
  {
    return _nameOutputFile;
  }
  
  /**
   * Get the frequency of saving the output file for divB errors
   */
  CFuint getOutputFileSaveRate() const 
  {
    return _saveRate;
  }

  /**
   * Get the name of the first input file containing the spherical harmonics coefficients for 
   * reconstructing the initial solar coronal magnetic field 
   */
  std::string getNameBeginPFSSCoeffFile() const
  {
    return _nameBeginPFSSCoeffFile;
  }

  /**
   * Get the name of the second input file containing the spherical harmonics coefficients for 
   * interpolating the B0 field during the unsteady simulation
   */
  std::string getNameEndPFSSCoeffFile() const
  {
    return _nameEndPFSSCoeffFile;
  }

  /**
   * Get the name of the method to increase the accuracy in global magnetosphere simulations 
   * especially in the inner magnetosphere: ISLND or Boris 
   */
  std::string getNameAccuracyIncreaserMethod() const
  {
    return _nameAccuracyIncreaserMethod;
  }
  
  /**
   * Configures this object by complementing the 
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );
    
  /**
   * Get \f$\gamma\f$ 
   */
  CFreal getGamma() const {return _gamma;}

  /**
   * Get the polytropic index 
   */
  CFreal getPolytropicIndex() const {return _n;}

  /**
   * Get the ISLND wave speed limit (non-dimensional) 
   */
  CFreal getISLNDLimit() const {return _ISLNDLimit;}  

  /**
   * Get radius of the source surface for the PFSS model 
   */
  CFreal getRSource() const {return _rSource;}

  /**
   * Set the interpolation flag when the PFSS magnetic field 
   * is to be interpolated during unsteady simulations
   */
  void setInterpolationFlag(bool interpolationFlag)
  {
    _interpolationFlag = interpolationFlag;
  }

  /**
   * Get the interpolation flag when the PFSS magnetic field 
   * is to be interpolated during unsteady simulations
   */
  bool getInterpolationFlag() const {return _interpolationFlag;}

  /**
   * Get the name
   */
  static std::string getName() 
  {
    return "MHDTerm";
  }
  
protected:
  
  /// potential B type
  PotentialBType    _pBtype;
  
  /// gamma
  CFreal            _gamma;  
    
  /// polytropic index
  CFreal            _n;
  
  /// ISLND wave speed limit (non-dimensional)
  CFreal _ISLNDLimit;

  /// radius of the source surface for the PFSS model
  CFreal _rSource;
 
  /// reference speed necessary for matching 
  //of the units for the projection scheme defined interactively
  CFreal            _refSpeed;

  /// mass of the external object to be specified if different than the Sun 
  CFreal _mass;

  /// Number of l modes to be utilized in PFSS reconstruction
  CFuint	    _nbLModes;

  /// reference length (i.e. radius of the external object) to be specified if different than the radius of the Sun 
  /// to non-dimensionalize certain source terms
  CFreal _lRef;

  /// reference magnetic field to non-dimensionalize certain source terms
  CFreal _BRef;

  /// reference proton density to non-dimensionalize certain source terms
  CFreal _nRef;

  /// reference temperature to non-dimensionalize the equations for the solar wind problem 
  CFreal _TRef;

  /// dissipation coefficient for damping divB errors 
  /// for the projection scheme
  CFreal            _dissipCoeff;

  /// interpolation flag for PFSS magnetic field
  bool _interpolationFlag;

  /// x-component of the magnetic dipole moment (mX)
  CFreal _mX;

  /// y-component of the magnetic dipole moment (mY)
  CFreal _mY;

  /// z-component of the magnetic dipole moment (mZ)
  CFreal _mZ;
    
  /// Storage for choosing when to save the divB error output file
  CFuint _saveRate;
 
  /// The potential magnetic field type to model the coronal magnetic field initially: Dipole or PFSS 
  std::string _potentialBType; 
 
  /// Name of correction for projection scheme
  std::string _correctionType;

  /// Name of the first input file containing PFSS spherical harmonics coefficients 
  std::string _nameBeginPFSSCoeffFile;

  /// Name of the second input file containing PFSS spherical harmonics coefficients 
  std::string _nameEndPFSSCoeffFile;

  /// Name of Accuracy Increaser Method in global magnetosphere simulations: ISLND or Boris
  std::string _nameAccuracyIncreaserMethod;

  /// Name of Output File where to write the coefficients
  std::string _nameOutputFile;
  
}; // end of class MHDTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD
 
  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDTerm_hh
