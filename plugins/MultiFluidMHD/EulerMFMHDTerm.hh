#ifndef COOLFluiD_Physics_MultiFluidMHD_EulerMFMHDTerm_hh
#define COOLFluiD_Physics_MultiFluidMHD_EulerMFMHDTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Maxwell/MaxwellProjectionTerm.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a multi-fluid (GPU-enabled)
 * Euler convective
 * physical term
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 * @author Isaac Alonso
 *
 */
class EulerMFMHDTerm : public Maxwell::MaxwellProjectionTerm {
  
  enum {START=Maxwell::MaxwellProjectionTerm::END};
  
public:



#ifdef CF_HAVE_CUDA
    
   //Nested class defining local options
   template <typename P = NOTYPE >   //Need to ask about this
   class DeviceConfigOptions{
   public:
       //Constructor
       HOST_DEVICE DeviceConfigOptions() {}
       //Destructor
       HOST_DEVICE virtual ~DeviceConfigOptions() {}

       //Initialize with another DeviceConfigOption object
       HOST_DEVICE void init(DeviceConfigOptions<P> *const in){
          CFLog(VERBOSE, "HOST_DEVICE EulerMFMHDTerm::DeviceConfigOptions::init(); \n");
          gamma = in->gamma;
	  omega = in->omega;
	  K = in->K;
	  molecularMass1 = in->molecularMass1;
	  molecularMass2 = in->molecularMass2;
	  molecularMass3 = in->molecularMass3;
	  epsilon = in->epsilon;
	  mu = in->mu;
	  lightSpeed = in->lightSpeed;

          //NonInducedEMField = in->NonInducedEMField;
          isLeake = in->isLeake;
          dim = in->dim;
          is2DHalf = in->is2DHalf;
          firstSpecies = in->firstSpecies;
          firstVelocity = in->firsrVelocity;
          firstTemperature = in->firstTemperature;
          nbSpecies = in->nbSpecies;
          nbMomentum = in->nbMomentum;
          nbEnergyEqs = in->nbEnergyEqs;
          divECleaningConst = in->divECleaningConst;
          divBCleaningConst = in->divBCleaningConst;
       }
       
       CFreal gamma;
       CFreal omega;
       CFreal K;
       CFreal molecularMass1;  
       CFreal molecularMass2;  
       CFreal molecularMass3; 
       CFreal epsilon;
       CFreal mu;
       CFreal lightSpeed;  
       CFreal	divECleaningConst; 
       CFreal	divBCleaningConst; 
       CFreal NonInducedEMField[6];  

       bool isLeake;
       
       CFuint dim;
       bool is2DHalf;
       
       CFuint firstSpecies;
       CFuint firstVelocity;
       CFuint firstTemperature;

       CFuint nbSpecies;
       CFuint nbMomentum;
       CFuint nbEnergyEqs;
    };
      

    //Copy the configuration to the protected variables (CPU)
    void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco)
    {
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions::copyConfigOptions()  START \n \n");



          _dim = Framework::PhysicalModelStack::getActive()->getDim();
           
          _is2DHalf = Framework::PhysicalModelStack::getActive()->getImplementor()->is2DHalf();

          _firstSpecies = getFirstScalarVar(0);
          _firstVelocity = getFirstScalarVar(1);
          _firstTemperature = getFirstScalarVar(2);

          _nbSpecies = getNbScalarVars(0);
          _nbMomentum = getNbScalarVars(1);
          _nbEnergyEqs = getNbScalarVars(2);



	  dco->gamma = _gamma;
	  dco->omega = _omega;
	  dco->K = _K;	 
	  dco->molecularMass1 = _molecularMass1;	 
	  dco->molecularMass2 = _molecularMass2; 
	  dco->molecularMass3 = _molecularMass3; 
	  dco->epsilon = _epsilon;
	  dco->mu = _mu;    
          dco->lightSpeed = _lightSpeed;    
          dco->divECleaningConst = _divECleaningConst;
          dco->divBCleaningConst = _divBCleaningConst;  
          dco->NonInducedEMField[0] = _NonInducedEMField[0];
          dco->NonInducedEMField[1] = _NonInducedEMField[1];
          dco->NonInducedEMField[2] = _NonInducedEMField[2];
          dco->NonInducedEMField[3] = _NonInducedEMField[3];
          dco->NonInducedEMField[4] = _NonInducedEMField[4];
          dco->NonInducedEMField[5] = _NonInducedEMField[5];
          dco->isLeake = _isLeake;  
          dco->dim = _dim;       
          dco->is2DHalf = _is2DHalf; 
          dco->firstSpecies = _firstSpecies;
          dco->firstVelocity = _firstVelocity;
          dco->firstTemperature = _firstTemperature;
          dco->nbSpecies = _nbSpecies;
          dco->nbMomentum = _nbMomentum;
          dco->nbEnergyEqs = _nbEnergyEqs;
    }

       //Copy the local configuration to the DEVICE
    void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco)
    { 
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions::copyConfigOptionsToDevice()  START \n \n");

          
          _dim = Framework::PhysicalModelStack::getActive()->getDim();
           
          _is2DHalf = Framework::PhysicalModelStack::getActive()->getImplementor()->is2DHalf();

          _firstSpecies = getFirstScalarVar(0);
          _firstVelocity = getFirstScalarVar(1);
          _firstTemperature = getFirstScalarVar(2);

          _nbSpecies = getNbScalarVars(0);
          _nbMomentum = getNbScalarVars(1);
          _nbEnergyEqs = getNbScalarVars(2);

          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions gamma = " << _gamma  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions omega = " << _omega  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions K = " << _K  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions molecularMass1 = " << _molecularMass1  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions molecularMass2 = " << _molecularMass2  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions molecularMass3 = " << _molecularMass3  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions epsilon = " << _epsilon  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions lighspeed = " << _lightSpeed  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions divBCleaningConst = " << _divBCleaningConst  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions divECleaningConst = " << _divECleaningConst  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions isLeake = " << _isLeake  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions dim = " << _dim  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions is2DHalf = " << _is2DHalf  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions firstSpecies = " << _firstSpecies  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions firstVelocity = " << _firstVelocity  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions firstTemperature = " << _firstTemperature  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions nbSpecies = " << _nbSpecies  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions nbEnergyEqs = " << _nbEnergyEqs  << "\n");
          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions NonInducedEM = " << _NonInducedEMField[0]  
                       << "\t" << _NonInducedEMField[1] << "\t" << _NonInducedEMField[2] << "\t" << _NonInducedEMField[3] 
                       << "\t" << _NonInducedEMField[4] << "\t" << _NonInducedEMField[5] << "\n");


          CudaEnv::copyHost2Dev(&dco->mu, &_mu, 1);
          CudaEnv::copyHost2Dev(&dco->gamma, &_gamma, 1);
          CudaEnv::copyHost2Dev(&dco->omega, &_omega, 1);
          CudaEnv::copyHost2Dev(&dco->K, &_K, 1);
          CudaEnv::copyHost2Dev(&dco->molecularMass1, &_molecularMass1, 1);
          CudaEnv::copyHost2Dev(&dco->molecularMass2, &_molecularMass2, 1);
          CudaEnv::copyHost2Dev(&dco->molecularMass3, &_molecularMass3, 1);
          CudaEnv::copyHost2Dev(&dco->epsilon, &_epsilon, 1);
          CudaEnv::copyHost2Dev(&dco->mu, &_mu, 1);
          CudaEnv::copyHost2Dev(&dco->lightSpeed, &_lightSpeed, 1);
          CudaEnv::copyHost2Dev(&dco->divBCleaningConst, &_divBCleaningConst, 1);
          CudaEnv::copyHost2Dev(&dco->divECleaningConst, &_divECleaningConst, 1);
          CudaEnv::copyHost2Dev(&dco->NonInducedEMField[0], &_NonInducedEMField[0], 6);
          CudaEnv::copyHost2Dev(&dco->isLeake, &_isLeake, 1);
          CudaEnv::copyHost2Dev(&dco->dim, &_dim, 1);
          CudaEnv::copyHost2Dev(&dco->is2DHalf, &_is2DHalf, 1);
          CudaEnv::copyHost2Dev(&dco->firstSpecies, &_firstSpecies, 1);
          CudaEnv::copyHost2Dev(&dco->firstVelocity, &_firstVelocity, 1);
          CudaEnv::copyHost2Dev(&dco->firstTemperature, &_firstTemperature, 1);
          CudaEnv::copyHost2Dev(&dco->nbSpecies, &_nbSpecies, 1);
          CudaEnv::copyHost2Dev(&dco->nbMomentum, &_nbMomentum, 1);
          CudaEnv::copyHost2Dev(&dco->nbEnergyEqs, &_nbEnergyEqs, 1);
	  CudaEnv::copyHost2Dev(&dco->NonInducedEMField[0], &_NonInducedEMField[0], 6);

          CFLog(VERBOSE, "EulerMFMHDTerm::DeviceConfigOptions::copyConfigOptionsToDevice()  END \n \n");

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
  enum {RHO=START,XP=START+1, YP=START+2, ZP=START+3};
  
  /**
   * Constructor without arguments
   */
  EulerMFMHDTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~EulerMFMHDTerm();

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
    return Maxwell::MaxwellProjectionTerm::getDataSize() + 4;   //Maxwell::MaxwellProjectionTerm::getDataSize() = 6
  }

  /**
   * Get the start of the scalar vars data (mass fractions)
   */
  virtual CFuint getFirstScalarVar(CFuint i) const
  {
    return getDataSize();
  }

  /**
   * Get the number of scalar vars
   */
  virtual CFuint getNbScalarVars(CFuint i) const
  {
    return 0;
  }

  /**
   * Get \f$\gamma\f$
   */
  CFreal getGamma() const
  {
    return _gamma;
  }
   
  /**
   * Get rotational velocity in radiants/sec
   */
  CFreal getOmega() const
  {
    return _omega;
  }
  
  /**
   * Get gas constant in [J/mol K]
   */ 
  CFreal getK() const
  {
    return _K;
  } 
  
  /**
   * Get electronMass in [kg/mol]
   */ 
  CFreal getMolecularMass1() const
  {
    return _molecularMass1;
  } 
  
   /**
   * Get molecularMass2 in [kg/mol]
   */ 
  CFreal getMolecularMass2() const
  {
    return _molecularMass2;
  } 
  
   /**
   * Get molecularMass3 in [kg/mol]
   */ 
  CFreal getMolecularMass3() const
  {
    return _molecularMass3;
  }   
  
  CFreal getPermittivity() const
  {
    return _epsilon;
  } 
  
  CFreal getPermeability() const
  {
    return _mu;
  }
  
  CFreal getLightSpeed() const
  {
    return _lightSpeed;
  }

  bool isLeake() const
  {
    return _isLeake;
  }
  
  /**
   * Get index of pressure
   */
  //static CFuint getPressureTerm()//should be removed?
  //{
  //	return EulerMFMHDTerm::P;
  //}
  
  /**
   * Get dimensional pressure from given pState which is stored in the state
   */
  CFreal getPressureDimFromState(CFreal pState)
  {
    return pState;
  }
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "EulerMFMHDTerm";
  }
  
   /**
   * Compute the Non Induced Electromagnetic field
   */
  void computeNonInducedEMField(CFreal xCoord, CFreal yCoord);

  /**
   * Get the magnetic dipole field and dipole moment values
   */
  RealVector& getNonInducedEMField(CFreal x, CFreal y, CFreal z)   //El ultimo argumento no estaba
  {
    computeNonInducedEMField(x,y);
    return _NonInducedEMField;
  }

protected:

  /// specific heat ratio
  CFreal _gamma;
  
  /// rotational velocity
  CFreal _omega;
   
  /// gas constant [J/ K]
  CFreal _K;

  /// electron mass [kg/mol]
  CFreal _molecularMass1;
  
  /// proton mass [kg/mol]  
  CFreal _molecularMass2;
  
  /// neutral mass(Atomic Hydrogen) [kg/mol]    
  CFreal _molecularMass3;
  
  ///permittivity of free space [F/m]
  CFreal _epsilon;
  
  /// permeabilitty of free space[H/m]
  CFreal _mu;
  
  ///Speed of light
  CFreal _lightSpeed;
  
  /// Non Induced part of Electromagnetic field
  RealVector _NonInducedEMField;
  
  ///introduced non induced electromagnetic field
  std::vector<CFreal> _nonInducedEMField;   

  /// Flag to use Leake's model in the convective term
  bool _isLeake;

  bool _is2DHalf;

  /// Dimension
  CFuint _dim;
       
  CFuint _firstSpecies;
  CFuint _firstVelocity;
  CFuint _firstTemperature;

  CFuint _nbSpecies;
  CFuint _nbMomentum;
  CFuint _nbEnergyEqs;



}; // end of class EulerMFMHDTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_EulerMFMHDTerm_hh
