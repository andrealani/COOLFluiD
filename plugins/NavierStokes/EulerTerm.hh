#ifndef COOLFluiD_Physics_NavierStokes_EulerTerm_hh
#define COOLFluiD_Physics_NavierStokes_EulerTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an Euler convective
 * physical term
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 */
class EulerTerm : public Framework::BaseTerm {
public:
    
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
      gamma = in->gamma;
      R = in->R;
    }
    
    /// gamma
    CFreal gamma; 
  
    /// R
    CFreal R;
  };
  
  /// copy the local configuration options to the Framework::DEVICE
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {
    CudaEnv::copyHost2Dev(&dco->gamma, &_gamma, 1);
    CudaEnv::copyHost2Dev(&dco->R, &_RDim, 1);
  }  
  
  /// copy the local configuration options to the Framework::DEVICE
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    dco->gamma = _gamma;
    dco->R = _RDim;
  }     
  #endif

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Enumerator defining the mapping between the variable name 
   * and its position in the physical data
   * @pre EulerTerm::P is the computed pressure, 
   *      either static pressure (compressible case) or
   *      pressure fluctuation (incompressible case) 
   */
  enum {RHO=0, P=1, H=2, E=3, A=4, T=5, V=6, VX=7, VY=8, VZ=9, 
	GAMMA=10, XP=11, YP=12, ZP=13};
  
  /**
   * Constructor without arguments
   */
  EulerTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~EulerTerm();

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
  virtual CFuint getDataSize() const {return 14;}
  
  /**
   * Get the start of the scalar vars data (mass fractions)
   */
  virtual CFuint getFirstScalarVar(CFuint i) const {return 14;}
  
  /**
   * Get the number of scalar vars
   */
  virtual CFuint getNbScalarVars(CFuint i) const {return 0;}
  
  /// Get the number of sets of scalar vars
  virtual CFuint getNbScalarVarSets() const {return 0;}
  
  /**
   * Get the constant of real gas
   */
  CFreal getR() const {return _RRef;}

  /**
   * Get the constant of real gas
   */
  CFreal getRdim() const {return _RDim;}

  /**
    * Get \f$\gamma\f$
    */
  CFreal getGamma() const {return _gamma;}

  /**
   * Get the Vel infinity
   */
  CFreal getVelInf() const {return _uInf;}
  
  /**
   * Get rotational velocity in radiants/sec
   */
  CFreal getOmega() const {return _omega;}
  
  /**
   * Get the Mach infinity
   */
  CFreal getMachInf() const {return _machInf;}

  /**
   * Get the specific heat at constant pressure (Cp) of gas
   */
  CFreal getCp() const {return getR()/(1.-(1./_gamma));}
  
  /**
   * Get the specific heat at constant volume (Cv) of gas
   */
  CFreal getCv() const {return getR()/(_gamma - 1.);}
  
  /**
   * Get the Mach infinity
   */
  CFreal getVelRef() const {return m_refPhysicalData[EulerTerm::V];}
    
  /**
   * Get the reference temperature
   */
  CFreal getTempRef() const {return _tempRef;}

  /**
   * Get the reference pressure
   */
  CFreal getPressRef() const {return _pRef;}
  
  /**
   * Get the pressure infinity
   */
  CFreal getPressInf() const {return _p0Inf;}
  
  /**
   * Get the pressure infinity (compressible case)
   */
  CFreal getPressInfComp() const {return _p0InfComp;}
  
  /**
   * Get dimensional pressure from given pState which is stored in the state
   */
  CFreal getPressureFromState(CFreal pState) 
  {
    return pState + _p0Inf/m_refPhysicalData[EulerTerm::P];
  }
  
  /**
   * Get the density
   * @pre the input pressure if the total static pressure, as from getPressureFromState()
   */
  CFreal getDensity(CFreal staticP, CFreal staticT) const
  {
    // pure incompressible vs. low Mach number case
    return (_isPureIncomp) ? _rho/m_refPhysicalData[EulerTerm::RHO] : staticP/(_RRef*staticT);
  }
  
  /**
   * Set the reference R constant
   */
  void setRRef(CFreal R) {_RRef = R;}

  /**
   * Set the reference pressure
   */
  void setPressRef(CFreal pRef) {_pRef = pRef;}

  /**
   * Set the reference temperature
   */
  void setTempRef(CFreal tempRef) {_tempRef = tempRef;}
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Get the name
   */
  static std::string getName() {return "EulerTerm";}
  
  /// Tell if the flow model is incompresible
  bool isIncompressible() const {return (_p0Inf > 0.);}
  
protected:
 
  /// reference R
  CFreal _RRef;
  
  /// dimensional R
  CFreal _RDim;
  
  /// specific heat ratio
  CFreal _gamma;

  /// reference static temperature
  CFreal _tempRef;

  /// reference static pressure
  CFreal _pRef;

  /// Mach infinity
  CFreal _machInf;
  
  /// u infinity
  CFreal _uInf;
    
  /// thermodynamic pressure infinity
  CFreal _p0Inf;
  
  /// thermodynamic pressure infinity (compressible case)
  CFreal _p0InfComp;
  
  /// density
  CFreal _rho;
  
  /// rotational velocity
  CFreal _omega;
  
  /// flag is true if we are using pure incompressible case,
  /// else we are using low Mach number case
  bool _isPureIncomp;

}; // end of class EulerTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_EulerTerm_hh
