#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/AUSMFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM flux
 *
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class AUSMPlusUpFlux : public AUSMFlux<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMPlusUpFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPlusUpFlux();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the data from the supplied arguments.
   * @param args configuration arguments
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    AUSMFlux<UPDATEVAR>::configure(args);
  }
  
  // /*Set up private data*/

  //virtual void setup();

protected:

  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFlux();
  
  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux();
  
  /**
   * Compute the flux respectively to variable Var
   */
  virtual void ComputedFluxdVar(CFuint iVar, RealMatrix &JacobMatrix);


  /// Correct the given Mach number
  virtual CFreal correctMachInf(CFreal oldMach) const
  {
    return oldMach;
  }

  /// Compute the left flux jacobian
  virtual void computeLeftJacobian();

  /// Compute the right flux jacobian
  virtual void computeRightJacobian();

  /// Compute the jacobian with respect to Puvt or PVt variables
  virtual void computeJacobianDPvt(CFuint side, RealMatrix &JacobMatrix);

  /// Compute the jacobian with respect to cons variables
  virtual void computeJacobianDCons(CFuint side, RealMatrix &JacobMatrix);

  /// compute flux derivative with respect to pressure
  void dFdP(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);
  
  /// compute flux derivative with respect to density
  void dFdRho(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);
  
  /// compute flux derivative with respect to temperature
  void dFdT(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);
  
  /// compute flux derivative with respect to velocity u
  void dFdU(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);
  
  /// compute flux derivative with respect to velocity v
  void dFdV(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);
  
  /// compute flux derivative with respect to velocity w
  void dFdW(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);
    
  /// compute flux derivative with respect to k
  void dFdK(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);

  /// compute flux derivative with respect to rhou
  void dFdRhoU(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);

      /// compute flux derivative with respect to rhov
  void dFdRhoV(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);

      /// compute flux derivative with respect to rhow
  void dFdRhoW(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);

      /// compute flux derivative with respect to rhoE
  void dFdRhoE(CFuint side, CFuint iVar, RealMatrix &JacobMatrix);
 
private:
  
  /// preconditioning coefficient 
  CFreal m_fa;
  
  /// P5 plus coefficient 
  CFreal m_P5Plus;
  
  /// P5 minus coefficient 
  CFreal m_P5Minus;
    
  /// user defined coefficient for Ku
  CFreal m_coeffKu;
  
  /// user defined coefficient for Kp
  CFreal m_coeffKp;
  
  /// user defined coefficient for sigma
  CFreal m_coeffSigma;

  /// mach infinity
  CFreal m_machInf;
  
  /// beta  coefficient
  CFreal m_beta;

  /// derivative of the speed (u,v,w) with respect to a certain variable
  // RealVector m_dVLdVar;
  // RealVector m_dVRdVar;
  CFreal m_duLdVar;
  CFreal m_duRdVar;
  CFreal m_dvLdVar;
  CFreal m_dvRdVar;
  CFreal m_dwLdVar;
  CFreal m_dwRdVar;

  /// derivative of some variable to compute the Jacobian
  CFreal m_dqnLdVar;
  CFreal m_dqnRdVar;

  CFreal m_dhLdVar;
  CFreal m_dhRdVar;

  CFreal m_drhoLdVar;
  CFreal m_drhoRdVar;

  CFreal m_dpLdVar;
  CFreal m_dpRdVar;


}; // end of class AUSMPlusUpFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMPlusUpFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_hh
