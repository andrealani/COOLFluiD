#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpIcpFlux_cp_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpIcpFlux_cp_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/AUSMFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    } 
  }
 
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
class AUSMPlusUpIcpFlux_cp : public AUSMFlux<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMPlusUpIcpFlux_cp(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPlusUpIcpFlux_cp();
  
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
 
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
protected:

  
  virtual void computeMassFlux();
  virtual void computePressureFlux();
  virtual void computeIncompCorrectionTerm();
  virtual CFreal correctMachInf(CFreal oldMach) const
  {
    return oldMach;
  }
  
private:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffusiveVarSet;

  ///Dummy vector for the gradients
  std::vector<RealVector*> m_DummyGradients;

  CFreal m_fa;						/// preconditioning coefficient
  CFreal m_coeffKu;					/// user defined coefficient for Ku
  CFreal m_coeffKp;					/// user defined coefficient for Kp
  CFreal m_coeffSigma;					/// user defined coefficient for sigma
  CFreal m_machInf; 					/// mach infinity
  CFreal m_beta;					/// beta  coefficient
  CFreal m_Vinf;					// V infinite
  CFreal m_Lref;					// Reference length, shortest cell size
  CFreal m_nu;						// dynamic viscosity
  CFuint m_ChoiceLref;					// Choice of method to compute Lref
  CFuint m_ChoiceMp;					// Choice of method for computing mp
  CFreal m_uco;						// cut-off speed
  CFreal m_umax;					// higher bound speed
  CFreal m_Ml;						// Mach limit of weighting function for mp
  CFreal m_ChoiceVisc;					// Choice of calculation of dynamic viscosity
  
}; // end of class AUSMPlusUpIcpFlux_cp

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMPlusUpIcpFlux_cp.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpIcpFlux_cp_hh
