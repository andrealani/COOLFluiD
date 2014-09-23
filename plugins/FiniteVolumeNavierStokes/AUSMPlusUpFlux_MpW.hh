#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_MpW_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_MpW_hh

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
class AUSMPlusUpFlux_MpW : public AUSMFlux<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMPlusUpFlux_MpW(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPlusUpFlux_MpW();
  
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

protected:

  
  virtual void computeMassFlux();
  virtual void computePressureFlux();
  virtual void computeIncompCorrectionTerm();
  virtual CFreal correctMachInf(CFreal oldMach) const
  {
    return oldMach;
  }
  
private:
     
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
  
}; // end of class AUSMPlusUpFlux_MpW

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMPlusUpFlux_MpW.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_MpW_hh
