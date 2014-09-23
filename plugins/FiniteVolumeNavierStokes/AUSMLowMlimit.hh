#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMLowMlimit_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMLowMlimit_hh

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
 * This class computes the AUSM low Mach limit flux
 * @author Lani
 * @author VDH
 *
 */
template <class UPDATEVAR>
class AUSMLowMlimit : public AUSMFlux<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMLowMlimit(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMLowMlimit();
  
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
 
  /**													//?
   * Set up private data to prepare the simulation		//?
   */													//?
  virtual void setup();									//?
  
protected:

  virtual void computeMassFlux();		//interface mass flux
  virtual void computeIncompCorrectionTerm();   //incompressible correction term mp
  virtual void computePressureFlux();		//interface pressure flux
  
private:
/// acquaintance of the concrete variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffusiveVarSet;

  ///Dummy vector for the gradients
  std::vector<RealVector*> m_DummyGradients;
    
  CFreal m_Vinf;					// V infinite
  CFreal m_Lref;					// Reference length, shortest cell size
  CFreal m_nu;						// dynamic viscosity
  CFuint m_ChoiceLref;					// Choice of method to compute Lref
  CFuint m_ChoiceMp;					// Choice of method for computing mp
  CFreal m_uco;						// cut-off speed
  CFreal m_umax;					// higher bound speed
  CFreal m_ChoiceVisc;					// Choice of calculation of dynamic viscosity
  
  
}; // end of class AUSMLowMlimit

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMLowMlimit.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMLowMlimit_hh
