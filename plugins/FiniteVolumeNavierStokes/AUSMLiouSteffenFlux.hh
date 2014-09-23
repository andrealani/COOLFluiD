#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMLiouSteffenFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMLiouSteffenFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/AUSMFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Lax-Friedrichs flux
 *
 * @author Andrea Lani
 * @author Sarp Yalim
 *
 */
template <class UPDATEVAR>
class AUSMLiouSteffenFlux : public AUSMFlux<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMLiouSteffenFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMLiouSteffenFlux();
  
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

  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFlux();
  
  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux();
  
private:
  
  /// user defined coefficient for the calculation of the pressure
  CFuint _choiceKp;
  
}; // end of class AUSMLiouSteffenFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMLiouSteffenFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMLiouSteffenFlux_hh
