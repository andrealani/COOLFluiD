#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMPWFluxMHD_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMPWFluxMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeMHD/AUSMFluxMHD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM flux
 *
 * @author Alejandro Alvarez Laguna
 *
 */
template <class UPDATEVAR>
class AUSMPWFluxMHD : public AUSMFluxMHD<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMPWFluxMHD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPWFluxMHD();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the data from the supplied arguments.
   * @param args configuration arguments
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data
   */
  virtual void setup();

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
  
  CFreal m_beta;					//beta coeficient
  CFreal m_alpha;					//alpha coeficient
  CFreal m_Vinf;					// V infinite
  CFreal m_Lref;					// Reference length, shortest cell size
  CFreal m_nu;						// dynamic viscosity
  CFuint m_ChoiceLref;					// Choice of method to compute Lref
  CFuint m_ChoiceMp;					// Choice of method for computing mp
  CFreal m_uco;						// cut-off speed
  CFreal m_umax;					// higher bound speed
  CFreal m_Ml;						//limit M for weighting of Mp

  
}; // end of class AUSMPWFluxMHD

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMPWFluxMHD.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPWFluxMHD_hh
