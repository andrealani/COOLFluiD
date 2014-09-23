#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMFluxPrec_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMFluxPrec_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/AUSMFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an extension of the AUSM scheme with preconditioning 
 * suitable for low and high speeds
 *
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class AUSMFluxPrec : public AUSMFlux<UPDATEVAR> {
public: // classes
  
  /**
   * Constructor
   */
  AUSMFluxPrec(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~AUSMFluxPrec();
  
protected:
  
  /**
   * Compute the flux
   */
  virtual void computeMassFluxImpl(const CFuint eulerID,
				   const CFuint nbEulerEqs,
				   const std::vector<CFuint>& eulerVarIDs,
				   RealVector& result);
  
protected:
  
  /// interface pressure diffusion mass flux
  CFreal m_mfluxp;
  
}; // end of class AUSMFluxPrec

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMFluxPrec.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMFluxPrec_hh
