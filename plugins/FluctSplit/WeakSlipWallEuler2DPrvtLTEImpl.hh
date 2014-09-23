#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DPrvtLTEImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DPrvtLTEImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakSlipWallEuler2DImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }
    


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for Euler2D
 *
 * @author Andrea Lani
 *
 */
class WeakSlipWallEuler2DPrvtLTEImpl : public WeakSlipWallEuler2DImpl {
public:

  /**
   * Constructor.
   */
  WeakSlipWallEuler2DPrvtLTEImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallEuler2DPrvtLTEImpl();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
  
protected:
  
  /**
   * Compute the normal flux and the corresponding jacobian 
   * dF/dP with P = [p rhoU rhoV T]
   */
  void computeNormalFluxAndJacob(const Framework::State& state,
				 const RealVector& normal,
				 RealVector& flux,
				 RealMatrix& fluxJacob);

private:
  
  /// pointer to the library computing the thermodinamics
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// tolerance for numerical derivation
  CFreal _tol;

  /// temporary array for density enthalpy energy
  RealVector _dhe;
  
}; // end of class WeakSlipWallEuler2DPrvtLTEImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2DPrvtLTEImpl_hh
