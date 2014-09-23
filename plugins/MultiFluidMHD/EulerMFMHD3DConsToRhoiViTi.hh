#ifndef COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD3DConsToRhoiViTi_hh
#define COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD3DConsToRhoiViTi_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarVarSetBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

      class EulerMFMHDTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative variable to primitive 
 * [rhoi vi ui Ti]
 * AAL: This one is NOT used when you select RhoiViTi as vars
 *
 * @author Alejandro Alvarez
 *
 */
class EulerMFMHD3DConsToRhoiViTi : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<EulerMFMHDTerm> PTERM;


  /**
   * Default constructor without arguments
   */
  EulerMFMHD3DConsToRhoiViTi(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~EulerMFMHD3DConsToRhoiViTi();
  
  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
private:

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerMFMHDTerm> _model;

  /// Array with the particle mass of the species
  RealVector			   _m_i;  
}; // end of class EulerMFMHD3DConsToRhoiViTi

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD3DConsToRhoiViTi_hh
