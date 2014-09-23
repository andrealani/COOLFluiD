#ifndef COOLFluiD_Numerics_FiniteVolume_NullDiffusiveFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_NullDiffusiveFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeDiffusiveFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class DiffusiveVarSet;
  }
  
  namespace Numerics {

    namespace FiniteVolume {
      
      class DerivativeComputer;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an object computing a diffusive flux
 * with FV method
 *
 * @author Andrea Lani
 *
 */
class NullDiffusiveFlux : public ComputeDiffusiveFlux {
public:
   
  /**
   * Constructor
   */
  NullDiffusiveFlux(const std::string& name);

  /**
   * Default destructor
   */
  ~NullDiffusiveFlux();
  
  /**
   * Tells that this object is a Null object.
   */
  bool isNull() const
  {
    return true;
  }
  
  /**
   * Set the diffusive variable set
   */
  void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> 
			  diffVar);
  
  /**
   * Compute the flux in the current face
   */
  void computeFlux(RealVector& result);
  
  /**
   * Set the FluxSplitterData
   */
  void setFluxData(Common::SafePtr<Framework::FluxSplitterData> fluxData);
  
}; // end of class NullDiffusiveFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NullDiffusiveFlux_hh
