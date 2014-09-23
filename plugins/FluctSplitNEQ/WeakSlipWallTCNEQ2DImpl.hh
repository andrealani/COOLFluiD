#ifndef COOLFluiD_Numerics_FluctSplitNEQ_WeakSlipWallTCNEQ2DImpl_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_WeakSlipWallTCNEQ2DImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakSlipWall2DImpl.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
    

    
    namespace FluctSplitNEQ {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an implicit weak slip wall bc (E. V. Weide) for Euler2D
 * with thermo-chemical NEQ
 *
 * @author Andrea Lani
 *
 */
class WeakSlipWallTCNEQ2DImpl : public FluctSplit::WeakSlipWall2DImpl {
public:
  
  typedef Physics::NavierStokes::MultiScalarVarSet
  <Physics::NavierStokes::Euler2DVarSet> TCNEQ2DVarSet;
  
  /**
   * Constructor.
   */
  WeakSlipWallTCNEQ2DImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallTCNEQ2DImpl();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected:
  
  /**
   * Compute the normal flux and the corresponding jacobian
   */
  void computeNormalFluxAndJacob(const Framework::State& state,
				 const RealVector& normal,
				 RealVector& flux,
				 RealMatrix& fluxJacob);
protected:
  
  /// Euler TCNEQ var set
  Common::SafePtr<TCNEQ2DVarSet> _varSet;
  
  /// physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// array to store the mass fractions
  RealVector _ys;
  
  /// species molar masses
  RealVector _Ms;
  
  /// vibrational temperatures
  RealVector _tvDim;

}; // end of class WeakSlipWallTCNEQ2DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_WeakSlipWallTCNEQ2DImpl_hh
