#ifndef COOLFluiD_Numerics_FluctSplitNEQ_StrongNoSlipWallIsothermalTCNEQ2DImpl_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_StrongNoSlipWallIsothermalTCNEQ2DImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/StrongImplBC.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  


    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong no slip wall BC with fixed wall
 * temperature for Navier Stokes 2D with thermo-chemical NEQ
 *
 * @author Andrea Lani
 */
class StrongNoSlipWallIsothermalTCNEQ2DImpl : public FluctSplit::StrongImplBC {

public:
  
  typedef Physics::NavierStokes::MultiScalarVarSet
  <Physics::NavierStokes::Euler2DVarSet> TCNEQ2DVarSet;
      
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  StrongNoSlipWallIsothermalTCNEQ2DImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongNoSlipWallIsothermalTCNEQ2DImpl();
  
  /**
   * Set up private data
   */
  void setup();
  
protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

private:
  
  /// dimensional wall temperature
  CFreal _TWall;
        
}; // end of class StrongNoSlipWallIsothermalTCNEQ2DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_StrongNoSlipWallIsothermalTCNEQ2DImpl_hh
