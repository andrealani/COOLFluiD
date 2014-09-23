#ifndef COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallIsothermalVTImpl_hh
#define COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallIsothermalVTImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/StrongImplBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong no slip wall BC with fixed wall
 * temperature for Navier Stokes 2D with [ ... V T] variables
 *
 * @author Andrea Lani
 */
class StrongNoSlipWallIsothermalVTImpl : public FluctSplit::StrongImplBC {

public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  StrongNoSlipWallIsothermalVTImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongNoSlipWallIsothermalVTImpl();

  /**
   * Set up private data
   */
  void setup();
  
protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Compute in  the 2D case
   */
  void compute2D();
  
  /**
   * Compute in  the 3D case
   */
  void compute3D();
  
private:
  
  /// dimensional wall temperature
  CFreal _TWall;
        
}; // end of class StrongNoSlipWallIsothermalVTImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallIsothermalVTImpl_hh
