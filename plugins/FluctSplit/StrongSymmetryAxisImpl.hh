#ifndef COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallIsothermalVTImpl_hh
#define COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallIsothermalVTImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/StrongImplBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong symmetry Axis BC 
 *  with [ ... V T] variables.
 * It is adapted from StrongNoSlipWallVTImpl (at Revision: 14439)
 * @author Andrea Lani
 * @author Jesus Garicano Mena
 */
class StrongSymmetryAxisImpl : public FluctSplit::StrongImplBC {

public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  StrongSymmetryAxisImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongSymmetryAxisImpl();

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

}; // end of class StrongSymmetryAxisImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSymmetryAxisImpl_hh
