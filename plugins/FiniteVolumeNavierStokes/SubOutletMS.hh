#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletMS_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletMS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic outlet command
 *
 * @author Andrea Lani
 *
 */
template <class BASE>
class SubOutletMS : public BASE {
public:

  /**
   * Constructor
   */
  SubOutletMS(const std::string& name);

  /**
   * Default destructor
   */
  ~SubOutletMS();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

}; // end of class SubOutletEulerPvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubOutletMS.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEulerPvt_hh
