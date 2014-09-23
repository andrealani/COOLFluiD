#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtPis_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtPis_hh

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
 * @author Radek Honzatko
 *
 */
class SubOutletEuler2DPuvtPis : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubOutletEuler2DPuvtPis(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~SubOutletEuler2DPuvtPis();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

private: // data

  /// isoentropic pressure
  CFreal P_is;

}; // end of class SubOutletEuler2DPuvtPis

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtPis_hh
