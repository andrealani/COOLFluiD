#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletRiga_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletRiga_hh

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
class SubOutletRiga : public BASE {
public:

  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubOutletRiga(const std::string& name);

  /**
   * Default destructor
   */
  ~SubOutletRiga();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);


private:

  /// user defined inlet scalar variables
  std::vector<CFreal> m_ye;
}; // end of class SubOutletEulerPvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubOutletRiga.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEulerPvt_hh
