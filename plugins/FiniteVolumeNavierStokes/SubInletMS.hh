#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletMS_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletMS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command with the initial conditions given for tTotal, pTotal and alpha
   *
   * @author Andrea Lani
   * @author Radek Honzatko
   *
   */
template <class BASE>
class SubInletMS : public BASE {
public:
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletMS(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletMS();
  
  /**
   * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );
  
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
  
}; // end of class SubInletMS

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubInletMS.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletMS_hh
