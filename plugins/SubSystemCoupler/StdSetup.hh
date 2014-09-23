#ifndef COOLFluiD_Numerics_SubSystemCoupler_StdSetup_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to setup the SubSystemCoupler Method
 *
 * @author Thomas Wuilbaut
 *
 */
class StdSetup : public CouplerCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdSetup();

  /**
   * Execute Processing actions
   */
  void execute();

protected: // data

}; // class StdSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StdSetup_hh

