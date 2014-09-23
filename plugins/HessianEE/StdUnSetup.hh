#ifndef COOLFluiD_Numerics_HessianEE_StdUnSetup_hh
#define COOLFluiD_Numerics_HessianEE_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "HessEEData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Jurek Majevksi
 */
class StdUnSetup : public HessEECom {
public:

  /**
   * Constructor.
   */
  explicit StdUnSetup(std::string name) : HessEECom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdUnSetup()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_StdUnSetup_hh

