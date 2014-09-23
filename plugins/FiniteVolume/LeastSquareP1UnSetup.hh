#ifndef COOLFluiD_Numerics_FiniteVolume_LeastSquareP1UnSetup_hh
#define COOLFluiD_Numerics_FiniteVolume_LeastSquareP1UnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Andrea Lani
 */
class LeastSquareP1UnSetup : public StdUnSetup {
public:

  /**
   * Constructor.
   */
  explicit LeastSquareP1UnSetup(const std::string& name) :
    StdUnSetup(name)
  {
  }

  /**
   * Destructor.
   */
  ~LeastSquareP1UnSetup()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LeastSquareP1UnSetup_hh

