#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_StdUnSetup_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "SpringAnalogyData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Thomas Wuilbaut
 */
class StdUnSetup : public SpringAnalogyCom {
public:

  /**
   * Constructor.
   */
  explicit StdUnSetup(std::string name);

  /**
   * Destructor.
   */
  ~StdUnSetup();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshAdapterSpringAnalogy_StdUnSetup_hh

