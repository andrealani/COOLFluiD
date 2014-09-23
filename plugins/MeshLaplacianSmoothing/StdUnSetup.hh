#ifndef COOLFluiD_Numerics_MeshLaplacianSmoothing_StdUnSetup_hh
#define COOLFluiD_Numerics_MeshLaplacianSmoothing_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "LaplacianSmoothingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Thomas Wuilbaut
 */
class StdUnSetup : public LaplacianSmoothingCom {
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

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshLaplacianSmoothing_StdUnSetup_hh

