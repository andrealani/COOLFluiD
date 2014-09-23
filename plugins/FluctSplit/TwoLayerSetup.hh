#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerSetup_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action
 * to be executed in order to setup the FluctSplit Method
 * for the space-time two layer variant.
 *
 * @author Thomas Wuilbaut
 *
 */
class TwoLayerSetup : public StdSetup {
public:

  /**
   * Constructor.
   */
  explicit TwoLayerSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~TwoLayerSetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Execute Processing actions
   */
  void execute();

private:

  /// socket for intermediate normals
  Framework::DataSocketSource<
                              InwardNormalsData*> socket_interNormals;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerSetup_hh

