#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerUnSetup_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Thomas Wuilbaut
 */
class TwoLayerUnSetup : public StdUnSetup {
public:

  /**
   * Constructor.
   */
  explicit TwoLayerUnSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~TwoLayerUnSetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

private:

  /// socket for intermediate normals
  Framework::DataSocketSink<
                            InwardNormalsData*> socket_interNormals;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerUnSetup_hh

