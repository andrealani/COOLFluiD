#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerStrongMirrorEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerStrongMirrorEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "StrongMirrorEuler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for space time two-layer schemes
 *
 * @author Thomas Wuilbaut
 *
 */
class TwoLayerStrongMirrorEuler2DCons : public StrongMirrorEuler2DCons {

public:

  /**
   * Constructor.
   */
  TwoLayerStrongMirrorEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~TwoLayerStrongMirrorEuler2DCons();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

private:

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_interRhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink<Framework::State*> socket_interStates;

}; // end of class TwoLayerStrongMirrorEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerStrongMirrorEuler2DCons_hh
