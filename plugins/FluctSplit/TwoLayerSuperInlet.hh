#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerSuperInlet_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerSuperInlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a supersonic inlet command
 *
 * @author Thomas Wuilbaut
 *
 */

class TwoLayerSuperInlet : public SuperInlet {
public:

  /**
   * Constructor
   */
  TwoLayerSuperInlet(const std::string& name);

  /**
   * Default destructor
   */
  ~TwoLayerSuperInlet();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on the current TRS
   */
  void executeOnTrs();

protected: // data

  /// socket for intermediate Rhs
  Framework::DataSocketSink<
                            CFreal> socket_interRhs;

  /// socket for intermediate State's
  Framework::DataSocketSink<
                            Framework::State*> socket_interStates;

}; // end of class TwoLayerSuperInlet

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerSuperInlet_hh
