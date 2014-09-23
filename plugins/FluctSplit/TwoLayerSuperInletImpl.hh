#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerSuperInletImpl_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerSuperInletImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/SuperInletImpl.hh"

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
class TwoLayerSuperInletImpl : public SuperInletImpl {
public:

  /**
   * Constructor
   */
  TwoLayerSuperInletImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~TwoLayerSuperInletImpl();

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

  // the socket to the data handle of the intermediate rhs
  Framework::DataSocketSink<CFreal> socket_interRhs;

  // the socket to the data handle of the intermediate state's
  Framework::DataSocketSink<Framework::State*> socket_interStates;

}; // end of class TwoLayerSuperInletImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerSuperInletImpl_hh
