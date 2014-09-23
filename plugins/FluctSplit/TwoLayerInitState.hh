#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerInitState_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerInitState_hh

//////////////////////////////////////////////////////////////////////////////

#include "InitState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a initalizing solution command
 *
 * @author Thomas Wuilbaut
 *
 */
class TwoLayerInitState : public InitState {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::VarSetTransformer VectorTransformer;

  /**
   * Constructor.
   */
  explicit TwoLayerInitState(const std::string& name);

  /**
   * Destructor.
   */
  ~TwoLayerInitState();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute Processing actions
   */
  void executeOnTrs();


protected: // data

  // the socket to the data handle of the state's
  Framework::DataSocketSink<Framework::State*> socket_interStates;

}; // class TwoLayerInitState

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerInitState_hh

