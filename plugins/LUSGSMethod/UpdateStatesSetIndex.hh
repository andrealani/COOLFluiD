#ifndef COOLFluiD_Numerics_LUSGSMethod_UpdateStatesSetIndex_hh
#define COOLFluiD_Numerics_LUSGSMethod_UpdateStatesSetIndex_hh

//////////////////////////////////////////////////////////////////////////////

#include "LUSGSMethod/LUSGSIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class updates the index of the states set which will be updated.
   * @author Matteo Parsani
   * @author Kris Van den Abeele
   */
class UpdateStatesSetIndex : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit UpdateStatesSetIndex(std::string name);

  /**
   * Destructor.
   */
  ~UpdateStatesSetIndex() {}

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks.
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

}; // class UpdateStatesSetIndex

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_UpdateStatesSetIndex_hh

