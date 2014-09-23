#ifndef COOLFluiD_Numerics_FiniteElement_StdPrepare_hh
#define COOLFluiD_Numerics_FiniteElement_StdPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FiniteElement/FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action tto prepare
   * the computation by clearing the isUpdated data.
   */
class StdPrepare : public FiniteElementMethodCom {
public:

  /**
   * Constructor.
   */
  explicit StdPrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~StdPrepare()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for isUpdated
  Framework::DataSocketSink<
                              bool> socket_isUpdated;

}; // class Prepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_StdPrepare_hh

