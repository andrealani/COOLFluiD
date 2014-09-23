#ifndef COOLFluiD_MarcoTest_StdUnSetup_hh
#define COOLFluiD_MarcoTest_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"
#include "MarcoTest/MarcoTestMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

/**
 * This command deletes data which needs to "manually" be deleted 
 * (e.g. arrays of pointers to objects)
 * 
 * @author Andrea Lani
 */
class StdUnSetup : public MarcoTestMethodCom {
public:

  /**
   * Constructor.
   */
  explicit StdUnSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdUnSetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  virtual void execute();

protected:
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// socket for proxy to be able to use the data handle of nodal states uniformly
  /// independently from the actual storage type being RealVector or State*
  Framework::DataSocketSink<
                            Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;
  
}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MarcoTest

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MarcoTest_StdUnSetup_hh
