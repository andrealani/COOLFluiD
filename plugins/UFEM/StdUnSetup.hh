#ifndef COOLFluiD_UFEM_StdUnSetup_hh
#define COOLFluiD_UFEM_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMSolverData.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Tiago Quintino
class UFEM_API StdUnSetup : public UFEMSolverCom {
public:

  /// Constructor.
  explicit StdUnSetup(const std::string& name);

  /// Destructor.
  ~StdUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

protected:

  /// Socket to store the proxy of states
  Framework::DataSocketSink<Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;


}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_StdUnSetup_hh

