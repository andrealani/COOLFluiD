#ifndef COOLFluiD_UFEM_StdClean_hh
#define COOLFluiD_UFEM_StdClean_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMSolverData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to setup the MeshData.
/// @author Tamas Banyai
class UFEM_API StdClean : public UFEMSolverCom {
public:

  /// Constructor.
  explicit StdClean(const std::string& name);

  /// Destructor.
  virtual ~StdClean();

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket with flags to check if a state's equation has been updated or not
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// socket with flags to check if a state's equation has been updated or not
  Framework::DataSocketSink<CFreal> socket_rhs;

}; // class Clean

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_StdClean_hh

