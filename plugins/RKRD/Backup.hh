#ifndef COOLFluiD_Numerics_RKRD_Backup_hh
#define COOLFluiD_Numerics_RKRD_Backup_hh

//////////////////////////////////////////////////////////////////////////////

#include "RKRD/RKRDData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Tiago Quintino
class Backup : public RKRDCom {
public:

  /// Constructor.
  explicit Backup(const std::string& name);

  /// Destructor.
  virtual ~Backup();

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for states at different k-steps
  Framework::DataSocketSink<RealMatrix> socket_kstates;

  /// socket to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RKRD_Backup_hh

