#ifndef COOLFluiD_Numerics_RKRD_Update_hh
#define COOLFluiD_Numerics_RKRD_Update_hh

//////////////////////////////////////////////////////////////////////////////

#include "RKRD/RKRDData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

/// Updates the states after the RKRD steps
/// @author Tiago Quintino
class Update : public RKRDCom {
public:

  /// Constructor.
  explicit Update(const std::string& name);

  /// Destructor.
  virtual ~Update();

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

#endif // COOLFluiD_Numerics_RKRD_Update_hh

