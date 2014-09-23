#ifndef COOLFluiD_Numerics_FluctSplit_UpdateCoeffUnSetup_hh
#define COOLFluiD_Numerics_FluctSplit_UpdateCoeffUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "Framework/ProxyDofIterator.hh"
#include "FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to deallocate theta socket
/// To be used in conjunction with ThetaSetup
/// @author Nadege Villedieu
/// @author Tiago Quintino
class FluctSplit_API UpdateCoeffUnSetup : public FluctuationSplitCom {
public: // member functions

  /// Constructor.
  explicit UpdateCoeffUnSetup(const std::string& name);

  /// Destructor.
  virtual ~UpdateCoeffUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /// Executes this command
  virtual void execute();

protected: // data

  /// socket for theta blending coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdUnSetup_hh

