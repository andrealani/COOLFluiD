#ifndef COOLFluiD_Numerics_FluctSplit_StdUnSetup_hh
#define COOLFluiD_Numerics_FluctSplit_StdUnSetup_hh

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
/// to be executed in order to deallocate data specific to the
/// method to which it belongs.
/// @author Tiago Quintino
class FluctSplit_API StdUnSetup : public FluctuationSplitCom {
public: // member functions

  /// Constructor.
  explicit StdUnSetup(const std::string& name);

  /// Destructor.
  ~StdUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /// Execute Processing actions
  virtual void execute();

protected: // data

  /// socket for inward normals
  Framework::DataSocketSink<
                            InwardNormalsData*> socket_normals;

  /// socket for the proxy of states
  Framework::DataSocketSink
    <Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdUnSetup_hh

