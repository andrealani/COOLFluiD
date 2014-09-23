#ifndef COOLFluiD_SpectralFV_StdExtrapolate_hh
#define COOLFluiD_SpectralFV_StdExtrapolate_hh

#include "SpectralFV/SpectralFVMethodData.hh"
#include "SpectralFV/SpectralFVElementData.hh"
#include "SpectralFV/LineSpectralFVElementData.hh"
#include "SpectralFV/TriagSpectralFVElementData.hh"

namespace COOLFluiD {
  namespace SpectralFV {


/// This is a standard command to setup the empty method
class StdExtrapolate : public SpectralFVMethodCom {

public:

  /// Constructor
  explicit StdExtrapolate(const std::string& name);

  /// Destructor
  ~StdExtrapolate();

  /// Execute processing actions
  void execute();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected:

  /// storage for interpolated values in the mesh vertices
  Framework::DataSocketSink<RealVector> socket_nstates;

  /**
   * socket with proxy to be able to use the data handle of nodal states
   * uniformly independently from the actual storage type being RealVector or
   * State*
   */
  Framework::DataSocketSink<Framework::ProxyDofIterator< RealVector >* > socket_nstatesProxy;

  /// socket for state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

};  // class StdExtrapolate


  }  // namespace SpectralFV
}  // namespace COOLFluiD

#endif // COOLFluiD_SpectralFV_StdExtrapolate_hh

