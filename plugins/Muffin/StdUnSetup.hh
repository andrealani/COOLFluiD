#ifndef COOLFluiD_Muffin_StdUnSetup_hh
#define COOLFluiD_Muffin_StdUnSetup_hh

#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {


/// This is a standard command to deallocate data specific to Muffin method
class StdUnSetup : public MuffinCom {

public:

  /// Constructor
  StdUnSetup(const std::string& name) :
      MuffinCom(name),
      s_nstatesproxy("nstatesProxy")
  {}

  /// Destructor
  ~StdUnSetup()
  {}

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_nstatesproxy);
    return r;
  }

  /// Execute processing actions
  void execute();


private:  // data

  /// Socket to store the proxy of states
  Framework::DataSocketSink< Framework::ProxyDofIterator< RealVector >* > s_nstatesproxy;

}; // class StdUnSetup


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_StdUnSetup_hh

