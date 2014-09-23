#ifndef COOLFluiD_Numerics_FiniteVolume_BDF2ALEPrepare_hh
#define COOLFluiD_Numerics_FiniteVolume_BDF2ALEPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents
 * @todo missing documentation
 */

class BDF2ALEPrepare : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit BDF2ALEPrepare(const std::string& name) : 
    CellCenterFVMCom(name),
    socket_nodes("nodes"),
    socket_pastNodes("pastNodes"),
    socket_pastPastNodes("pastPastNodes")
  {
    std::cout << "BDF2ALEPrepare" << std::endl;
  }
  
  /**
   * Destructor.
   */
  ~BDF2ALEPrepare()
  {
  }

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:
  
  // mesh nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  // mesh past nodes
  Framework::DataSocketSink< Framework::Node*> socket_pastNodes;
  
  // mesh past past nodes
  Framework::DataSocketSink< Framework::Node*> socket_pastPastNodes;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_BDF2ALEPrepare_hh

