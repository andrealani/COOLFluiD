#ifndef COOLFluiD_Numerics_FiniteVolume_StdALEPrepare_hh
#define COOLFluiD_Numerics_FiniteVolume_StdALEPrepare_hh

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
class StdALEPrepare : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit StdALEPrepare(const std::string& name) : 
    CellCenterFVMCom(name),
    socket_nodes("nodes"),
    socket_pastNodes("pastNodes")
  {
    std::cout << "StdALEPrepare" << std::endl;
  
  }

  /**
   * Destructor.
   */
  ~StdALEPrepare()
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
  
  // handle to nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  // handle to past nodes
  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StdALEPrepare_hh

