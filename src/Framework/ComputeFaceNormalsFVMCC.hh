#ifndef COOLFluiD_Framework_ComputeFaceNormalsFVMCC_hh
#define COOLFluiD_Framework_ComputeFaceNormalsFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeNormals.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/**
 * Base class for the Face Normals Computers
 *
 * @author Andrea Lani
 *
 */
class Framework_API ComputeFaceNormalsFVMCC : public Framework::ComputeNormals {

public:

  /**
   * Constructor
   */
  ComputeFaceNormalsFVMCC() : Framework::ComputeNormals()
  {
  }

  /**
   * Destructor
   */
  virtual ~ComputeFaceNormalsFVMCC()
  {
  }
  
  /**
   * Set some useful sockets
   */
  void setSockets
  (Common::SafePtr<Framework::DataSocketSink<CFreal> > normals,
   Common::SafePtr<Framework::DataSocketSink<CFint> > isOutward)
  {
    socket_normals = normals;
    socket_isOutward = isOutward;
  }
  
protected:
  
  /// socket for normals storage
  Common::SafePtr<Framework::DataSocketSink<CFreal> > socket_normals;
  
  /// socket storing the IDs of the element which the normal is outward
  Common::SafePtr<Framework::DataSocketSink<CFint> > socket_isOutward;
  
}; // end of class ComputeFaceNormalsFVMCCHexaP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeFaceNormalsFVMCC_hh
