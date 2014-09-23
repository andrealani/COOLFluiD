#ifndef COOLFluiD_Numerics_FiniteVolume_BDF2ALEUpdate_hh
#define COOLFluiD_Numerics_FiniteVolume_BDF2ALEUpdate_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdALEUpdate.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command to be executed after
   * the mesh has been updated
   */

//////////////////////////////////////////////////////////////////////////////

class BDF2ALEUpdate : public StdALEUpdate {
public:

  /**
   * Constructor.
   */
  explicit BDF2ALEUpdate(const std::string& name) : 
    StdALEUpdate(name),
    socket_avNormals("avNormals")
  {
  }
  
  /**
   * Destructor.
   */
  ~BDF2ALEUpdate()
  {
  }

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:

  /**
   * Update the normals (at intermediate time)
   */
  void updateNormalsData();

protected: // data
  
  // average normals
  Framework::DataSocketSink< CFreal> socket_avNormals;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_BDF2ALEUpdate_hh

