#ifndef COOLFluiD_Numerics_MeshFEMMove_BasePrepare_hh
#define COOLFluiD_Numerics_MeshFEMMove_BasePrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "FEMMoveData.hh"
#include "Common/SafePtr.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */

class BasePrepare : public FEMMoveCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  BasePrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~BasePrepare()
  {
  }

  /**
   * Set up member data
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected: // functions

  /**
   * Move the boundary nodes
   */
  virtual void moveBoundaries()
  {
  }

protected: // data

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the node's
//  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;

  /// Temporary Vector
  RealVector _tmpVector;

}; // class BasePrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshFEMMove_BasePrepare_hh

