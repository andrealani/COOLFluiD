#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_CoupledPrepare_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_CoupledPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpringAnalogyData.hh"
#include "MeshTools/QualityCalculator.hh"
#include "Common/SafePtr.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */

class CoupledPrepare : public SpringAnalogyCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  CoupledPrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~CoupledPrepare()
  {
  }

 /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // functions

  /**
   * For each node, compute the quality of the worse cell connected to it
   */
//   void computeWorstQualityCells();

  /**
   * Move the boundary nodes
   */
  void moveBoundaries();

  /**
   * For each node, compute the distance to the closest boundary
   */
  void computeDistanceToBoundary();

private: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// storage of displacements vector
  Framework::DataSocketSink<RealVector> socket_nodalDisplacements;

  /// name of the interface
  std::string _interfaceName;

  /// type of data received
  std::string _dataTypeReceived;

  /// Functor that computes the requested norm
  Common::SelfRegistPtr<MeshTools::QualityCalculator> _computeQuality;

}; // class CoupledPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshAdapterSpringAnalogy_CoupledPrepare_hh

