#ifndef COOLFluiD_Numerics_RemeshMeandros_WriteMetricCom_hh
#define COOLFluiD_Numerics_RemeshMeandros_WriteMetricCom_hh

//////////////////////////////////////////////////////////////////////////////

#include "RMeshMeData.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/Node.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// sent to Domain to be executed in order to ComputeSpaceResidual the MeshData.
/// @author Jurek Majewski
class WriteMetricCom : public RMeshMeCom
{
public:

  /// Constructor.
  explicit WriteMetricCom(const std::string& name);

  /// Destructor.
  ~WriteMetricCom();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // helper methods

  void  WriteControlSpc( const std::string& fname);

  void  writeTEC( const std::string& fname);

private: // data

  // the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  // the socket to the data handle of the global metric
  Framework::DataSocketSink< RealMatrix> socket_glob_metric;

  /// reference length in mesh
  CFreal _refH;

}; // class WriteMetricCom

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RemeshMeandros_WriteMetricCom_hh

