#ifndef COOLFluiD_UFEM_SetupPLaS_hh
#define COOLFluiD_UFEM_SetupPLaS_hh

#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "UFEM/UFEMSolverData.hh"

namespace COOLFluiD {
  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/// This class is a NumericalCommand to setup PLaS required data
class UFEM_API SetupPLaS : public UFEMSolverCom {

 public:

  /// Constructor
  explicit SetupPLaS(const std::string& name);

  /// Destructor
  ~SetupPLaS();

  /// Execute
  void execute();

  /// Returns the DataSocketSources's this command provides
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets() {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > r;
    r.push_back(&s_nvfraction);
    r.push_back(&s_nvolume);
    r.push_back(&s_evolume);
    r.push_back(&s_ienormals);
    r.push_back(&s_benormals);
    r.push_back(&s_faceNeighCell);
    return r;
  }

  /// Returns the DataSocketSink's this command needs
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets() {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_nodes);
    return r;
  }

private:

  /// Set the mapping between each boundary Face and its neighbour cell
  void setFaceNeighCell();

protected:

  /// Socket to access nodes
  Framework::DataSocketSink< Framework::Node*,Framework::GLOBAL > s_nodes;

  /// Socket to provide node-wise void fraction
  Framework::DataSocketSource< CFreal > s_nvfraction;

  /// Socket to provide node-wise volume
  Framework::DataSocketSource< CFreal > s_nvolume;

  /// Socket to provide element-wise volume
  Framework::DataSocketSource< CFreal > s_evolume;

  /// Socket to provide 'inner' elements normals
  Framework::DataSocketSource< std::vector< RealVector > > s_ienormals;

  /// Socket to provide 'boundary' elements faces normals
  Framework::DataSocketSource< std::vector< RealVector > > s_benormals;

  /// socket for GeometricEntity mapper to corresponding TopologicalRegionSet
  Framework::DataSocketSource<std::pair<CFuint,CFuint> > s_faceNeighCell;

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace UFEM
}  // namespace COOLFluiD

#endif // COOLFluiD_UFEM_SetupPLaS_hh

