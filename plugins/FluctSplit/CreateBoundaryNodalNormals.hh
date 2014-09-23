#ifndef COOLFluiD_Numerics_FluctSplit_CreateBoundaryNodalNormals_hh
#define COOLFluiD_Numerics_FluctSplit_CreateBoundaryNodalNormals_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/Trio.hh"

#include "MathTools/RealVector.hh"

#include "Framework/Storage.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/DataSocketSink.hh"

#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {  class TopologicalRegionSet;  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class creates the normals for strong BC treatment
/// @author Andrea Lani
class FluctSplit_API CreateBoundaryNodalNormals {
public:

  /// Constructor
  CreateBoundaryNodalNormals
  (Common::SafePtr<Framework::GeometricEntityPool
  <Framework::StdTrsGeoBuilder> > geoBuilder);

  /// Default destructor
  ~CreateBoundaryNodalNormals();

  /// Set the data sockets
  void setDataSockets(Framework::DataSocketSink<InwardNormalsData* > normalsSocket,
                      Framework::DataSocketSink< Common::Trio<CFuint, CFuint,
                        Common::SafePtr<Framework::TopologicalRegionSet> > > faceNeighCellSocket );

  /// Create the bc normals on TopologicalRegionSet's
  void create(const std::vector<Common::SafePtr<Framework::TopologicalRegionSet> >& trsList,
      std::vector< std::vector<RealVector> >& bcNormals);

private: // helper method

  /// Create the bc normals on TopologicalRegionSet's
  void createOnTrs(Common::SafePtr<Framework::TopologicalRegionSet> trs,
                   std::vector<RealVector>& bcNormals);

private:

  /// pointer to the GeometricEntity builder
  Common::SafePtr<Framework::GeometricEntityPool
  <Framework::StdTrsGeoBuilder> > _geoBuilder;

  /// socket for normals
  Framework::DataSocketSink<InwardNormalsData* > socket_normals;

  /// socket for face neighbour cell
  Framework::DataSocketSink<
                              Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
                      socket_faceNeighCell;


}; // class CreateBoundaryNodalNormals

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_CreateBoundaryNodalNormals_hh
