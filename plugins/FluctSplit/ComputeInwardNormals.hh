#ifndef COOLFluiD_Numerics_FluctSplit_ComputeInwardNormals_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeInwardNormals_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeNormals.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/DataSocketSink.hh"

#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// Base class for the Inward Normals Computers
/// @author Thomas Wuilbaut
class FluctSplit_API ComputeInwardNormals : public Framework::ComputeNormals {

public: // functions

  /// Constructor
  ComputeInwardNormals();

  /// Destructor
  ~ComputeInwardNormals();

  /// setup and return the GeometricEntity Builder
  /// @return the GeometricEntity builder
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> >
  getStdTrsGeoBuilder()
  {
    _stdTrsGeoBuilder.setup();
    return &_stdTrsGeoBuilder;
  }

  /// Set the DataSocket's that corresponds to the normals
  /// @param normals SafePtr to the normals datasocket
  void setNormalsSockets(Common::SafePtr<Framework::DataSocketSink< InwardNormalsData*> > normals);

  /// Set the DataSocket's that corresponds to the normalsData
  /// @param normalsData SafePtr to the normals datasocket
  void setNormalsDataSockets(Common::SafePtr<Framework::DataSocketSink< CFreal> > normalsData);

  /// Set the DataSocket's that corresponds to the tempSize
  /// @param tempSize SafePtr to the normals datasocket
  void setTempSizeSockets(Common::SafePtr<Framework::DataSocketSink< CFuint> > tempSize);

  /// Set the DataSocket's that corresponds to the past Normals
  /// @param pastNormals SafePtr to the normals datasocket
  void setPastNormalsSockets(Common::SafePtr<Framework::DataSocketSink<InwardNormalsData*> > pastNormals)
  {
    socket_pastNormals = pastNormals;
  }

  /// Set the DataSocket's that corresponds to the intermediate Normals
  /// @param interNormals SafePtr to the normals datasocket
  void setInterNormalsSockets(Common::SafePtr<Framework::DataSocketSink<InwardNormalsData*> > interNormals)
  {
    socket_interNormals = interNormals;
  }

protected: // data

  /// socket for the inward normals storage
  Common::SafePtr< Framework::DataSocketSink<InwardNormalsData*> > socket_normals;

  /// @todo missing documentation
  Common::SafePtr< Framework::DataSocketSink<CFreal> > socket_normalsData;

  /// @todo missing documentation
  Common::SafePtr< Framework::DataSocketSink< CFuint> > socket_tempSize;

  /// socket for the inward normals storage
  Common::SafePtr< Framework::DataSocketSink<InwardNormalsData*> > socket_pastNormals;

  /// socket for the inward normals storage
  Common::SafePtr< Framework::DataSocketSink<InwardNormalsData*> > socket_interNormals;

private: // data

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;

}; // end of class ComputeInwardNormalsHexaP1

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeInwardNormals_hh
