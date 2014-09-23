#ifndef COOLFluiD_Numerics_MeshLaplacianSmoothing_LaplacianSmoothingData_hh
#define COOLFluiD_Numerics_MeshLaplacianSmoothing_LaplacianSmoothingData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "Framework/MeshAdapterData.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/FaceTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * MeshLaplacianSmoothingCom 's that compose the MeshLaplacianSmoothing.
   *
   * @see LaplacianSmoothingCom
   *
   * @author Thomas Wuilbaut
   */
class LaplacianSmoothingData : public Framework::MeshAdapterData {
public:

  /**
   * Default constructor without arguments
   */
  LaplacianSmoothingData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~LaplacianSmoothingData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Sets up the data
   */
   void setup();

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "LaplacianSmoothing";
  }

  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> >
  getGeoWithNodesBuilder()
  {
    return &_geoWithNodesBuilder;
  }

  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> >
  getStdTrsGeoBuilder()
  {
    return &_stdTrsGeoBuilder;
  }

  /**
   * @return the GeometricEntity builder for faces in FVMCC
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
  getFaceTrsGeoBuilder()
  {
    return &_faceTrsGeoBuilder;
  }

private:

  // builder of GeometricEntity's with Node's
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> _geoWithNodesBuilder;

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;

  /// builder for face TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceTrsGeoBuilder;

}; // end of class LaplacianSmoothingData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for MeshLaplacianSmoothing
typedef Framework::MethodCommand<LaplacianSmoothingData> LaplacianSmoothingCom;

/// Definition of a command provider for MeshLaplacianSmoothing
typedef Framework::MethodCommand<LaplacianSmoothingData>::PROVIDER LaplacianSmoothingComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshLaplacianSmoothing_LaplacianSmoothingData_hh
