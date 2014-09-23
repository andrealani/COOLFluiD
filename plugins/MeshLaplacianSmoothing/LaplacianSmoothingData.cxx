#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModelImpl.hh"

#include "MeshLaplacianSmoothing/MeshLaplacianSmoothing.hh"
#include "MeshLaplacianSmoothing/LaplacianSmoothingData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<LaplacianSmoothingData>, LaplacianSmoothingData, MeshLaplacianSmoothingModule> nullLaplacianSmoothingComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

LaplacianSmoothingData::LaplacianSmoothingData(Common::SafePtr<Framework::Method> owner)
  : MeshAdapterData(owner),
    _geoWithNodesBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

LaplacianSmoothingData::~LaplacianSmoothingData()
{
}

//////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothingData::configure ( Config::ConfigArgs& args )
{
  MeshAdapterData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothingData::setup()
{
  CFAUTOTRACE;

  _geoWithNodesBuilder.setup();
  _stdTrsGeoBuilder.setup();
  _faceTrsGeoBuilder.setup();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

