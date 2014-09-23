#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "TwoLayerSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "InwardNormalsData.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerSetup, FluctuationSplitData, FluctSplitSpaceTimeModule> twoLayerSetupProvider("TwoLayerSetup");

//////////////////////////////////////////////////////////////////////////////

TwoLayerSetup::TwoLayerSetup(const std::string& name) : StdSetup(name),
  socket_interNormals("interNormals")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerSetup::~TwoLayerSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
TwoLayerSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
    StdSetup::providesSockets();

  result.push_back(&socket_interNormals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSetup::execute()
{
  CFAUTOTRACE;

  StdSetup::execute();

  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();

  /// @todo try to find a way to avoid this because it takes lots of memory
  // If moving mesh, then allocate data to store the past Normals
  if (SubSystemStatusStack::getActive()->isMovingMesh()) {

    DataHandle<InwardNormalsData*> interNormals = socket_interNormals.getDataHandle();

    interNormals.resize(normals.size());

    for (CFuint i=0; i < normals.size(); ++i){
      interNormals[i] = new InwardNormalsData(*normals[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
