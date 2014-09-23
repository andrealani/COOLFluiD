#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "TwoLayerUnSetup.hh"
#include "Framework/MeshData.hh"
#include "InwardNormalsData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerUnSetup, FluctuationSplitData, FluctSplitSpaceTimeModule> TwoLayerUnSetupProvider("TwoLayerUnSetup");


//////////////////////////////////////////////////////////////////////////////

TwoLayerUnSetup::TwoLayerUnSetup(const std::string& name) :
  StdUnSetup(name),
  socket_interNormals("interNormals")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerUnSetup::~TwoLayerUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StdUnSetup::needsSockets();

  result.push_back(&socket_interNormals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerUnSetup::execute()
{
  StdUnSetup::execute();

  if (SubSystemStatusStack::getActive()->isMovingMesh()){
    deleteAllPtr(socket_interNormals);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
