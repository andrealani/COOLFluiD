#include "MeshFEMMove/MeshFEMMove.hh"
#include "BasePrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/SelfRegistPtr.hh"
#include "MeshTools/ComputeShortestDistance.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BasePrepare, FEMMoveData, MeshFEMMoveModule> BasePrepareProvider("BasePrepare");

//////////////////////////////////////////////////////////////////////////////

void BasePrepare::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

BasePrepare::BasePrepare(const std::string& name) :
FEMMoveCom(name),
  socket_nodes("nodes"),
//  socket_pastNodes("pastNodes"),
  _tmpVector()
{
   addConfigOptionsTo(this);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
BasePrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
//  result.push_back(&socket_pastNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BasePrepare::setup()
{
  FEMMoveCom::setup();

  _tmpVector.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void BasePrepare::execute()
{

  moveBoundaries();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD
