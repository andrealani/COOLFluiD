#include "MeshAdapterSpringAnalogy/MeshAdapterSpringAnalogy.hh"
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

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BasePrepare, SpringAnalogyData, MeshAdapterSpringAnalogyModule> BasePrepareProvider("BasePrepare");

//////////////////////////////////////////////////////////////////////////////

void BasePrepare::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::string >("DataType","Type of Data Received Coordinates or Displacements.");
}

//////////////////////////////////////////////////////////////////////////////

BasePrepare::BasePrepare(const std::string& name) :
SpringAnalogyCom(name),
  socket_nodes("nodes"),
  socket_pastNodes("pastNodes"),
  socket_nodalDisplacements("nodalDisplacements"),
  _tmpVector()
{
   addConfigOptionsTo(this);

  _interfaceName = "";
   setParameter("Interface",&_interfaceName);

  _dataTypeReceived = "";
   setParameter("DataType",&_dataTypeReceived);

}

//////////////////////////////////////////////////////////////////////////////

void BasePrepare::setup()
{
  SpringAnalogyCom::setup();

  _tmpVector.resize(PhysicalModelStack::getActive()->getDim());
}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
BasePrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_nodalDisplacements);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BasePrepare::execute()
{
  CFAUTOTRACE;

  moveBoundaries();

 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD
