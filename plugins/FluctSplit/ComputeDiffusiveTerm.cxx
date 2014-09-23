#include "Framework/MeshData.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("StorePeCell", "Stores the cell Peclet number");
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeDiffusiveTerm::ComputeDiffusiveTerm(const std::string& name) :
  Framework::MethodStrategy<FluctuationSplitData>(name),
  socket_updateCoeff("updateCoeff"),
  socket_normals("normals"),
  socket_volumes("volumes"),
  socket_Pe_cell("Pe_cell")
{
  addConfigOptionsTo(this);
  
  _store_Pe_cell = false;
  this->setParameter("StorePeCell", &_store_Pe_cell);
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeDiffusiveTerm::~ComputeDiffusiveTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveTerm::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveTerm::setMeshData()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeDiffusiveTerm::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_Pe_cell);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveTerm::setup()
{
  const CFuint nbcells =
    MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();
  
  DataHandle<CFreal> Pe_cell = socket_Pe_cell.getDataHandle();
  Pe_cell.resize( nbcells );
}
      
//////////////////////////////////////////////////////////////////////////////
  
} // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
