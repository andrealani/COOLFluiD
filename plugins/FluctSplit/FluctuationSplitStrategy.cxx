#include "Framework/MeshData.hh"
#include "Framework/JacobianLinearizer.hh"

#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

FluctuationSplitStrategy::FluctuationSplitStrategy(const std::string& name) :
  Framework::MethodStrategy<FluctuationSplitData>(name),
  socket_normals("normals"),
  socket_isBState("isBState"),
  socket_volumes("volumes"),
  _cells(CFNULL),
  m_linearStates(CFNULL),
  _stdTrsGeoBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

FluctuationSplitStrategy::~FluctuationSplitStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitStrategy::configure ( Config::ConfigArgs& args )
{
  Framework::MethodStrategy<FluctuationSplitData>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitStrategy::unsetup()
{  
  CFAUTOTRACE;

  // unsetup the geometric entity builder
  _stdTrsGeoBuilder.unsetup();

  // resset a pointer to the cells
  _cells.reset( CFNULL );

  Framework::MethodStrategy<FluctuationSplitData>::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitStrategy::setup()
{
  CFAUTOTRACE;
  
  Framework::MethodStrategy<FluctuationSplitData>::setup();

  // set a pointer to the cells
  _cells.reset(MeshDataStack::getActive()->getTrs("InnerCells"));

  // set up the geometric entity builder
  _stdTrsGeoBuilder.setup();

}

//////////////////////////////////////////////////////////////////////////////

void FluctuationSplitStrategy::prepare()
{
  // default do nothing
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
FluctuationSplitStrategy::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;

   result.push_back(&socket_normals);
   result.push_back(&socket_isBState);
   result.push_back(&socket_volumes);

   return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
