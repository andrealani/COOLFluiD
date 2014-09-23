#include "FluctSplit/Splitter.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Splitter::Splitter(const std::string& name) :
  Framework::MethodStrategy<FluctuationSplitData>(name),
  m_dim(0),
  _nbEquations(0),
  _firstVarID(0),
  _lastVarID(0),
  _nbStatesInCell(0),
  _blockSeparator(0),
  m_normals(CFNULL),
  _adimNormal(),
  socket_isBState("isBState"),
  socket_normals("normals"),
  socket_volumes("volumes"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

Splitter::~Splitter()
{
}

//////////////////////////////////////////////////////////////////////////////

void Splitter::setup()
{
  MethodStrategy<FluctuationSplitData>::setup();
  
  m_dim = PhysicalModelStack::getActive()->getDim();
  _adimNormal.resize(PhysicalModelStack::getActive()->getDim());

  _blockSeparator = getMethodData().getBlockSeparator();
  setBlockData();
}

//////////////////////////////////////////////////////////////////////////////

void Splitter::unsetup()
{
  // unsetup actions should go here

  MethodStrategy<FluctuationSplitData>::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
Splitter::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  
  result.push_back(&socket_isBState);
  result.push_back(&socket_normals);
  result.push_back(&socket_volumes);
  result.push_back(&socket_states);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Splitter::configure ( Config::ConfigArgs& args )
{
  Framework::MethodStrategy<FluctuationSplitData>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

