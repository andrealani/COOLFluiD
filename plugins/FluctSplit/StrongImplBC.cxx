#include "FluctSplit/StrongImplBC.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

StrongImplBC::StrongImplBC(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_discardTimeJacob("discardTimeJacob"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  _useForInitialization(false)
{
}

//////////////////////////////////////////////////////////////////////////////

StrongImplBC::~StrongImplBC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongImplBC::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_discardTimeJacob);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongImplBC::setup()
{
  CFAUTOTRACE;
  
  FluctuationSplitCom::setup();
}
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
