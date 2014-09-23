#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "StdALEUnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdALEUnSetup, FluctuationSplitData, FluctSplitSpaceTimeModule> StdALEUnSetupProvider("StdALEUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdALEUnSetup::StdALEUnSetup(const std::string& name) : FluctuationSplitCom(name),
  socket_pastNormals("pastNormals"),
  socket_pastNodes("pastNodes")
{
  cout << "StdALEUnSetup" << endl;
}

//////////////////////////////////////////////////////////////////////////////

StdALEUnSetup::~StdALEUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdALEUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_pastNormals);
  result.push_back(&socket_pastNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdALEUnSetup::execute()
{
  CFAUTOTRACE;

  deleteAllPtr(socket_pastNormals);
  deleteAllPtr(socket_pastNodes);
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluctSplit



} // namespace COOLFluiD
