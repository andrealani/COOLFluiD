#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/QRadSetup.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<QRadSetup, CellCenterFVMData, FiniteVolumeModule> 
qRadSetupProvider("QRadSetup");

//////////////////////////////////////////////////////////////////////////////

QRadSetup::QRadSetup(const std::string& name) :
  CellCenterFVMCom(name),
  socket_qrad("qrad"),
  socket_states("states")

{
}

//////////////////////////////////////////////////////////////////////////////

QRadSetup::~QRadSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > QRadSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_qrad);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void QRadSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle <CFreal> qrad = socket_qrad.getDataHandle();
  qrad.resize(states.size());
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > QRadSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
